#!/usr/bin/env python3
"""
Eneza: Thyroid Carcinoma Histologic Subtype Classifier
=========================================================
End-to-end, reproducible machine-learning pipeline for cancer subtype
classification from tumour gene expression.

Data source (local)
-------------------
- TCGA THCA RNA-seq v2 RSEM: ``data/data_mrna_seq_v2_rsem.txt``
- TCGA THCA sample phenotype: ``data/data_clinical_sample.txt``
- TCGA THCA patient survival: ``data/data_clinical_patient.txt``

Clinical task
-------------
Classify thyroid papillary carcinoma (THCA) histologic subtypes from the
``TUMOR_TYPE`` column of the TCGA clinical sample file:

  * Classical / usual type
  * Follicular variant
  * Tall cell variant
  * Other (dropped because n < 10)


How to run
----------
    python scripts/run.py

All outputs are written to ``output/`` by default.
"""

from __future__ import annotations
from sklearn.svm import LinearSVC
from sklearn.preprocessing import LabelBinarizer, StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.model_selection import (
    StratifiedKFold,
    cross_validate,
    train_test_split,
)
from sklearn.metrics import (
    ConfusionMatrixDisplay,
    auc,
    balanced_accuracy_score,
    classification_report,
    f1_score,
    precision_recall_curve,
    precision_score,
    recall_score,
    roc_auc_score,
    roc_curve,
)
from sklearn.linear_model import LogisticRegression
from sklearn.feature_selection import SelectFdr, SelectKBest, f_classif, mutual_info_classif
from sklearn.ensemble import RandomForestClassifier
from sklearn.base import BaseEstimator, TransformerMixin
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

import json
import logging
import os
import warnings
from dataclasses import asdict, dataclass
from functools import partial
from pathlib import Path
from typing import Any, Dict, Iterable, List, Literal, Sequence, Tuple

# Suppress environment / library noise before importing plotting libraries
warnings.filterwarnings("ignore")

# non-interactive backend; prevents GUI warnings on headless servers
matplotlib.use("Agg")

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

LOG = logging.getLogger("eneza")


@dataclass(frozen=True)
class Config:
    """Immutable configuration container.

    Notes
    -----
    * ``selector`` controls the univariate feature-selection step.
    * ``selected_features`` is used when ``selector == "kbest"``.
    * ``alpha_threshold`` is used when ``selector == "fdr"``.
    * All paths are relative to the working directory by default.
    """

    # Files
    expression_path: str = "./data/data_mrna_seq_v2_rsem.txt"
    clinical_sample_path: str = "./data/data_clinical_sample.txt"
    patient_path: str = "./data/data_clinical_patient.txt"

    # Phenotype
    sample_id_column: str = "SAMPLE_ID"
    label_column: str = "TUMOR_TYPE"
    min_class_size: int = 20  # drop subtypes with fewer samples

    # Pre-processing
    top_variable_genes: int = 5000  # unsupervised variance filter for EDA / speed
    log_transform: bool = True  # RSEM is positive continuous -> log2(x+1)

    # Feature selection inside CV
    selector: Literal["kbest", "fdr"] = "kbest"
    selected_features: int = 100
    alpha_threshold: float = 0.05

    # Model
    classifier: Literal["logistic", "random_forest", "linear_svc"] = "logistic"
    cv_splits: int = 5
    test_size: float = 0.20
    random_state: int = 42
    class_weight: str = "balanced"
    max_iter: int = 5000

    # DESeq/edgeR differential-expression feature selection
    deseq_results_dir: str = "./data/deseq_data"
    deseq_k: int = 100
    deseq_score_col: str = "FDR"  # or "PValue"
    deseq_rank_method: str = "min"  # min FDR across pairwise comparisons

    # Outputs
    output_dir: str = "./output"

    # Runtime
    verbosity: int = logging.INFO

    def model_id(self) -> str:
        return f"{self.classifier}_{self.selector}_k{self.selected_features}"


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------


def setup_logging(level: int = logging.INFO) -> None:
    """Configure root logger for readable console output."""
    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(name)s | %(levelname)s | %(message)s",
        datefmt="%H:%M:%S",
    )


def load_expression_matrix(path: str) -> pd.DataFrame:
    """Load raw RSEM expression matrix (genes x samples)."""
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"Expression file not found: {p.resolve()}")
    LOG.info("Loading expression matrix from %s", p)
    df = pd.read_csv(p, sep="\t")
    LOG.info("Raw expression shape: %s", df.shape)
    if df.shape[0] < 10 or df.shape[1] < 10:
        raise ValueError(
            f"Expression matrix seems empty or malformed: {df.shape}")
    return df


def clean_expression_matrix(df: pd.DataFrame) -> pd.DataFrame:
    """Set Hugo_Symbol index, drop NA/duplicate genes, remove Entrez ID."""
    required = {"Hugo_Symbol", "Entrez_Gene_Id"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Expression matrix missing required columns: {missing}")

    cleaned = (
        df.dropna(subset=["Hugo_Symbol"])
        .drop_duplicates(subset=["Hugo_Symbol"], keep="first")
        .set_index("Hugo_Symbol")
        .drop(columns="Entrez_Gene_Id")
    )
    # Force numeric; non-numeric sample columns become NaN and are dropped
    cleaned = cleaned.apply(pd.to_numeric, errors="coerce")
    cleaned = cleaned.dropna(axis=1, how="any")  # drop samples with any NaN
    cleaned = cleaned.loc[cleaned.var(axis=1) > 0]  # keep variable genes
    LOG.info("Cleaned expression shape: %s", cleaned.shape)
    return cleaned


def log2_transform(matrix: pd.DataFrame) -> pd.DataFrame:
    """Apply log2(x + 1) to stabilise variance and reduce skewness."""
    return np.log2(matrix + 1.0)


def select_top_variable_genes(matrix: pd.DataFrame, n_genes: int) -> pd.DataFrame:
    """Return the ``n_genes`` rows with highest variance across samples.

    This is an unsupervised, global filter used only for EDA and to speed up
    training.  Strictly speaking it is a mild form of information leakage;
    it is acceptable for high-level exploration but should be documented.
    The supervised SelectKBest step is inside cross-validation to avoid
    leaking labels.
    """
    n = min(n_genes, matrix.shape[0])
    top = matrix.var(axis=1).nlargest(n).index
    return matrix.loc[top]


def load_clinical_sample(path: str) -> pd.DataFrame:
    """Load TCGA clinical sample sheet and strip header comment rows."""
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(
            f"Clinical sample file not found: {p.resolve()}")
    LOG.info("Loading clinical sample sheet from %s", p)
    # TCGA clinical files have 4 header rows: title, description, type, priority
    df = pd.read_csv(p, sep="\t", skiprows=4, low_memory=False)
    return df


def load_labels(
    clinical: pd.DataFrame,
    sample_id_column: str,
    label_column: str,
    min_class_size: int,
) -> pd.Series:
    """Extract sample-id -> label mapping, dropping missing/rare classes."""
    if sample_id_column not in clinical.columns:
        raise ValueError(f"Sample ID column '{sample_id_column}' not found")
    if label_column not in clinical.columns:
        raise ValueError(f"Label column '{label_column}' not found")

    labels = clinical[[sample_id_column, label_column]].dropna()
    labels = labels.set_index(sample_id_column)[label_column].astype(str)
    labels = labels.loc[labels.index != ""]

    counts = labels.value_counts()
    keep = counts[counts >= min_class_size].index
    labels = labels.loc[labels.isin(keep)]

    LOG.info("Loaded %d labels across %d classes",
             len(labels), labels.nunique())
    LOG.info("Class distribution:\n%s", labels.value_counts())
    return labels


# ---------------------------------------------------------------------------
# Preprocessing orchestration
# ---------------------------------------------------------------------------


def align_expression_and_labels(
    expression: pd.DataFrame, labels: pd.Series
) -> Tuple[pd.DataFrame, pd.Series]:
    """Return expression (samples x genes) and labels with matching sample IDs."""
    common = expression.columns.intersection(labels.index)
    if len(common) == 0:
        # Clinical IDs may be patient IDs without the -01 sample suffix.
        # Try matching by stripping the last segment.
        patient_ids = expression.columns.str.replace(
            r"-[^-]+$", "", regex=True)
        common = patient_ids.intersection(labels.index)
        if len(common) == 0:
            raise ValueError(
                "No matching sample IDs between expression and clinical data")
        patient_to_sample = dict(zip(patient_ids, expression.columns))
        aligned_expr = expression.loc[:, [
            patient_to_sample[c] for c in common]]
        aligned_labels = labels.loc[common]
        aligned_expr.columns = aligned_labels.index
    else:
        aligned_expr = expression.loc[:, common]
        aligned_labels = labels.loc[common]

    # transpose to samples x genes
    X = aligned_expr.T
    y = aligned_labels.loc[X.index]
    return X, y


def prepare_features(config: Config) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load and transform expression into ML-ready features.

    Returns
    -------
    raw_df, cleaned_df, transformed_df, filtered_df
    """
    raw_df = load_expression_matrix(config.expression_path)
    cleaned_df = clean_expression_matrix(raw_df)
    transformed_df = log2_transform(
        cleaned_df) if config.log_transform else cleaned_df
    filtered_df = select_top_variable_genes(
        transformed_df, config.top_variable_genes)
    return raw_df, cleaned_df, transformed_df, filtered_df


# ---------------------------------------------------------------------------
# Exploratory analysis: PCA
# ---------------------------------------------------------------------------


def plot_pca(
    X: pd.DataFrame,
    y: pd.Series,
    title: str = "PCA of tumours by subtype",
    save_path: str | None = None,
    palette: str = "Set2",
) -> Tuple[pd.DataFrame, Any]:
    """Compute PCA on X (samples x genes) and return PC dataframe + model."""
    from sklearn.decomposition import PCA

    # Center and scale before PCA
    Xs = StandardScaler().fit_transform(X)
    pca = PCA(n_components=min(10, X.shape[1]), random_state=42)
    pcs = pca.fit_transform(Xs)

    pc_df = pd.DataFrame(
        pcs[:, :3],
        columns=["PC1", "PC2", "PC3"],
        index=X.index,
    )
    pc_df["subtype"] = y.loc[pc_df.index].values

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    sns.scatterplot(
        data=pc_df,
        x="PC1",
        y="PC2",
        hue="subtype",
        palette=palette,
        s=70,
        alpha=0.8,
        ax=axes[0],
        edgecolor="black",
        linewidth=0.3,
    )
    axes[0].set_title(f"{title} (PC1 vs PC2)")
    axes[0].set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    axes[0].set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")

    sns.scatterplot(
        data=pc_df,
        x="PC1",
        y="PC3",
        hue="subtype",
        palette=palette,
        s=70,
        alpha=0.8,
        ax=axes[1],
        legend=False,
        edgecolor="black",
        linewidth=0.3,
    )
    axes[1].set_title(f"{title} (PC1 vs PC3)")
    axes[1].set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    axes[1].set_ylabel(f"PC3 ({pca.explained_variance_ratio_[2]*100:.1f}%)")

    plt.tight_layout()
    if save_path:
        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        LOG.info("PCA plot saved to %s", save_path)
    else:
        plt.show()
    plt.close(fig)
    return pc_df, pca


# ---------------------------------------------------------------------------
# Model pipeline
# ---------------------------------------------------------------------------


def build_selector(config: Config) -> Any:
    """Return a feature-selection transformer based on config."""
    if config.selector == "kbest":
        return SelectKBest(score_func=f_classif, k=min(config.selected_features, 10000))
    return SelectFdr(score_func=f_classif, alpha=config.alpha_threshold)


def build_classifier(config: Config) -> Any:
    """Return the classifier specified by config."""
    if config.classifier == "logistic":
        # multinomial is the default behaviour for LogisticRegression in modern
        # scikit-learn; omit multi_class to avoid deprecation warnings.
        return LogisticRegression(
            max_iter=config.max_iter,
            class_weight=config.class_weight,
            solver="lbfgs",
        )
    if config.classifier == "random_forest":
        return RandomForestClassifier(
            n_estimators=300,
            class_weight=config.class_weight,
            random_state=config.random_state,
            n_jobs=-1,
        )
    if config.classifier == "linear_svc":
        return LinearSVC(
            max_iter=config.max_iter,
            class_weight=config.class_weight,
            multi_class="crammer_singer",
            dual="auto",
        )
    raise ValueError(f"Unknown classifier: {config.classifier}")


def build_pipeline(config: Config) -> Pipeline:
    """Build sklearn Pipeline with selection, scaling, and classification.

    IMPORTANT: the ``select`` step is fit inside each CV fold, so the model
    never peeks at the test-set genes.
    """
    return Pipeline(
        [
            ("select", build_selector(config)),
            ("scale", StandardScaler()),
            ("clf", build_classifier(config)),
        ]
    )


def multi_metric_scorers() -> Dict[str, Any]:
    """Scorers used for cross-validation, robust to missing predicted classes."""
    from sklearn.metrics import make_scorer

    return {
        "balanced_accuracy": "balanced_accuracy",
        "f1_macro": "f1_macro",
        "precision_macro": make_scorer(
            precision_score, average="macro", zero_division=0
        ),
        "recall_macro": make_scorer(recall_score, average="macro", zero_division=0),
        "roc_auc_ovr_weighted": "roc_auc_ovr_weighted",
    }


def evaluate_cv(
    pipeline: Pipeline, X: pd.DataFrame, y: pd.Series, config: Config
) -> Dict[str, Any]:
    """Run stratified cross-validation and return scores."""
    cv = StratifiedKFold(
        n_splits=config.cv_splits, shuffle=True, random_state=config.random_state
    )
    scoring = multi_metric_scorers()
    scores = cross_validate(
        pipeline,
        X,
        y,
        cv=cv,
        scoring=scoring,
        n_jobs=-1,
        return_train_score=False,
    )
    summary = {k.replace("test_", ""): float(np.mean(v))
               for k, v in scores.items()}
    summary["fold_scores"] = {
        k.replace("test_", ""): [float(x) for x in v] for k, v in scores.items()
    }
    LOG.info("Cross-validation results:")
    for metric, value in summary.items():
        if metric == "fold_scores":
            continue
        LOG.info("  %-25s %.4f", metric, value)
    return summary


def evaluate_holdout(
    pipeline: Pipeline,
    X_train: pd.DataFrame,
    y_train: pd.Series,
    X_test: pd.DataFrame,
    y_test: pd.Series,
    classes: np.ndarray,
    output_dir: Path,
) -> Dict[str, Any]:
    """Fit on training data and report comprehensive test-set metrics."""
    pipeline.fit(X_train, y_train)
    y_pred = pipeline.predict(X_test)

    report = classification_report(
        y_test, y_pred, target_names=classes, output_dict=True, zero_division=0
    )

    # ROC-AUC one-vs-rest
    prob = pipeline.predict_proba(X_test) if hasattr(
        pipeline, "predict_proba") else None
    roc_auc_ovr = None
    roc_auc_pr = None
    if prob is not None:
        lb = LabelBinarizer().fit(y_train)
        y_test_bin = lb.transform(y_test)
        roc_auc_ovr = float(roc_auc_score(
            y_test_bin, prob, average="weighted", multi_class="ovr"))
        # PR-AUC one-vs-rest
        pr_aucs = []
        for i, cls in enumerate(lb.classes_):
            precision, recall, _ = precision_recall_curve(
                y_test_bin[:, i], prob[:, i])
            pr_aucs.append(auc(recall, precision))
        roc_auc_pr = float(np.average(
            pr_aucs, weights=np.sum(y_test_bin, axis=0)))

    metrics = {
        "balanced_accuracy": float(balanced_accuracy_score(y_test, y_pred)),
        "f1_macro": float(f1_score(y_test, y_pred, average="macro", zero_division=0)),
        "precision_macro": float(precision_score(y_test, y_pred, average="macro", zero_division=0)),
        "recall_macro": float(recall_score(y_test, y_pred, average="macro", zero_division=0)),
        "roc_auc_ovr_weighted": roc_auc_ovr,
        "pr_auc_ovr_weighted": roc_auc_pr,
        "classification_report": report,
    }

    LOG.info("Hold-out test metrics:")
    for k, v in metrics.items():
        if k == "classification_report":
            continue
        LOG.info("  %-25s %s", k, "N/A" if v is None else f"{v:.4f}")

    # Confusion matrix
    fig, ax = plt.subplots(figsize=(8, 7))
    ConfusionMatrixDisplay.from_predictions(
        y_test, y_pred, display_labels=classes, ax=ax, xticks_rotation=45
    )
    ax.set_title("Confusion matrix (hold-out test set)")
    cm_path = output_dir / "figures" / "confusion_matrix.png"
    cm_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(cm_path, dpi=300, bbox_inches="tight")
    LOG.info("Confusion matrix saved to %s", cm_path)
    plt.close(fig)

    return metrics


# ---------------------------------------------------------------------------
# Interpretation: top genes
# ---------------------------------------------------------------------------


def extract_top_genes(
    pipeline: Pipeline,
    X: pd.DataFrame,
    classes: np.ndarray,
    top_n: int = 50,
) -> pd.DataFrame:
    """Extract selected genes and their mean absolute coefficients/importance.

    For linear classifiers we average absolute multinomial coefficients.
    For tree-based models we use feature_importances_.
    """
    clf = pipeline.named_steps["clf"]
    selector = pipeline.named_steps["select"]

    selected_mask = selector.get_support()
    selected_genes = np.array(X.columns)[selected_mask]

    if hasattr(clf, "coef_"):
        # shape (n_classes, n_features) for multinomial logistic regression
        weights = np.asarray(clf.coef_)
        if weights.ndim == 2:
            importance = np.mean(np.abs(weights), axis=0)
        else:
            importance = np.abs(weights)
    elif hasattr(clf, "feature_importances_"):
        importance = np.asarray(clf.feature_importances_)
    else:
        importance = np.ones(len(selected_genes))

    df = pd.DataFrame(
        {
            "gene": selected_genes,
            "mean_abs_weight": importance,
        }
    ).sort_values("mean_abs_weight", ascending=False)

    # Add per-class weight direction for the top gene (optional detail)
    if hasattr(clf, "coef_") and np.asarray(clf.coef_).ndim == 2:
        coefs = np.asarray(clf.coef_)
        for i, cls in enumerate(classes):
            df[f"coef_{cls}"] = coefs[i, np.argsort(importance)[::-1]]

    return df.head(top_n).reset_index(drop=True)


def annotate_genes(gene_df: pd.DataFrame) -> pd.DataFrame:
    """Add lightweight biological notes for known thyroid / cancer genes.

    This is a curated, hackathon-friendly annotation. In production, use
    MyGene.info, MSigDB, or a pathway database (gseapy / ReactomePA).
    """
    notes = {
        "BRAF": "MAPK driver; BRAF V600E defines aggressive PTC / tall-cell variant.",
        "RAS": "MAPK/PI3K signalling; common in follicular-variant PTC and N-RAS.",
        "HRAS": "RAS family oncogene; follicular-variant and anaplastic disease.",
        "NRAS": "RAS family oncogene; associated with follicular variant PTC.",
        "KRAS": "RAS family oncogene; less common in classical PTC.",
        "TERT": "Telomerase reverse transcriptase; promoter mutations linked to poor prognosis.",
        "TP53": "Tumour suppressor; mutated in aggressive / anaplastic thyroid cancer.",
        "CDKN2A": "Cell-cycle inhibitor; loss linked to progression and poor outcome.",
        "PIK3CA": "PI3K/AKT pathway; implicated in aggressive thyroid carcinoma.",
        "AKT1": "PI3K/AKT signalling; promotes survival and therapy resistance.",
        "PTEN": "PI3K/AKT negative regulator; loss associated with aggressiveness.",
        "CTNNB1": "Beta-catenin; Wnt pathway activation in poorly differentiated disease.",
        "RET": "Rearrangements (RET/PTC) common in radiation-induced PTC.",
        "NTRK1": "Receptor tyrosine kinase fusions in subset of PTC.",
        "TG": "Thyroglobulin; thyroid differentiation marker, used for disease monitoring.",
        "TPO": "Thyroid peroxidase; thyroid differentiation / iodine metabolism.",
        "SLC5A5": "Sodium-iodide symporter (NIS); iodine-uptake therapy target.",
        "PAX8": "Thyroid lineage transcription factor; essential for differentiation.",
        "FOXE1": "Thyroid transcription factor; PTC susceptibility and differentiation.",
        "HMGA2": "Chromatin remodeler; over-expressed in PTC and follicular neoplasms.",
        "MET": "HGF receptor; linked to invasion and metastasis in thyroid cancer.",
        "MUC1": "Mucin; epithelial marker, can be elevated in tall-cell variant.",
        "CDH1": "E-cadherin; epithelial marker, loss associated with EMT.",
        "VIM": "Vimentin; mesenchymal marker, EMT / de-differentiation signal.",
        "FN1": "Fibronectin; stromal / EMT marker.",
        "S100A4": "Metastasis-promoting calcium-binding protein.",
        "KRT19": "Cytokeratin 19; diagnostic marker in thyroid nodules.",
        "HBEGF": "Heparin-binding EGF; EGFR ligand in aggressive PTC.",
        "EGFR": "Epidermal growth factor receptor; therapy target in anaplastic disease.",
        "VEGFA": "Angiogenesis driver; target of anti-angiogenic thyroid-cancer therapy.",
        "PDGFRA": "Receptor tyrosine kinase; implicated in thyroid tumourigenesis.",
        "KIT": "Stem-cell factor receptor; expressed in medullary thyroid cancer.",
        "ATM": "DNA-damage response; linked to radio-sensitivity.",
        "CHEK2": "DNA-damage checkpoint; tumour suppressor role.",
        "BRCA1": "DNA repair; rarely mutated but part of genomic instability signature.",
        "BRCA2": "DNA repair; homologous recombination deficiency context.",
    }
    gene_df["biological_note"] = gene_df["gene"].map(notes).fillna(
        "No curated note; investigate in MSigDB / PubMed / GeneCards."
    )
    return gene_df


def plot_top_genes(gene_df: pd.DataFrame, save_path: Path) -> None:
    """Bar plot of the top 20 most predictive genes."""
    top20 = gene_df.head(20).sort_values("mean_abs_weight", ascending=True)
    fig, ax = plt.subplots(figsize=(8, 10))
    sns.barplot(
        data=top20,
        y="gene",
        x="mean_abs_weight",
        hue="gene",
        palette="viridis",
        ax=ax,
        edgecolor="black",
        legend=False,
    )
    ax.set_title("Top 20 predictive genes (mean |coefficient|)")
    ax.set_xlabel("Mean absolute weight")
    ax.set_ylabel("Gene")
    plt.tight_layout()
    save_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    LOG.info("Top-genes plot saved to %s", save_path)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------


def report_data_shapes(
    raw_df: pd.DataFrame,
    cleaned_df: pd.DataFrame,
    transformed: pd.DataFrame,
    filtered_df: pd.DataFrame,
    X: pd.DataFrame,
    y: pd.Series,
) -> None:
    """Log dimension and shape diagnostics."""
    LOG.info("--- Data shape report ---")
    LOG.info("Raw expression shape:        %s (genes x samples)", raw_df.shape)
    LOG.info("Cleaned expression shape:    %s", cleaned_df.shape)
    LOG.info("Log-transformed shape:       %s", transformed.shape)
    LOG.info("Top-variable features shape: %s (genes x samples)",
             filtered_df.shape)
    LOG.info("ML matrix X shape:           %s (samples x features)", X.shape)
    LOG.info("Labels y length:             %d", len(y))
    LOG.info("Number of classes:           %d", y.nunique())


def save_results(
    config: Config,
    cv_scores: Dict[str, Any],
    test_metrics: Dict[str, Any],
    top_genes: pd.DataFrame,
    output_dir: Path,
) -> None:
    """Persist JSON report, tables, and figures."""
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "tables").mkdir(exist_ok=True)
    (output_dir / "figures").mkdir(exist_ok=True)
    (output_dir / "gene_lists").mkdir(exist_ok=True)

    # JSON report
    report: Dict[str, Any] = {
        "config": asdict(config),
        "cv_scores": cv_scores,
        "test_metrics": {
            k: v for k, v in test_metrics.items() if k != "classification_report"
        },
        "n_top_genes": len(top_genes),
    }
    report_path = output_dir / "tables" / "results_report.json"
    with open(report_path, "w") as fh:
        json.dump(report, fh, indent=2)
    LOG.info("JSON report saved to %s", report_path)

    # Top genes CSV
    genes_path = output_dir / "gene_lists" / "top_predictive_genes.csv"
    top_genes.to_csv(genes_path, index=False)
    LOG.info("Top genes saved to %s", genes_path)

    # Test classification report as CSV
    if "classification_report" in test_metrics:
        cr = test_metrics["classification_report"]
        cr_df = pd.DataFrame(cr).transpose()
        cr_path = output_dir / "tables" / "classification_report.csv"
        cr_df.to_csv(cr_path)
        LOG.info("Classification report saved to %s", cr_path)


def print_summary(
    config: Config,
    cv_scores: Dict[str, Any],
    test_metrics: Dict[str, Any],
    top_genes: pd.DataFrame,
) -> None:
    """Print a concise, hackathon-friendly summary to the console."""
    print("\n" + "=" * 70)
    print(" ENEZA: Thyroid carcinoma subtype classifier — summary")
    print("=" * 70)
    print(f"Model:        {config.model_id()}")
    print(
        f"Features:     top {config.top_variable_genes} variable genes → {config.selected_features} selected")
    print(f"CV splits:    {config.cv_splits}-fold stratified")
    print(f"Test split:   {config.test_size:.0%} held-out")
    print("\nCross-validation macro metrics:")
    for metric in ["balanced_accuracy", "f1_macro", "roc_auc_ovr_weighted"]:
        print(f"  {metric:25s} {cv_scores.get(metric, 0.0):.4f}")
    print("\nHold-out test metrics:")
    for metric in ["balanced_accuracy", "f1_macro", "roc_auc_ovr_weighted", "pr_auc_ovr_weighted"]:
        val = test_metrics.get(metric)
        print(f"  {metric:25s} {'N/A' if val is None else f'{val:.4f}'}")
    print("\nTop 10 predictive genes:")
    for _, row in top_genes.head(10).iterrows():
        print(f"  {row['gene']:15s} weight={row['mean_abs_weight']:.4f}")
    print("=" * 70 + "\n")


# ---------------------------------------------------------------------------
# Main workflow
# ---------------------------------------------------------------------------


def run_workflow(config: Config) -> Dict[str, Any]:
    """Run the full data → EDA → model → interpretation → report pipeline."""
    output_dir = Path(config.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # ---------------- Phase 1: load and clean ----------------
    raw_df, cleaned_df, transformed, filtered_df = prepare_features(config)

    # ---------------- Phase 2: labels ----------------
    clinical = load_clinical_sample(config.clinical_sample_path)
    labels = load_labels(
        clinical,
        config.sample_id_column,
        config.label_column,
        config.min_class_size,
    )

    # ---------------- Phase 3: align ----------------
    X, y = align_expression_and_labels(filtered_df, labels)
    report_data_shapes(raw_df, cleaned_df, transformed, filtered_df, X, y)

    # ---------------- Phase 4: EDA (PCA) ----------------
    LOG.info("Running PCA for exploratory visualisation ...")
    pc_df, pca_model = plot_pca(
        X,
        y,
        title="TCGA THCA histologic subtypes",
        save_path=str(output_dir / "figures" / "pca_subtype.png"),
    )
    var_explained = {
        f"PC{i+1}": float(ratio)
        for i, ratio in enumerate(pca_model.explained_variance_ratio_[:5])
    }
    LOG.info("Variance explained by top 5 PCs: %s", var_explained)

    # ---------------- Phase 5: train / test split ----------------
    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=config.test_size,
        random_state=config.random_state,
        stratify=y,
    )
    LOG.info(
        "Train size: %d | Test size: %d | Classes: %d",
        len(y_train),
        len(y_test),
        y_train.nunique(),
    )

    # ---------------- Phase 6: cross-validation ----------------
    pipeline = build_pipeline(config)
    LOG.info("Pipeline: %s", pipeline)
    cv_scores = evaluate_cv(pipeline, X_train, y_train, config)

    # ---------------- Phase 7: hold-out test ----------------
    classes = np.array(sorted(y.unique()))
    test_metrics = evaluate_holdout(
        pipeline,
        X_train,
        y_train,
        X_test,
        y_test,
        classes,
        output_dir,
    )

    # ---------------- Phase 8: interpretation ----------------
    top_genes = extract_top_genes(
        pipeline, X_train, classes, top_n=config.selected_features
    )
    top_genes = annotate_genes(top_genes)
    plot_top_genes(top_genes, output_dir / "figures" / "top_genes.png")

    # ---------------- Phase 9: save and report ----------------
    save_results(config, cv_scores, test_metrics, top_genes, output_dir)
    print_summary(config, cv_scores, test_metrics, top_genes)

    # ---------------- Stretch goal 1: feature-selection comparison ----------------
    LOG.info("\n--- Stretch goal 1: feature-selection comparison ---")
    selector_comparison = compare_feature_selectors(
        config, X_train, y_train, output_dir
    )

    # ---------------- Stretch goal 1b: classifier comparison ----------------
    LOG.info("\n--- Stretch goal 1b: classifier comparison ---")
    classifier_comparison = compare_models(
        config, X_train, y_train, output_dir)

    # ---------------- Stretch goal 2: survival analysis ----------------
    LOG.info("\n--- Stretch goal 2: survival analysis ---")
    survival_results = run_survival_analysis(
        config, X, y, pipeline, output_dir
    )

    return {
        "config": config,
        "cv_scores": cv_scores,
        "test_metrics": test_metrics,
        "top_genes": top_genes,
        "selector_comparison": selector_comparison,
        "classifier_comparison": classifier_comparison,
        "survival_results": survival_results,
        "pca_model": pca_model,
        "pc_df": pc_df,
        "pipeline": pipeline,
    }


def compare_models(
    config: Config,
    X_train: pd.DataFrame,
    y_train: pd.Series,
    output_dir: Path,
) -> pd.DataFrame:
    """Compare multiple classifiers with honest stratified CV.

    Only classifiers that provide probability estimates are included so that
    ROC-AUC and PR-AUC can be compared on the same footing. LinearSVC is
    excluded because it lacks ``predict_proba``.
    """
    classifiers: List[Literal["logistic", "random_forest"]] = [
        "logistic",
        "random_forest",
    ]
    cv = StratifiedKFold(
        n_splits=config.cv_splits, shuffle=True, random_state=config.random_state
    )
    rows: List[Dict[str, Any]] = []
    for clf_name in classifiers:
        cfg = Config(**{**asdict(config), "classifier": clf_name})
        pipe = build_pipeline(cfg)
        scores = cross_validate(
            pipe,
            X_train,
            y_train,
            cv=cv,
            scoring=multi_metric_scorers(),
            n_jobs=-1,
        )
        rows.append(
            {
                "classifier": clf_name,
                "balanced_accuracy": float(np.mean(scores["test_balanced_accuracy"])),
                "f1_macro": float(np.mean(scores["test_f1_macro"])),
                "precision_macro": float(np.mean(scores["test_precision_macro"])),
                "recall_macro": float(np.mean(scores["test_recall_macro"])),
                "roc_auc_ovr_weighted": float(
                    np.mean(scores["test_roc_auc_ovr_weighted"])
                ),
            }
        )

    df = pd.DataFrame(rows).sort_values("f1_macro", ascending=False)
    df.to_csv(output_dir / "tables" / "classifier_comparison.csv", index=False)
    LOG.info("Classifier comparison:\n%s", df.to_string(index=False))

    fig, ax = plt.subplots(figsize=(8, 5))
    plot_df = df.melt(
        id_vars="classifier",
        value_vars=["balanced_accuracy", "f1_macro", "roc_auc_ovr_weighted"],
        var_name="metric",
        value_name="score",
    )
    sns.barplot(data=plot_df, x="classifier", y="score", hue="metric", ax=ax)
    ax.set_ylim(0, 1)
    ax.set_title("Classifier comparison (5-fold CV)")
    ax.set_ylabel("Score")
    ax.set_xlabel("Classifier")
    plt.tight_layout()
    save_path = output_dir / "figures" / "classifier_comparison.png"
    save_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    LOG.info("Classifier comparison plot saved to %s", save_path)
    plt.close(fig)
    return df


# ---------------------------------------------------------------------------
# Stretch goal 1: compare feature-selection strategies
# ---------------------------------------------------------------------------


def build_feature_selection_pipeline(
    selector_name: Literal["kbest", "fdr", "mutual_info"],
    k_features: int = 100,
    alpha_fdr: float = 0.05,
    class_weight: str = "balanced",
    random_state: int = 42,
) -> Pipeline:
    """Build a pipeline that differs only in the feature-selection step.

    This lets us fairly compare selection strategies while keeping the final
    classifier identical (multinomial logistic regression).
    """
    if selector_name == "kbest":
        select_step = SelectKBest(score_func=f_classif, k=k_features)
    elif selector_name == "fdr":
        select_step = SelectFdr(score_func=f_classif, alpha=alpha_fdr)
    elif selector_name == "mutual_info":
        # Captures non-linear associations without the convergence cost of LASSO.
        # Setting a random_state makes the k-NN entropy estimate reproducible.
        select_step = SelectKBest(
            score_func=partial(mutual_info_classif, random_state=random_state),
            k=k_features,
        )
    else:
        raise ValueError(f"Unknown selector: {selector_name}")

    return Pipeline(
        [
            ("select", select_step),
            ("scale", StandardScaler()),
            (
                "clf",
                LogisticRegression(
                    max_iter=10_000,
                    class_weight=class_weight,
                    solver="lbfgs",
                    random_state=random_state,
                ),
            ),
        ]
    )


def compare_feature_selectors(
    config: Config,
    X_train: pd.DataFrame,
    y_train: pd.Series,
    output_dir: Path,
) -> pd.DataFrame:
    """Compare univariate ANOVA (k-best), FDR-corrected ANOVA, mutual
    information, and DESeq/edgeR differential-expression ranking.

    All pipelines use the same final classifier (multinomial logistic
    regression) and the same 5-fold stratified CV, so differences are due to
    the feature-selection method alone.

    Returns a leaderboard DataFrame ordered by macro F1. The plot is saved to
    ``output_dir/figures/feature_selector_comparison.png``.

    Interpretation of the statistics
    ----------------------------------
    * **Balanced accuracy** — average per-class recall; robust to the tall-cell
      class being small.
    * **Macro F1** — harmonic mean of precision and recall averaged across the
      three PTC variants; our primary comparison metric.
    * **ROC-AUC OvR weighted** — ability to rank samples correctly within each
      class, weighted by class size.
    * **n_selected** — median number of genes retained by the selector across
      CV folds. ``fdr`` may vary; ``kbest``, ``mutual_info`` and ``deseq`` are
      fixed at the configured ``k``.
    """
    selectors: List[Literal["kbest", "fdr", "mutual_info", "deseq"]] = [
        "kbest",
        "fdr",
        "mutual_info",
    ]
    if Path(config.deseq_results_dir).is_dir():
        selectors.append("deseq")
    else:
        LOG.warning(
            "DESeq results directory not found at %s; omitting from comparison",
            config.deseq_results_dir,
        )
    cv = StratifiedKFold(
        n_splits=config.cv_splits, shuffle=True, random_state=config.random_state
    )
    rows: List[Dict[str, Any]] = []

    for sel_name in selectors:
        if sel_name == "deseq":
            pipe = build_deseq_pipeline(config)
            n_selected = config.deseq_k
        else:
            pipe = build_feature_selection_pipeline(
                selector_name=sel_name,
                k_features=config.selected_features,
                alpha_fdr=config.alpha_threshold,
            )
            n_selected = _median_selected_features(pipe, X_train, y_train, cv)
        scores = cross_validate(
            pipe,
            X_train,
            y_train,
            cv=cv,
            scoring=multi_metric_scorers(),
            n_jobs=-1,
        )
        rows.append(
            {
                "selector": sel_name,
                "balanced_accuracy": float(np.mean(scores["test_balanced_accuracy"])),
                "f1_macro": float(np.mean(scores["test_f1_macro"])),
                "precision_macro": float(np.mean(scores["test_precision_macro"])),
                "recall_macro": float(np.mean(scores["test_recall_macro"])),
                "roc_auc_ovr_weighted": float(np.mean(scores["test_roc_auc_ovr_weighted"])),
                "n_selected": n_selected,
            }
        )

    df = pd.DataFrame(rows).sort_values("f1_macro", ascending=False)
    df.to_csv(output_dir / "tables" /
              "feature_selector_comparison.csv", index=False)
    LOG.info("Feature-selection comparison:\n%s", df.to_string(index=False))

    # Plot
    fig, ax = plt.subplots(figsize=(10, 5))
    plot_df = df.melt(
        id_vars="selector",
        value_vars=["balanced_accuracy", "f1_macro", "roc_auc_ovr_weighted"],
        var_name="metric",
        value_name="score",
    )
    sns.barplot(data=plot_df, x="selector", y="score", hue="metric", ax=ax)
    ax.set_ylim(0, 1)
    ax.set_title("Feature-selection strategy comparison (5-fold CV)")
    ax.set_ylabel("Score")
    ax.set_xlabel("Feature selector")
    plt.tight_layout()
    save_path = output_dir / "figures" / "feature_selector_comparison.png"
    save_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    LOG.info("Feature-selector comparison plot saved to %s", save_path)
    plt.close(fig)
    return df


def _median_selected_features(
    pipeline: Pipeline, X: pd.DataFrame, y: pd.Series, cv: StratifiedKFold
) -> int:
    """Return the median number of features selected by ``pipeline`` across CV folds."""
    counts = []
    for train_idx, _ in cv.split(X, y):
        pipeline.fit(X.iloc[train_idx], y.iloc[train_idx])
        counts.append(int(pipeline.named_steps["select"].get_support().sum()))
    return int(np.median(counts))


# ---------------------------------------------------------------------------
# DESeq / edgeR differential-expression feature selection (used inside the
# unified feature-selector comparison above)
# ---------------------------------------------------------------------------


def load_deseq_results(deseq_dir: str) -> Dict[str, pd.DataFrame]:
    """Load all CSV files of pairwise DESeq/edgeR results from a directory.

    Expected columns include ``Hugo_Symbol``, ``logFC``, ``logCPM``, ``LR``,
    ``PValue``, ``FDR``. Each file should represent one pairwise comparison
    (e.g., TallCell vs Classical).
    """
    p = Path(deseq_dir)
    if not p.is_dir():
        raise FileNotFoundError(
            f"DESeq results directory not found: {p.resolve()}")
    files = sorted(p.glob("*.csv"))
    if not files:
        raise FileNotFoundError(f"No CSV files found in {p.resolve()}")

    results: Dict[str, pd.DataFrame] = {}
    for f in files:
        df = pd.read_csv(f)
        if "Hugo_Symbol" not in df.columns:
            raise ValueError(f"DESeq file {f} missing 'Hugo_Symbol' column")
        results[f.stem] = df.set_index("Hugo_Symbol")
        LOG.info("Loaded DESeq results: %s (%d genes)",
                 f.name, len(results[f.stem]))
    return results


def rank_deseq_genes(
    deseq_results: Dict[str, pd.DataFrame],
    expression_genes: pd.Index,
    k: int = 100,
    score_col: str = "FDR",
    method: str = "min",
) -> List[str]:
    """Rank genes by their most significant appearance across pairwise comparisons.

    Parameters
    ----------
    deseq_results
        Mapping from comparison name to DESeq result DataFrame indexed by
        ``Hugo_Symbol``.
    expression_genes
        Gene symbols present in the expression matrix (after cleaning).
    k
        Number of top genes to return.
    score_col
        Column to use for ranking (typically ``FDR`` or ``PValue``). Lower is
        more significant.
    method
        How to combine scores across comparisons:

        * ``min`` - use the smallest (most significant) score per gene. This
          selects genes that are significant in at least one comparison.
        * ``max`` - use the largest (least significant) score per gene. This
          is more stringent: a gene must be significant in all comparisons.
        * ``mean`` - average the scores across comparisons.

    Returns
    -------
    List of the top ``k`` gene symbols, ordered by the combined score.
    """
    gene_scores: Dict[str, List[float]] = {}
    for name, df in deseq_results.items():
        if score_col not in df.columns:
            raise ValueError(
                f"DESeq result '{name}' missing column '{score_col}'")
        for gene, row in df.iterrows():
            score = float(row[score_col])
            gene_scores.setdefault(str(gene), []).append(score)

    if not gene_scores:
        raise ValueError("No gene scores found in DESeq results")

    if method == "min":
        combined = {g: min(scores) for g, scores in gene_scores.items()}
    elif method == "max":
        combined = {g: max(scores) for g, scores in gene_scores.items()}
    elif method == "mean":
        combined = {g: float(np.mean(scores))
                    for g, scores in gene_scores.items()}
    else:
        raise ValueError(f"Unknown DESeq rank method: {method}")

    score_series = pd.Series(combined)
    score_series = score_series.loc[score_series.index.isin(expression_genes)]
    score_series = score_series.dropna()
    score_series = score_series.sort_values(ascending=True)
    selected = score_series.index[:k].tolist()
    LOG.info(
        "DESeq selected %d genes (from %d available in expression) using %s/%s",
        len(selected),
        len(score_series),
        score_col,
        method,
    )
    return selected


class DESeqFeatureSelector(BaseEstimator, TransformerMixin):
    """sklearn-compatible selector that picks the top k genes from DESeq results.

    Notes
    -----
    The selected genes are determined entirely from the pre-computed DESeq CSV
    files, not from the training labels. This is fast and consistent across CV
    folds, but it assumes the DESeq results were computed on a representative,
    batch-corrected dataset. If the DESeq results were derived from the same
    data that is being cross-validated, this is a mild form of information
    leakage; we document this in the report.
    """

    def __init__(
        self,
        deseq_results_dir: str,
        k: int = 100,
        score_col: str = "FDR",
        rank_method: str = "min",
    ):
        self.deseq_results_dir = deseq_results_dir
        self.k = k
        self.score_col = score_col
        self.rank_method = rank_method

    def fit(self, X: pd.DataFrame, y=None) -> "DESeqFeatureSelector":
        deseq_results = load_deseq_results(self.deseq_results_dir)
        self.selected_genes_ = rank_deseq_genes(
            deseq_results,
            expression_genes=X.columns,
            k=self.k,
            score_col=self.score_col,
            method=self.rank_method,
        )
        self.feature_names_in_ = np.array(X.columns)
        self.support_ = np.isin(X.columns, self.selected_genes_)
        return self

    def transform(self, X: pd.DataFrame) -> np.ndarray:
        return X.loc[:, self.selected_genes_].values

    def get_support(self, indices: bool = False):
        if indices:
            return np.where(self.support_)[0]
        return self.support_

    def get_feature_names_out(self, input_features=None):
        return np.array(self.selected_genes_)


def build_deseq_pipeline(config: Config) -> Pipeline:
    """Pipeline using DESeq/edgeR-derived genes instead of ANOVA F-test."""
    return Pipeline(
        [
            (
                "select",
                DESeqFeatureSelector(
                    deseq_results_dir=config.deseq_results_dir,
                    k=config.deseq_k,
                    score_col=config.deseq_score_col,
                    rank_method=config.deseq_rank_method,
                ),
            ),
            ("scale", StandardScaler()),
            (
                "clf",
                LogisticRegression(
                    max_iter=10_000,
                    class_weight=config.class_weight,
                    solver="lbfgs",
                    random_state=config.random_state,
                ),
            ),
        ]
    )


def compare_deseq_vs_anova(
    config: Config,
    X_train: pd.DataFrame,
    y_train: pd.Series,
    output_dir: Path,
) -> pd.DataFrame:
    """Compare ANOVA k-best vs. DESeq/edgeR pre-ranked feature selection.

    Both pipelines use the same final classifier (multinomial logistic
    regression) and the same stratified CV, so the comparison isolates the
    effect of the feature-selection method.
    """
    cv = StratifiedKFold(
        n_splits=config.cv_splits, shuffle=True, random_state=config.random_state
    )
    rows: List[Dict[str, Any]] = []

    # ANOVA baseline: same k as DESeq for a fair comparison
    anova_cfg = Config(
        **{**asdict(config), "selector": "kbest", "selected_features": config.deseq_k}
    )
    anova_pipe = build_pipeline(anova_cfg)
    anova_scores = cross_validate(
        anova_pipe,
        X_train,
        y_train,
        cv=cv,
        scoring=multi_metric_scorers(),
        n_jobs=-1,
    )
    rows.append(
        {
            "method": "anova_kbest",
            "k": config.deseq_k,
            "balanced_accuracy": float(np.mean(anova_scores["test_balanced_accuracy"])),
            "f1_macro": float(np.mean(anova_scores["test_f1_macro"])),
            "precision_macro": float(np.mean(anova_scores["test_precision_macro"])),
            "recall_macro": float(np.mean(anova_scores["test_recall_macro"])),
            "roc_auc_ovr_weighted": float(
                np.mean(anova_scores["test_roc_auc_ovr_weighted"])
            ),
        }
    )

    # DESeq pipeline
    deseq_pipe = build_deseq_pipeline(config)
    deseq_scores = cross_validate(
        deseq_pipe,
        X_train,
        y_train,
        cv=cv,
        scoring=multi_metric_scorers(),
        n_jobs=-1,
    )
    rows.append(
        {
            "method": "deseq",
            "k": config.deseq_k,
            "balanced_accuracy": float(np.mean(deseq_scores["test_balanced_accuracy"])),
            "f1_macro": float(np.mean(deseq_scores["test_f1_macro"])),
            "precision_macro": float(np.mean(deseq_scores["test_precision_macro"])),
            "recall_macro": float(np.mean(deseq_scores["test_recall_macro"])),
            "roc_auc_ovr_weighted": float(
                np.mean(deseq_scores["test_roc_auc_ovr_weighted"])
            ),
        }
    )

    df = pd.DataFrame(rows).sort_values("f1_macro", ascending=False)
    df.to_csv(output_dir / "tables" / "deseq_vs_anova.csv", index=False)
    LOG.info("DESeq vs ANOVA comparison:\n%s", df.to_string(index=False))

    fig, ax = plt.subplots(figsize=(8, 5))
    plot_df = df.melt(
        id_vars="method",
        value_vars=["balanced_accuracy", "f1_macro", "roc_auc_ovr_weighted"],
        var_name="metric",
        value_name="score",
    )
    sns.barplot(data=plot_df, x="method", y="score", hue="metric", ax=ax)
    ax.set_ylim(0, 1)
    ax.set_title("DESeq vs ANOVA feature selection (5-fold CV)")
    ax.set_ylabel("Score")
    ax.set_xlabel("Feature selection method")
    plt.tight_layout()
    save_path = output_dir / "figures" / "deseq_vs_anova.png"
    save_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    LOG.info("DESeq vs ANOVA plot saved to %s", save_path)
    plt.close(fig)
    return df


# ---------------------------------------------------------------------------
# Stretch goal 2: survival analysis
# ---------------------------------------------------------------------------


def load_patient_survival(path: str) -> pd.DataFrame:
    """Load overall-survival (OS) data from the patient clinical file.

    TCGA encodes status as ``0:LIVING`` or ``1:DECEASED`` and time in months.
    We parse these into ``event`` (1 = death, 0 = censored) and ``time_months``.
    """
    if not Path(path).is_file():
        raise FileNotFoundError(f"Patient file not found: {path}")
    #  this loads the patients data skipping rows which are ideally metadata
    df = pd.read_csv(path, sep="\t", skiprows=4, low_memory=False)

    # data  isnt that big try to avoid error types so pd recommends  low_memory
    required = {"PATIENT_ID", "OS_STATUS", "OS_MONTHS"}
    # TODO: potential bug in this is if we  casing issues
    # safety net convert all columns to uppercase
    df.columns = df.columns.str.upper()
    # diff to get the missing required items and raise errors
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Patient file missing columns: {missing}")

    surv = (
        df[["PATIENT_ID", "OS_STATUS", "OS_MONTHS"]]
        .dropna()
        .copy()
    )
    # Status format examples: "0:LIVING", "1:DECEASED"
    surv["event"] = surv["OS_STATUS"].str.split(":").str[0].astype(int)
    surv["time_months"] = pd.to_numeric(surv["OS_MONTHS"], errors="coerce")
    surv = surv.dropna(subset=["event", "time_months"])
    surv = surv.set_index("PATIENT_ID")[["event", "time_months"]]
    LOG.info("Loaded survival data for %d patients", len(surv))
    return surv


def align_survival_to_samples(
    survival: pd.DataFrame, sample_ids: pd.Index
) -> pd.DataFrame:
    """Map patient-level survival to sample IDs by stripping the sample suffix."""
    patient_map = pd.DataFrame({"sample_id": sample_ids})
    patient_map["patient_id"] = patient_map["sample_id"].str.replace(
        r"-[^-]+$", "", regex=True
    )
    merged = patient_map.merge(
        survival.reset_index().rename(columns={"PATIENT_ID": "patient_id"}),
        on="patient_id",
        how="inner",
        validate="many_to_one",
    )
    merged = merged.set_index("sample_id")[["event", "time_months"]]
    return merged


def plot_kaplan_meier(
    survival_df: pd.DataFrame,
    group_col: str,
    title: str,
    save_path: Path,
    palette: str = "Set2",
) -> Dict[str, float]:
    """Plot Kaplan-Meier curves by ``group_col`` and return log-rank p-values.

    Statistics explained
    --------------------
    * **Kaplan-Meier estimate** — non-parametric estimate of the survival
      function S(t) = P(alive at time t). Censored observations contribute
      up to their last follow-up time and are then removed from the risk set.
    * **Log-rank test** — compares the observed number of events in each group
      to the number expected if all groups had the same survival curve. A small
      p-value (< 0.05) indicates that survival differs significantly between
      groups.

    documentation for using kaplanmeir
    https://lifelines.readthedocs.io/en/latest/fitters/univariate/KaplanMeierFitter.html
    """
    from lifelines import KaplanMeierFitter
    from lifelines.statistics import logrank_test

    fig, ax = plt.subplots(figsize=(9, 6))
    kmf = KaplanMeierFitter()
    groups = survival_df[group_col].unique()
    colors = sns.color_palette(palette, n_colors=len(groups))

    for i, grp in enumerate(sorted(groups)):
        sub = survival_df[survival_df[group_col] == grp]
        kmf.fit(
            sub["time_months"],
            event_observed=sub["event"],
            label=f"{grp} (n={len(sub)})",
        )
        kmf.plot_survival_function(ax=ax, ci_show=True, color=colors[i])

    ax.set_title(title)
    ax.set_xlabel("Time (months)")
    ax.set_ylabel("Overall survival probability")
    ax.set_ylim(0, 1.05)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    save_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    LOG.info("Kaplan-Meier plot saved to %s", save_path)
    plt.close(fig)

    # Pairwise log-rank tests

    # mhh more like  a chi test observed vs expected interesting!!!
    # The log-rank test is vital because it prevents biased conclusions when analyzing time-vent data
    logrank_results: Dict[str, float] = {}
    groups = sorted(survival_df[group_col].unique())
    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):
            g1 = survival_df[survival_df[group_col] == groups[i]]
            g2 = survival_df[survival_df[group_col] == groups[j]]
            result = logrank_test(
                g1["time_months"],
                g2["time_months"],
                event_observed_A=g1["event"],
                event_observed_B=g2["event"],
            )
            key = f"{groups[i]} vs {groups[j]}"
            logrank_results[key] = float(result.p_value)
    return logrank_results


def run_survival_analysis(
    config: Config,
    X: pd.DataFrame,
    y: pd.Series,
    pipeline: Pipeline,
    output_dir: Path,
) -> pd.DataFrame:
    """Link predicted/annotated subtypes to overall survival.

    Returns a DataFrame of log-rank p-values and writes KM curves for both the
    pathologist-assigned subtype and the model-predicted subtype.
    """
    LOG.info("Running survival analysis ...")
    survival = load_patient_survival(config.patient_path)
    surv_aligned = align_survival_to_samples(survival, X.index)

    # Build a master survival table
    surv_df = surv_aligned.copy()
    surv_df["annotated_subtype"] = y.loc[surv_df.index]
    surv_df["predicted_subtype"] = pipeline.predict(X.loc[surv_df.index])

    # KM by annotated subtype
    p_annotated = plot_kaplan_meier(
        surv_df,
        "annotated_subtype",
        "Kaplan-Meier overall survival by annotated subtype",
        output_dir / "figures" / "km_annotated_subtype.png",
    )

    # KM by predicted subtype
    p_predicted = plot_kaplan_meier(
        surv_df,
        "predicted_subtype",
        "Kaplan-Meier overall survival by predicted subtype",
        output_dir / "figures" / "km_predicted_subtype.png",
    )

    # Summarise
    summary = []
    for key, p in p_annotated.items():
        summary.append(
            {"comparison": key, "type": "annotated_subtype", "logrank_p": p}
        )
    for key, p in p_predicted.items():
        summary.append(
            {"comparison": key, "type": "predicted_subtype", "logrank_p": p}
        )
    summary_df = pd.DataFrame(summary)
    summary_df["significant_0.05"] = summary_df["logrank_p"] < 0.05
    summary_df.to_csv(output_dir / "tables" /
                      "survival_logrank.csv", index=False)
    LOG.info("Survival log-rank results:\n%s",
             summary_df.to_string(index=False))
    return summary_df


def main() -> None:
    """Entry point: base pipeline + stretch goals."""
    config = Config(
        selector="kbest",
        selected_features=100,
        classifier="logistic",
        cv_splits=5,
        test_size=0.20,
        random_state=42,
        top_variable_genes=5000,
        log_transform=True,
    )
    setup_logging(config.verbosity)
    run_workflow(config)


if __name__ == "__main__":
    main()
