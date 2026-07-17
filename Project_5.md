# Clinical text mining: insight from medical notes

A large share of health information never becomes neat, analysable data: it lives in free text (clinical notes, discharge summaries, referral letters, and the biomedical literature). Turning that unstructured text into usable, structured information is the job of clinical natural language processing (NLP), a skill covered directly in the training and increasingly important as health systems digitise their records. In this project you will mine free-text clinical records to recover useful structure eg example predicting the medical specialty of a note, extracting the conditions or medications it mentions, or classifying the type of note, and then examine where the model fails. Clinical text is full of traps: rare conditions with few examples, heavy use of abbreviations, and negations such as 'no evidence of malaria' that completely flip meaning. Confronting these challenges, and the ethics of working with health text, is central to the project.

## Your task

1. Recover structure from clinical text; predict specialty/diagnosis, extract medications/conditions, or classify note type.
2. Compare classical NLP (using TF-IDF and linear models) and see where it breaks.
3. Surface biases/safety risks (rare conditions; negation like 'no evidence of').
4. Discuss what responsible clinical NLP deployment would require.

## Datasets

- **Medical Transcriptions (mtsamples)** — <https://raw.githubusercontent.com/socd06/medical-nlp/master/data/mtsamples.csv>. **Subset:** a subset of specialties.
- **PubMed abstracts (E-utilities)** — <https://www.ncbi.nlm.nih.gov/pubmed>. **Subset:** biomedical text for NER practice.

## Deliverables

Reproducible notebook (data → clean → model → error analysis), a working NLP component, a frank error/bias analysis, report with a data-governance/ethics section and per-member roles.

## Stretch goals

Negation/uncertainty handling; a clinical-entity extractor; (HPC, stretch only) fine-tune a clinical transformer and compare.

---
*Work in a shared Git repo, split tasks via issues, and make your analysis fully
reproducible (pinned environment + scripted data download). Every member must
understand the whole project and state their role in the report.*
