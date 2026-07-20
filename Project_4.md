# HIV drug-resistance prediction from viral sequences

Kenya runs one of the largest HIV treatment programmes in the world, and keeping patients on effective therapy is central to controlling the epidemic. HIV mutates rapidly, and drug-resistance mutations can render first-line regimens ineffective, leading to treatment failure and onward transmission of resistant virus. The Stanford HIV Drug Resistance Database provides open, expertly curated data linking viral mutations to resistance for each antiretroviral drug. In this project you will learn from that data to identify which mutations drive resistance and to build an interpretable model that predicts resistance from a mutation profile. A key part of the work is checking that what the model learns matches established biology, and discussing how such a tool could support regimen selection in Kenya.

## Your task

1. Identify which RT/protease/integrase mutations are most associated with resistance to specific antiretrovirals.
2. Predict resistance (resistant/susceptible or a score) from the mutation profile.
3. Compare classical models and check whether the drivers match the literature.
4. Discuss implications for regimen selection in Kenya.

## Datasets

- **Stanford HIV Drug Resistance DB** — <https://hivdb.stanford.edu>. **Subset:** one drug class.
- **Los Alamos HIV Sequence DB** — <https://www.hiv.lanl.gov>. **Subset:** context.
- **NCBI GenBank (E-utilities)** — <https://www.ncbi.nlm.nih.gov/nuccore>. **Subset:** extra seqs.

## Deliverables

Reproducible notebook, mutation–resistance association analysis, an interpretable classifier, report with clinical framing + per-member roles.

## Stretch goals

Multiple drug classes; compare to Stanford rules-based scores; a small 'enter a genotype → predicted resistance' demo.

---
*Work in a shared Git repo, split tasks via issues, and make your analysis fully
reproducible (pinned environment + scripted data download). Every member must
understand the whole project and state their role in the report.*
