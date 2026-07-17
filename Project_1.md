# Malaria parasite genomic surveillance

Malaria remains one of Kenya's heaviest disease burdens, causing millions of cases and thousands of deaths every year, with children under five and pregnant women mostly at risk. Decades of progress (insecticide-treated nets, rapid diagnostic tests, and artemisinin-based combination therapies) are now threatened by *Plasmodium falciparum* parasites that evolve resistance to the drugs/chemicals used against them. Because resistance can spread silently across regions before it is noticed in the clinic, health systems need fast, scalable ways to detect and track it. Parasite genomics allows tracking of these mutations. By reading the parasite's DNA we can identify the specific mutations that confer drug resistance and watch how their frequency rises and spreads across time and geography (genomic surveillance). MalariaGEN has openly released whole-genome data for tens of thousands of parasites from across Africa and beyond. In this project you will use that resource to map resistance markers, uncover how parasite populations are structured, and build a model that reads resistance directly from the genome.

## Your task

1. Quantify the frequency of known resistance markers (e.g. *pfkelch13*, *pfcrt*, *pfmdr1*, *dhfr/dhps*) across African samples.
2. Compare marker frequencies by region and over time.
3. Recover population structure with PCA / clustering.
4. Predict a resistance-associated label from genotype features and interpret it.

## Datasets

- **MalariaGEN Pf data** — <https://www.malariagen.net> via the `malariagen_data` Python client. **Subset:** a sample subset + specific resistance loci.
- **PlasmoDB (annotation)** — <https://plasmodb.org>. **Subset:** genes of interest.

## Deliverables

Reproducible notebook (data access → resistance profiling → analysis), region/time frequency plots, a classifier or clustering result with honest evaluation, report with limitations (sampling bias) and per-member roles.

## Stretch goals

Time-trend of *pfkelch13* (artemisinin resistance); East vs West Africa comparison; link to treatment-policy implications.

## useful links
https://malariagen.github.io/malariagen-data-python/latest/ 

https://www.who.int/teams/global-malaria-programme 

https://academic.oup.com/nar/article/50/D1/D898/6413610

https://veupathdb.org


---
*Work in a shared Git repo, split tasks via issues, and make your analysis fully
reproducible (pinned environment + scripted data download). Every member must
understand the whole project and state their role in the report.*
