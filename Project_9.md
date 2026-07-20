#  Biomedical literature mining (NLP) 

## Abstract
It is hard to know everything on a topic in biomedical science at a specific time. As scientists, we usually have to just read and hope to grasp all the amount of knowledge humanly possible. This leaves critical gaps in healthcare interventions in the sense that information is not uniformly spread. Even health practitioners do not all have the same level of information. Although search tools (Google) can alleviate this problem, they happen more than once to fail to give reliable and specific information. In particular, the AMR cris would probably be better managed with an accurate knowledge of current resistance profiles and effective treatments. This impact could last longer in time if our knowledge was cumulated with the power of predictive analytics. Will building a knowledge base cataloguing all we know about AMR crisis in East Africa help?

## Objectives

1. Mine biomedical literature available on AMR in East Africa up to 2026
2. Build a knowledge base from the literature gathered on AMR
3. Develop a RAG (Retrieval Augmented Generation) to predict (risk of) resistance of bacterial strains to common classes of antiobiotics based on the knowledge base
4. Test performance of the RAG system over traditional Google search

## Deliverables:
- Knowledge base (database) of AMR in parseable format (json, etc.)
- Working RAG scripts that takes antibiotic and predicts its current resistance profile to specific drug based on knowledge base
- Performance evaluation results

---
*Work in a shared Git repo, split tasks via issues, and make your analysis fully
reproducible (pinned environment + scripted data download). Every member must
understand the whole project and state their role in the report.*
