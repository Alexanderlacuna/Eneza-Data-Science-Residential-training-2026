# Child health outcomes

Maternal and child health is a national priority in Kenya, and the difference between a safe and a dangerous pregnancy often comes down to factors that can be measured and acted on early (a mother's blood pressure and blood sugar, and her access to skilled care). Data science can help both by flagging who is most at risk and by revealing where health services are falling short for particular communities. In this project you will work with open health data at two levels: an individual-level clinical dataset to build a maternal-risk classifier from routine measurements, and open Kenyan survey indicators to study how outcomes such as child stunting, skilled birth attendance, and immunisation coverage vary across regions and wealth groups. A central theme is doing this *responsibly*, by checking that the model is fair across groups and well-calibrated, and not merely accurate.

## Your task

1. Describe how Kenyan maternal/child indicators (stunting, skilled birth attendance, immunisation) vary by region and wealth quintile.
2. Discuss where to target interventions.

## Datasets

- **DHS API (STATcompiler)** — <https://api.dhsprogram.com>. **Subset:** Kenya indicators by region/wealth.
- **UNICEF / WHO GHO** — <https://data.unicef.org>. **Subset:** supporting indicators.

## Deliverables

Reproducible notebook, disaggregated EDA, report with ethics/Data Protection notes + roles.

## Stretch goals

Multidimensional index; compare DHS indicator waves; small-area estimation.

---
*Work in a shared Git repo, split tasks via issues, and make your analysis fully
reproducible (pinned environment + scripted data download). Every member must
understand the whole project and state their role in the report.*
