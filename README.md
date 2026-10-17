# Pan-cancer immune-related adverse event reporting

R analysis code for **“Five latent patterns of immune checkpoint inhibitor-related adverse event reporting across cancer types: an integrated pharmacovigilance and transcriptomic analysis.”**

## Data sources

- **FAERS:** US FDA Adverse Event Reporting System, reports from 2015–2025 ([FDA](https://www.fda.gov/drugs/questions-and-answers-fdas-adverse-event-reporting-system-faers/fda-adverse-event-reporting-system-faers-quarterly-data-extract-files)).
- **JADER:** Japanese Adverse Drug Event Report database, reports from 2015–2025 ([PMDA](https://www.pmda.go.jp/safety/info-services/drugs/adr-info/suspected-adr/0003.html)).
- **TCGA:** tumor RNA-seq, clinical and mutation data from the [NCI Genomic Data Commons](https://portal.gdc.cancer.gov/).
- **GTEx:** Analysis v11 normal-tissue expression data and sample annotations from the [GTEx Portal](https://gtexportal.org/home/downloads/adult-gtex).

Additional published TCGA PanImmune T-cell receptor data and gene sets are described in the manuscript and analysis scripts. Original database downloads are not included in this repository.

## Study-specific reference files

| File | Purpose |
| --- | --- |
| `Data/irAE列表.xlsx` | Study irAE preferred-term (PT) list and system organ class (SOC) assignments. |
| `Results/癌种映射.xlsx` | Reference mapping between primary organs, histological categories and TCGA entities. This is an input file, retained at the path expected by the scripts. |
| `Data/JADER/药物名查询表格.xlsx` | Japanese–English drug-name mapping used for JADER preprocessing, with Chinese labels. The release copy preserves the original name columns and row order. Unused external-lookup formulas and the redundant working sheet have been removed; the unused fourth column is retained for script compatibility. |

Additional study-specific cancer and drug classification rules are included in `Code/mapping_rules.R`. These reference files describe terminology and classification, not patient-level records.

The full MedDRA Japanese–English dictionary (`Data/JADER/DIC.RData`) is not distributed here. JADER preprocessing requires a dictionary obtained under the applicable MedDRA licence, saved as an R object named `DIC` with columns `pt` (Japanese PT) and `pt_name_en` (English PT).

## Code

`Code/` contains data preprocessing, disproportionality analyses, NMF pattern discovery, sensitivity analyses, transcriptomic analyses and figure/table scripts. Run scripts from the repository root after installing their required R packages and preparing the external inputs. Raw FAERS and JADER files belong in `Data/FAERS/` and `Data/JADER/`; TCGA/GTEx inputs belong in `Data/TCGA_GTEx/`, as specified in the scripts. Downstream scripts require outputs from earlier analyses; script numbers alone do not define execution order (for example, `6.2.15` requires outputs from `7.7`). See the manuscript and Supplementary Methods for analytical details.
