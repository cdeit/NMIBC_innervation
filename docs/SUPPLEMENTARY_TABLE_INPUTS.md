# Supplementary-table inputs

Updated 2026-09-22. Source workbooks are read without modification.

- `cfg$paths$table_s1`: `data/supp_tables/Table_S1_cohort_demographics.xlsx`, sheet `Table_S1_Cohort_demographics`.
- `cfg$paths$table_s4`: `data/supp_tables/Table_S4_segmented_PGP95_SYP_objects.xlsx`.
- `R/read_table_s1.R`: reads patient data and S4F image annotations. S1 source column names remain available; explicit aliases provide the existing analysis names. `Not sequenced` maps to missing Cohort, matching all 54 corresponding records in the earlier analysis table. No cohort filter is applied by the reader.
- `R/read_table_s4.R`: reads S4B–S4F, preserving the original stage-specific populations. It does not globally filter to final imaging or spatial flags.

Script 04 now uses S1. Its established helper still assigns within-cohort STaN groups from continuous abundance; the initial empty factor is only an input placeholder, not a patient classification. All relevant input values match the prior RDS. Counts remain 97 RNA-matched patients, A: 23 low/23 high and B: 10 low/10 high.

Script 08 now uses S4F image-to-patient mapping and spatial inclusion flags. It reads the separate spatial config from its current location, `config/spatial_config.R`. Its patient-level clinical values and stored STaN groups still come from the configured RDS. The same 45 adequate-BCG spatial images and their stored strata are retained. No preprocessing/measurement kernels were changed.

S4 reader validation returned 97,568 preliminary objects, 53,979 candidate/reference objects, 53,979 quantification rows, 32,501 tissue-region rows and 778 image annotations. Script 01 still requires its historical matched/compact clinical tables; S1 is not silently substituted for those upstream inputs. Full phenotyping was not rerun.

## Pending clinical input decision

Scripts 02, 03 and 06 retain their existing RDS input until the remaining publication-table fields are available. S1 omits stored STaN groups and some molecular subtype fields used by the existing figures. It also changes LVI representation: the nine `NA; TaHG` rows correspond to three No and six missing values in the existing RDS. These cannot be translated back using a single label mapping. The RDS provides these supplemental values; replacing it requires equivalent fields in the publication tables. No LVI definition, missingness, subtype or stratum has been inferred.

Validation covers input values, cohort counts and spatial membership; statistical models and figures were not rerun for this input-routing change.
