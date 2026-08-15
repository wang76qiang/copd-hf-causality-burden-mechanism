# MedComm reproducibility archive update addendum

Date: 2026-08-15

This folder contains the files that must be uploaded to GitHub and used to create a new Zenodo version before submission.

Required public-archive changes:

1. Remove the previously listed fourth author from README, CITATION.cff, and Zenodo creators.
2. Replace “Submitted to EClinicalMedicine (2026)” with the MedComm submission status; the earlier submission has formally ended.
3. Replace the stale 210-country wording with 204 countries and territories.
4. Replace Tables S1, S2, S3, S6, and S7 with the corrected CSVs in `tables/`; add S9 and S10.
5. Replace Figures 5–7 with the corrected PNG/PDF files in `figures/` and upload the revised plotting scripts.
6. Create a new Zenodo version; do not silently overwrite the published v1.0.2 record.

Known limitation still requiring source-data access: SII/CI and frontier re-execution require the GHDx-restricted GBD 2021 SDI file.
