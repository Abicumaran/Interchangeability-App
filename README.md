# PROXIMA Trueness + Bland–Altman App

This Streamlit application is the app version of the final validated notebook
`PROXIMA_Trueness_and_BlandAltman_v2.ipynb`.

The App is available at: https://interchangeability-app.streamlit.app/

## Statistical pipeline

### 1. APT-style trueness regression

For every selected specimen group and analyte:

1. Parse donor, specimen type, and replicate from `bloodSampleId`, or use explicit columns selected in the app.
2. Retain only rows for which `global_flag == FALSE`.
3. Apply the validated generalized-ESD workflow to raw linked MHS/Sysmex replicate rows before donor averaging.
4. Remove detected rows for that analyte only; do not remove the donor or other analytes.
5. Average retained replicate measurements within donor.
6. Fit Huber regression with Sysmex reference on x and MHS on y.
7. Report Pearson correlation, p-value, and Fisher-z 95% CI; Huber slope, intercept, and equation; reference/MHS ranges; RMSE; bias; and optional Huber bootstrap CIs.
8. Save individual regression PNGs, a combined panel, CSV audits, and Excel sheets.

### 2. Paired-specimen Bland–Altman analysis

For every selected analyte:

1. Retain only `global_flag == FALSE` rows and apply analyte-specific flags when available.
2. Calculate raw-replicate residuals as `100 × (MHS − Sysmex) / Sysmex`.
3. Use Shapiro–Wilk to select the outlier branch:
   - normal residuals: two-sided Grubbs test with manually calculated sample-size-specific `Gcrit`;
   - non-normal residuals: MAD modified-Z threshold.
4. Remove at most one linked MHS/Sysmex replicate per analyte before donor averaging.
5. Compute donor means for paired specimen A and specimen B.
6. Compute all three methods:
   - **M02:** direct signed MHS A-versus-B comparison;
   - **M11:** signed Sysmex reference-adjusted comparison and the preferred reference-adjusted option;
   - **M05:** absolute-gap magnitude-only sensitivity analysis.
7. Report percentage and exact native-unit mean bias, 95% CI, normal-theory LoA, endpoint uncertainty, donor bootstrap CIs, empirical LoA, Shapiro results, normal-range midpoint and 50%-midpoint context, and provisional AC/CLIA screens.
8. Save individual and combined percentage/native Bland–Altman plots and donor-profile plots.


## Local run

```bash
python -m venv .venv

# Windows
.venv\Scripts\activate

# macOS/Linux
source .venv/bin/activate

pip install -r requirements.txt
streamlit run app.py
```

Windows users can also double-click `run_local.bat`; macOS/Linux users can run
`bash run_local.sh`.

## Streamlit Community Cloud

1. Extract this ZIP.
2. Put the files in the repository root.
3. Select `app.py` as the Streamlit main file.
4. Use Python 3.11 or 3.12.
5. Deploy or reboot the existing app.

## Important acceptance note

The default AC and CLIA values are editable contextual (prvisional) criteria for current reporting purposes. They are not automatically established
capillary-versus-venous acceptance criteria. The app asks the analyst to verify
them; otherwise every Pass/Fail result is labelled provisional/unverified.

## Privacy

Use de-identified research IDs only. Do not upload names, dates of birth,
medical-record numbers, addresses, or other direct identifiers.
