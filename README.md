
# Interchangeability App

A Streamlit application for paired clinical-method interchangeability analysis.
https://interchangeability-app.streamlit.app/

## Main analysis

- Upload CSV/XLSX data.
- Supports:
  - **Long format**: donor/subject + method `Type` + analyte columns.
  - **Wide format**: explicit Method A and Method B columns.
- Donor/subject pairing and replicate averaging.
- Optional normal-range filtering.
- Bland–Altman mean percent bias:
  - `100 × (Method A − Method B) / Method B`
  - 95% confidence interval for mean bias
  - limits of agreement
- P/F against an analyte-specific allowable trueness-bias interval.
- Built-in trueness profiles:
  - MHS screenshot profile
  - Initial EP09 profile
  - custom values
  - uploaded bias table
- Normality diagnostics:
  - Shapiro–Wilk on paired differences and percent bias
  - Levene variance check
  - method recommendation
- Optional paired t-test or Wilcoxon signed-rank test with Holm correction.
- Optional correction/sensitivity analyses:
  - direct log ratio
  - leave-one-out median log-ratio recalibration
  - donor-specific Sysmex absolute difference-in-differences
  - donor-specific Sysmex log-ratio correction
- Optional outlier sensitivity:
  - manual Grubbs critical value
  - robust MAD
  - generalized ESD
  - automatic normality-guided recommendation
- Optional sample-size estimates:
  - mean-bias equivalence power
  - mean-bias CI precision
  - Bland–Altman LoA precision
- Optional SD-distance coverage table.
- Downloadable Excel, CSV, PNG plots, and ZIP package.

## Default analytes

RBC, WBC, PLT, HCT, HGB, MCV, RDW, MCH, MCHC, NEUT, LYMPH, and MXD.

The app automatically recognizes the column names used in the EP09/APT analyses, including `WBC_2`, `NEUT_2`, `LYMPH_2`, `MXD_2`, and corresponding `_ref` columns.

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

## Deploy on Streamlit Community Cloud

1. Create a public or private GitHub repository.
2. Upload all files in this folder.
3. In Streamlit Community Cloud, choose **Create app**.
4. Select the repository and branch.
5. Set the main file path to `app.py`.
6. Deploy.

## Statistical interpretation

The primary P/F rule is:

- **PASS** only when the entire mean BA percent-bias 95% CI lies inside the selected `± allowable bias` interval.

CI overlap, paired-test p-values, and LoA are useful supporting statistics but are not substitutes for the pre-specified P/F rule.

Outlier removal and recalibration should be treated as sensitivity analyses unless independently justified and validated.

## Privacy

Do not upload names, medical record numbers, dates of birth, or other direct identifiers. Use de-identified study IDs only.

