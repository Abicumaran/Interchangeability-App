# PROXIMA Trueness + Bland–Altman — final automated single-output build

https://interchangeability-app.streamlit.app/

This build keeps the validated trueness and Bland–Altman engine while simplifying the reportable workflow.

- Select MHS model columns independently (e.g. PLT_2, PLT_3, RDW_2); model variants share a genuine reference column when appropriate and generate separate result rows.
- Global-flag and analyte-specific flag exclusions are applied before reportable analysis and retained in dedicated audit sheets.
- Validated replicate-level outlier handling is always enabled for the reportable run: trueness uses generalized ESD; Bland–Altman uses Shapiro-Wilk → Grubbs when normal or robust MAD when non-normal.
- Trueness remains Huber regression with Pearson/Fisher-z reporting.
- For paired-specimen Bland–Altman, the single reportable method is selected automatically as M11 (preferred signed Sysmex-reference-adjusted comparison). M02/M05 may be evaluated internally but are not exported as competing reportable outcomes.
- For single-specimen MHS-versus-reference Bland–Altman, the reportable method is REF.
- The only downloadable results file is one Excel workbook with: `Trueness Results`, `Bland-Altman Results`, `Outliers`, `global flag TRUE`, `analyte flag TRUE`.

Run with:

```bash
pip install -r requirements.txt
streamlit run app.py
```

## September 17, 2026: deployment fixes (V6.3)

Replace **app.py, analysis_core.py, and requirements.txt together**, then reboot/redeploy
Streamlit so requirements are installed. The requirements now include `xlsxwriter>=3.2`.
If the package is still missing on a cached installation, the engine transparently
writes the same five Excel report tabs with the already-required `openpyxl` engine;
the main analysis is not interrupted by a missing optional formatting library.

The donor review control is only mandatory when the **selected identity rule** can
actually merge/split donors ambiguously. A complete and internally consistent
explicit `Donor` column resolves repeated source tokens such as `D03` and `D03b`
automatically; the identity audit remains visible. An invalid/incomplete donor
column is never silently replaced with parsed tokens.

**Input-data caveat:** `P-007_Capillary blood 1_PLT3_again_for analysis_2.xlsx`
contains MHS analytes (including `PLT_3`) and a usable `Donor` column, but **no
matched Sysmex/reference analyte columns** (such as `PLT_ref`). A valid trueness
or reference-adjusted Bland–Altman analysis of those data requires a separately
measured, matched reference dataset joined into the input workbook. The app
warns explicitly and never fabricates reference values or treats device values
as their own reference. After merging the genuine reference data, select `PLT_3` directly to obtain its own results, independently of `PLT` or `PLT_2`.


## September 17, 2026: final two Excel-output fixes (V6.4)

1. The `global flag TRUE` worksheet now lists the actual selected input rows with
   `global_flag = TRUE`, with an explicit exclusion reason. It is populated even
   when Bland–Altman is disabled or when flags occur outside the BA arm.
   No flagged rows are included in the reportable analyses.
2. The `analyte` field in both `Trueness Results` and `Bland-Altman Results`
   reports the **selected source MHS model column**, e.g. `PLT_2`, `PLT_3`,
   `RDW_2`. Model variants are analyzed independently, can share the same genuine
   reference column, and are never merged into one generic `PLT` output row.

All other statistical analysis branches and reportable workbook tabs are
unchanged. Deploy `app.py` and `analysis_core.py` together (the supplied
`requirements.txt` is unchanged from V6.3).
