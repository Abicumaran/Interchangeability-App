# PROXIMA Trueness + Bland–Altman — final automated single-output build

https://interchangeability-app.streamlit.app/

This build keeps the validated trueness and Bland–Altman engine while simplifying the reportable workflow.

- Select exactly which supported analytes to run, displayed in the same suggested MHS-column order as the Short-Term app.
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
