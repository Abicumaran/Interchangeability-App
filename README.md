# PROXIMA Trueness + Bland–Altman — matched reference build

This is the matched V6.1 app/core reference used by the other app updates. Its statistical engine is unchanged.

- Existing global-flag handling and locked trueness/Bland–Altman methods are preserved.
- UI is simplified to one downloadable results file: the combined Excel workbook.
- Plot previews remain visible in the app, but separate PNG/CSV/ZIP download buttons are removed to maintain the one-results-file workflow.

Run with:

```bash
pip install -r requirements.txt
streamlit run app.py
```
