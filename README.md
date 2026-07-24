# PROXIMA Interchangeability App 
https://interchangeability-app.streamlit.app/
## Required repository files

Main Files:

- `app.py`
- `analysis_core.py`
- `requirements.txt`

Set the Streamlit Cloud main file path to `app.py`.

## Replace the deployed files

1. Delete or overwrite the old `app.py` and `analysis_core.py`.
2. Upload all files from this package without renaming them.
3. Commit the changes.
4. In Streamlit Cloud, open **Manage app** and reboot the app.
5. If the old build remains cached, use **Reboot app** once more after the new
   commit is visible.

The fixed app checks the core API version at startup. A future mixed-file
upgrade will now display a clear version-mismatch message instead of the
redacted ImportError page.

## Local launch

```bash
pip install -r requirements.txt
streamlit run app.py
```
