# Deploy to Streamlit Community Cloud

## Important: do not leave the app inside a ZIP in GitHub

Streamlit does not unzip an uploaded archive. Extract the downloaded package first, then upload or commit the individual files so the GitHub repository root contains:

```text
app.py
analysis_core.py
requirements.txt
README.md
.streamlit/config.toml
```

The ZIP may remain outside the repository, but it is not the deployed app.

## Streamlit settings

- Repository: your GitHub repository
- Branch: `main`
- Main file path: `app.py`
- Python: choose Python 3.12 in Advanced settings for the most conservative deployment environment

After replacing the files, reboot the app. If changing the Python version, delete and redeploy the Streamlit app because Community Cloud does not change Python in place.
