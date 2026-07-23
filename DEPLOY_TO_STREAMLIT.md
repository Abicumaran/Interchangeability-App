# Deploy to Streamlit Community Cloud

Extract the ZIP before uploading it to GitHub. The repository root should contain:

```text
app.py
analysis_core.py
requirements.txt
README.md
.streamlit/config.toml
```

In Streamlit Community Cloud:

- choose the repository and branch;
- set the main file path to `app.py`;
- select Python 3.11 or 3.12;
- deploy, or reboot the existing app after replacing the files.

The app writes each run to a temporary server directory and places the complete
results in memory for download. Uploaded data are not bundled with the source code.
