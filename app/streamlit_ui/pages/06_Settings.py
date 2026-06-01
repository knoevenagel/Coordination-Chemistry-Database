from __future__ import annotations

import platform
import subprocess
import sys

import streamlit as st

from app.storage import repositories
from app.streamlit_ui import state
from app.streamlit_ui.components.sidebar import render_sidebar

st.set_page_config(page_title="Settings - ChemDB Web", layout="wide")
render_sidebar()
st.title("Settings")

status = state.get_db_status()
st.markdown(f"**Workspace Root:** `{status['workspace_root']}`")
st.markdown(f"**DB Path:** `{status['db_path']}`")
st.markdown(f"**DB Exists:** `{status['db_exists']}`")
st.markdown(f"**Python:** `{sys.version.split()[0]}` (`{platform.python_implementation()}`)")

conn = state.get_db_connection(read_only=True)
if conn is None:
    st.warning("SQLite DB not found or not readable.")
else:
    try:
        summary = repositories.get_registry_summary(conn)
    finally:
        conn.close()
    st.subheader("SQLite Summary")
    st.json(summary)

st.subheader("Rebuild SQLite Index")
st.caption("Allowed write operation in Phase 3A: refresh registry by calling app.storage.cli rebuild-index.")

cmd = [
    sys.executable,
    "-m",
    "app.storage.cli",
    "rebuild-index",
    "--workspace-root",
    status["workspace_root"],
    "--db-path",
    status["db_path"],
]

st.code(" ".join(cmd), language="bash")

if st.button("Run Rebuild SQLite Index"):
    with st.spinner("Rebuilding SQLite index..."):
        proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode == 0:
        st.success("Rebuild finished successfully. Refresh page data if needed.")
    else:
        st.error(f"Rebuild failed (exit={proc.returncode})")
    st.subheader("stdout")
    st.code(proc.stdout or "(empty)", language="text")
    st.subheader("stderr")
    st.code(proc.stderr or "(empty)", language="text")
