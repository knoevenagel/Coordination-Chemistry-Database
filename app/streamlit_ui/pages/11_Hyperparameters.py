from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import streamlit as st

from app.services import hyperparameter_manager
from app.storage import repositories
from app.streamlit_ui import state
from app.streamlit_ui.components.sidebar import render_sidebar
from app.streamlit_ui.components.tables import show_table

st.set_page_config(page_title="Hyperparameters - ChemDB Web", layout="wide")
render_sidebar()
st.title("Hyperparameter Library")

status = state.get_db_status()
with st.expander("Upload Hyperparameter YAML", expanded=False):
    with st.form("hp_upload_form"):
        hp_set_id = st.text_input("Hyperparameter Set ID", value="")
        hp_notes = st.text_area("Notes", value="")
        hp_file = st.file_uploader("YAML file", type=["yaml", "yml"])
        hp_submit = st.form_submit_button("Create Hyperparameter Version")
    if hp_submit:
        if not hp_set_id.strip() or hp_file is None:
            st.error("Hyperparameter set id and yaml file are required.")
        else:
            try:
                out = hyperparameter_manager.create_hyperparameter_version_from_upload(
                    workspace_root=status["workspace_root"],
                    hp_set_id=hp_set_id.strip(),
                    yaml_content=hp_file.getvalue(),
                    notes=hp_notes.strip(),
                )
            except Exception as exc:
                st.error(f"Upload failed: {exc}")
            else:
                st.success(f"Created {out['hp_set_id']}/{out['hp_version_id']}")
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
                with st.spinner("Rebuilding SQLite index..."):
                    proc = subprocess.run(cmd, capture_output=True, text=True)
                if proc.returncode == 0:
                    st.success("Rebuild index completed.")
                else:
                    st.error(f"Rebuild failed (exit={proc.returncode})")
                    st.code(proc.stderr or "(empty)", language="text")

conn = state.get_db_connection(read_only=True)
if conn is None:
    st.warning("SQLite DB not found.")
else:
    try:
        hp_sets = repositories.list_hyperparameter_sets(conn)
    finally:
        conn.close()

    show_table(hp_sets, empty_message="No hyperparameter sets found.")
    if hp_sets:
        selected_set = st.selectbox("Select Hyperparameter Set", [h["hp_set_id"] for h in hp_sets])
        conn = state.get_db_connection(read_only=True)
        if conn is not None:
            try:
                versions = repositories.list_hyperparameter_versions(conn, selected_set)
            finally:
                conn.close()
            show_table(versions, empty_message="No versions found.")
            if versions:
                selected_version = st.selectbox("Select Version", [v["hp_version_id"] for v in versions])
                conn = state.get_db_connection(read_only=True)
                if conn is not None:
                    try:
                        version = repositories.get_hyperparameter_version(conn, selected_set, selected_version)
                    finally:
                        conn.close()
                    yaml_path = Path(version["yaml_path"]) if version and version.get("yaml_path") else None
                    if yaml_path and yaml_path.is_file():
                        st.subheader("YAML Preview")
                        st.code(yaml_path.read_text(encoding="utf-8"), language="yaml")
