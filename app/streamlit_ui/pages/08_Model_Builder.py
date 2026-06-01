from __future__ import annotations

import streamlit as st

from app.services import model_build_manager
from app.storage import repositories
from app.streamlit_ui import state
from app.streamlit_ui.components.sidebar import render_sidebar

st.set_page_config(page_title="Model Builder - ChemDB Web", layout="wide")
render_sidebar()
st.title("Model Builder")
st.caption("Submit async run build pipeline job: init -> bind-ga -> core -> training -> rebuild-index")

status = state.get_db_status()
conn = state.get_db_connection(read_only=True)
projects = []
ga_sets = []
hp_sets = []
if conn is not None:
    try:
        projects = repositories.list_projects(conn)
        ga_sets = repositories.list_ga_sets(conn)
        hp_sets = repositories.list_hyperparameter_sets(conn)
    finally:
        conn.close()

with st.form("model_builder_form"):
    project_id = st.text_input("Project ID", value=(projects[0]["project_id"] if projects else "p001"))
    run_id = st.text_input("Run ID", value="")
    pubchem_source = st.text_input("PubChem Source (optional path in workspace)", value="")
    if ga_sets:
        ga_set_id = st.selectbox("GA Set", [g["ga_set_id"] for g in ga_sets])
    else:
        ga_set_id = st.text_input("GA Set", value="")
    versions = []
    if ga_set_id:
        conn = state.get_db_connection(read_only=True)
        if conn is not None:
            try:
                versions = repositories.list_ga_versions(conn, ga_set_id)
            finally:
                conn.close()
    if versions:
        ga_version_id = st.selectbox("GA Version", [v["ga_version_id"] for v in versions])
    else:
        ga_version_id = st.text_input("GA Version", value="")
    hp_set_options = ["(none)"] + [h["hp_set_id"] for h in hp_sets]
    hp_set_id = st.selectbox("Hyperparameter Set (optional)", hp_set_options)
    hp_versions = []
    if hp_set_id and hp_set_id != "(none)":
        conn = state.get_db_connection(read_only=True)
        if conn is not None:
            try:
                hp_versions = repositories.list_hyperparameter_versions(conn, hp_set_id)
            finally:
                conn.close()
    hp_version_id = st.selectbox(
        "Hyperparameter Version (optional)",
        ["(none)"] + [v["hp_version_id"] for v in hp_versions],
    )
    submitted = st.form_submit_button("Submit Model Build Job")

if submitted:
    if not run_id.strip():
        st.error("Run ID is required.")
    elif not ga_set_id or not ga_version_id:
        st.error("Please select GA set/version.")
    else:
        try:
            job = model_build_manager.submit_model_build_job(
                workspace_root=status["workspace_root"],
                db_path=status["db_path"],
                project_id=project_id.strip(),
                run_id=run_id.strip(),
                ga_set_id=ga_set_id,
                ga_version_id=ga_version_id,
                pubchem_source=(pubchem_source.strip() or None),
                hp_set_id=(hp_set_id if hp_set_id != "(none)" else None),
                hp_version_id=(hp_version_id if hp_version_id != "(none)" else None),
            )
        except Exception as exc:
            st.error(f"Submit failed: {exc}")
        else:
            st.success(f"Submitted job: {job['job_id']}")
            st.json(job)
            st.info("Track progress in `07_Jobs_and_Logs` page.")
