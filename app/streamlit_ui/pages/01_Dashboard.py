from __future__ import annotations

import streamlit as st

from app.storage import repositories
from app.streamlit_ui import state
from app.streamlit_ui.components.cards import render_metric_cards
from app.streamlit_ui.components.sidebar import render_sidebar
from app.streamlit_ui.components.tables import show_table

st.set_page_config(page_title="Dashboard - ChemDB Web", layout="wide")
render_sidebar()
st.title("Dashboard")

conn = state.get_db_connection(read_only=True)
if conn is None:
    st.warning("SQLite DB not found. Go to Settings and rebuild index.")
else:
    try:
        summary = repositories.get_registry_summary(conn)
        runs = repositories.list_runs(conn)
        models = repositories.list_models(conn)
        batches = repositories.list_all_batches(conn)
    finally:
        conn.close()

    render_metric_cards(
        {
            "projects": summary.get("projects", 0),
            "runs": summary.get("runs", 0),
            "ga_sets": summary.get("ga_sets", 0),
            "ga_versions": summary.get("ga_versions", 0),
            "models": summary.get("models", 0),
            "tasks": summary.get("tasks", 0),
            "batches": summary.get("batches", 0),
            "results": summary.get("results", 0),
        }
    )

    st.subheader("Recent Runs")
    show_table(runs[:20], empty_message="No runs found.")

    st.subheader("Recent Models")
    show_table(models[:20], empty_message="No models found.")

    st.subheader("Recent Batches")
    show_table(batches[:20], empty_message="No batches found.")
