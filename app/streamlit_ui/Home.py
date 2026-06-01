from __future__ import annotations

import streamlit as st

from app.streamlit_ui import state
from app.streamlit_ui.components.sidebar import render_sidebar

st.set_page_config(page_title="ChemDB Web", layout="wide")

render_sidebar()

st.title("ChemDB Web")
st.caption("Phase 3B MVP dashboard (resource library + async execution)")

status = state.get_db_status()
if status["db_exists"]:
    st.success(f"DB ready: {status['db_path']}")
else:
    st.warning("SQLite DB not found. Open Settings page and run Rebuild SQLite Index.")

st.info("Use the pages in the sidebar to view Dashboard, GA Library, Task Library, Model Library, Results, and Settings.")
st.markdown(
    """
**Recommended flow**
1. Upload/prepare resources in `GA Library`, `Task Library`, `Hyperparameters`.
2. Start long tasks in `Model Builder`, `Model Selector`, `Recommendation`.
3. Track execution in `Jobs & Logs`.
4. Review outputs in `Model Selector Results` and download rankings.
"""
)
