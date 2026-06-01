from __future__ import annotations

import streamlit as st


def render_metric_cards(metrics: dict[str, int]) -> None:
    if not metrics:
        st.info("No metrics available.")
        return
    keys = list(metrics.keys())
    cols = st.columns(min(4, len(keys)))
    for idx, key in enumerate(keys):
        cols[idx % len(cols)].metric(key.replace("_", " ").title(), metrics.get(key, 0))
