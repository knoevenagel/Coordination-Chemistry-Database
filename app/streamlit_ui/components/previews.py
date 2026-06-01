from __future__ import annotations

from typing import Any, Dict

import streamlit as st


def show_csv_preview(payload: Dict[str, Any], title: str) -> None:
    st.subheader(title)
    if not payload.get("ok"):
        st.warning(payload.get("error", "preview unavailable"))
        return
    st.dataframe(payload.get("data"), use_container_width=True)


def show_json_preview(payload: Dict[str, Any], title: str) -> None:
    st.subheader(title)
    if not payload.get("ok"):
        st.warning(payload.get("error", "preview unavailable"))
        return
    st.json(payload.get("data", {}))
