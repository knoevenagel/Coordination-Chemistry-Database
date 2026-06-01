from __future__ import annotations

from typing import Iterable

import pandas as pd
import streamlit as st


def _to_dataframe(rows: Iterable[dict]) -> pd.DataFrame:
    return pd.DataFrame(list(rows))


def show_table(rows: Iterable[dict], empty_message: str = "No data.") -> None:
    df = _to_dataframe(rows)
    if df.empty:
        st.info(empty_message)
        return
    path_cols = [c for c in df.columns if c.endswith("_path") or c.endswith("_root")]
    normal_cols = [c for c in df.columns if c not in path_cols]
    df = df[normal_cols + path_cols]
    st.dataframe(df, use_container_width=True)
