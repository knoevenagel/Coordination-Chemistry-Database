#!/usr/bin/env python3
"""Z-score normalize element_features.csv; write element_features_zscore.csv."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
DEFAULT_IN = ROOT / "element_features.csv"
DEFAULT_OUT = ROOT / "element_features_zscore.csv"


def zscore_features(in_path: Path, out_path: Path) -> None:
    df = pd.read_csv(in_path)
    exclude = {"element", "symbol"}
    numeric_cols = [
        c for c in df.columns
        if c not in exclude and pd.api.types.is_numeric_dtype(df[c])
    ]
    for c in df.columns:
        if c in exclude:
            continue
        if c not in numeric_cols:
            df[c] = pd.to_numeric(df[c], errors="coerce")
            if df[c].dtype in ("float64", "int64"):
                numeric_cols.append(c)
    numeric_cols = [c for c in numeric_cols if c in df.columns]
    out = df.copy()
    for c in numeric_cols:
        mu = out[c].mean()
        sigma = out[c].std()
        if sigma == 0 or pd.isna(sigma):
            out[c] = 0.0
        else:
            out[c] = (out[c] - mu) / sigma
    out[numeric_cols] = out[numeric_cols].fillna(0)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False)
    print(f"Written {out_path} ({len(out)} rows, {len(numeric_cols)} z-scored columns)")


def main() -> None:
    parser = argparse.ArgumentParser(description="Z-score metal element features CSV.")
    parser.add_argument("--input", type=str, default=None, help="Input element_features.csv")
    parser.add_argument("--output", type=str, default=None, help="Output element_features_zscore.csv")
    args = parser.parse_args()
    in_path = Path(args.input).resolve() if args.input else DEFAULT_IN
    out_path = Path(args.output).resolve() if args.output else DEFAULT_OUT
    if not in_path.is_file():
        raise FileNotFoundError(f"Input not found: {in_path}")
    zscore_features(in_path, out_path)


if __name__ == "__main__":
    main()
