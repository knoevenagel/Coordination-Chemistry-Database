from __future__ import annotations

import importlib
import importlib.util
from pathlib import Path

import pytest


if importlib.util.find_spec("streamlit") is None:
    pytestmark = pytest.mark.skip(reason="streamlit is not installed")


def _load_module_from_file(module_name: str, path: Path) -> None:
    spec = importlib.util.spec_from_file_location(module_name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)


def test_import_app_modules() -> None:
    importlib.import_module("app.streamlit_ui.Home")
    importlib.import_module("app.streamlit_ui.state")
    importlib.import_module("app.streamlit_ui.components.sidebar")
    importlib.import_module("app.streamlit_ui.components.cards")
    importlib.import_module("app.streamlit_ui.components.tables")
    importlib.import_module("app.streamlit_ui.components.previews")


def test_import_page_files() -> None:
    pages_dir = Path(__file__).resolve().parents[2] / "app" / "streamlit_ui" / "pages"
    page_files = sorted(pages_dir.glob("[0-9][0-9]_*.py"))
    assert page_files, "No streamlit page files found"
    for idx, page_file in enumerate(page_files):
        _load_module_from_file(f"phase3a_page_{idx}", page_file)
