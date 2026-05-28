#!/usr/bin/env python3
"""
ChemDB CLI tool with subcommands.

Currently supported:
- did: calculate DID from SMILES

Usage examples:
  python tools/chemtool.py did --smiles "CCO"
  python tools/chemtool.py did --file smiles.txt
  echo "CCO" | python tools/chemtool.py did --stdin
  python tools/chemtool.py did --smiles "CCO" --json
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from typing import Iterable, Optional
from pathlib import Path

# Ensure we can import from project src path
PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_PATH = PROJECT_ROOT / "src"
if str(SRC_PATH) not in sys.path:
    sys.path.insert(0, str(SRC_PATH))

try:
    # Import the DID calculator from the research module
    from tools.DID_calculate import calculate_canonical_did  # type: ignore
except Exception as e:  # pragma: no cover
    print(f"Error: failed to import DID calculator: {e}", file=sys.stderr)
    sys.exit(2)


def _iter_smiles_from_stdin() -> Iterable[str]:
    for line in sys.stdin:
        smi = line.strip()
        if smi:
            yield smi


def _iter_smiles_from_file(path: Path) -> Iterable[str]:
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            smi = line.strip()
            if smi and not smi.startswith("#"):
                yield smi


def cmd_did(args: argparse.Namespace) -> int:
    # Determine input source
    inputs: Iterable[str]
    if args.smiles:
        inputs = [args.smiles]
    elif args.stdin:
        if sys.stdin.isatty():
            print("Error: --stdin specified but no piped input detected", file=sys.stderr)
            return 2
        inputs = _iter_smiles_from_stdin()
    elif args.file:
        file_path = Path(args.file)
        if not file_path.exists():
            print(f"Error: file not found: {file_path}", file=sys.stderr)
            return 2
        inputs = _iter_smiles_from_file(file_path)
    else:
        print("Error: one of --smiles, --file, or --stdin must be provided", file=sys.stderr)
        return 2

    exit_code = 0
    for smi in inputs:
        try:
            did = calculate_canonical_did(smi)
            if args.json:
                print(json.dumps({"smiles": smi, "did": did}, ensure_ascii=False))
            else:
                # default: TSV for easy piping
                if args.quiet:
                    print(did)
                else:
                    print(f"{smi}\t{did}")
        except Exception as e:
            exit_code = 1
            if args.json:
                print(json.dumps({"smiles": smi, "error": str(e)}, ensure_ascii=False))
            else:
                print(f"{smi}\tERROR: {e}", file=sys.stderr)
    return exit_code


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="chemtool", description="ChemDB CLI tool")
    subparsers = parser.add_subparsers(dest="subcommand", required=True)

    # did subcommand
    p_did = subparsers.add_parser("did", help="Calculate DID from SMILES")
    group = p_did.add_mutually_exclusive_group(required=True)
    group.add_argument("--smiles", "-s", type=str, help="SMILES string")
    group.add_argument("--file", "-f", type=str, help="Input text file, one SMILES per line")
    group.add_argument("--stdin", action="store_true", help="Read SMILES from stdin (one per line)")
    p_did.add_argument("--json", action="store_true", help="Output JSON lines {smiles, did}")
    p_did.add_argument("--quiet", "-q", action="store_true", help="Only print DID (no SMILES)")
    p_did.set_defaults(func=cmd_did)

    return parser


def main(argv: Optional[list[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if not hasattr(args, "func"):
        parser.print_help()
        return 2
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
