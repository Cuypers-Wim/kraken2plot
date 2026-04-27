"""
loader.py — Read Kraken2 per-read output files into a tidy pandas DataFrame.

Kraken2 output columns
----------------------
0  status      C (classified) or U (unclassified)
1  contig_id   read / contig identifier
2  taxid       NCBI taxid assigned by Kraken2
3  seq_length  length of the sequence
4  kmer_hits   k-mer hit detail string (unused here)
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

from krakenplot.taxonomy import get_lineage

logger = logging.getLogger(__name__)

# Regex to pull a "sample name" from a filename.
# Tries to grab everything up to the first underscore; falls back to the stem.
_SAMPLE_RE = re.compile(r"^([^_]+)")


def _sample_name(path: Path) -> str:
    """Derive a short sample label from a file path."""
    # Remove '_kraken_output.txt' suffix if present
    name = path.stem
    if name.endswith('_kraken_output'):
        name = name[:-14]  # remove '_kraken_output'
    # Find the last part that starts with 'barcode' followed by digits, or 'unclassified'
    parts = name.split('_')
    for part in reversed(parts):
        if (part.startswith('barcode') and part[7:].isdigit()) or part == 'unclassified':
            return part
    # For regular names, return the full name without '_kraken_output'
    return name


def load_kraken_files(
    paths: List[Path],
    *,
    classified_only: bool = True,
) -> pd.DataFrame:
    """Read one or more Kraken2 output files into a single DataFrame.

    Parameters
    ----------
    paths:
        List of paths to Kraken2 per-read output files.
    classified_only:
        When ``True`` (default) only rows with ``status == "C"`` are kept.

    Returns
    -------
    pd.DataFrame
        Columns: ``sample``, ``status``, ``contig_id``, ``taxid``,
        ``seq_length``.

    Raises
    ------
    ValueError
        If *paths* is empty or none of the files contain any rows.
    """
    if not paths:
        raise ValueError("No input files provided.")

    frames: List[pd.DataFrame] = []
    for p in paths:
        logger.info("Reading %s …", p)
        df_tmp = pd.read_csv(
            p,
            sep="\t",
            header=None,
            usecols=[0, 1, 2, 3],
            names=["status", "contig_id", "taxid", "seq_length"],
            dtype={"status": str, "contig_id": str, "seq_length": str},
        )
        df_tmp.insert(0, "sample", _sample_name(p))
        if classified_only:
            df_tmp = df_tmp[df_tmp["status"] == "C"]
        df_tmp["taxid"] = pd.to_numeric(df_tmp["taxid"], errors="coerce")
        df_tmp["seq_length"] = pd.to_numeric(df_tmp["seq_length"], errors="coerce")
        frames.append(df_tmp)
        logger.info("  → %d classified reads", len(df_tmp))

    if not frames:
        raise ValueError("No data loaded from the provided files.")

    df = pd.concat(frames, ignore_index=True)
    logger.info("Total classified reads loaded: %d", len(df))
    return df


def add_lineage(
    df: pd.DataFrame,
    parent: Dict[int, int],
    rank: Dict[int, str],
    name: Dict[int, str],
) -> pd.DataFrame:
    """Append taxonomy lineage columns to *df* in-place (returns the same df).

    Parameters
    ----------
    df:
        DataFrame produced by :func:`load_kraken_files`.
    parent, rank, name:
        Taxonomy dictionaries from :func:`~krakenplot.taxonomy.load_taxonomy`.

    Returns
    -------
    pd.DataFrame
        The input DataFrame with additional columns:
        ``domain``, ``phylum``, ``family``, ``genus``, ``species``.
    """
    logger.info("Resolving lineages for %d unique taxids …", df["taxid"].nunique())
    unique_taxids: Dict[Optional[int], Dict[str, str]] = {}
    for tid in df["taxid"].unique():
        key = None if pd.isna(tid) else int(tid)
        unique_taxids[key] = get_lineage(key, parent, rank, name)

    def _lookup(tid):
        key = None if pd.isna(tid) else int(tid)
        return unique_taxids.get(key, {r: "unclassified" for r in ("domain", "phylum", "family", "genus", "species")})

    lineage_df = pd.DataFrame(list(df["taxid"].map(_lookup)))
    lineage_df.index = df.index
    df = pd.concat([df, lineage_df], axis=1)
    return df
