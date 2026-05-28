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

from krakenplot.taxonomy import LINEAGE_RANKS, get_lineage

logger = logging.getLogger(__name__)
UNRESOLVED_LABEL = "unresolved"
UNCLASSIFIED_READ_LABEL = "Unclassified"

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
    lineage_cols = list(LINEAGE_RANKS)

    # Split Kraken2 unclassified reads from classified reads with unresolved ranks.
    # get_lineage uses "unclassified" as its generic fallback; this is the
    # point where we still know each read's Kraken2 status.
    classified_mask = df["status"] == "C"
    unclassified_mask = df["status"] == "U"
    lineage_df.loc[classified_mask, lineage_cols] = lineage_df.loc[
        classified_mask, lineage_cols
    ].replace("unclassified", UNRESOLVED_LABEL)
    lineage_df.loc[unclassified_mask, lineage_cols] = UNCLASSIFIED_READ_LABEL

    df = pd.concat([df, lineage_df], axis=1)
    return df


def taxon_count_table(df: pd.DataFrame, rank: str) -> pd.DataFrame:
    """Return read counts per sample and taxon for one lineage rank.

    Kraken2 unclassified read rows are included whenever they are present in
    *df*. In the normal CLI flow this means they appear only when the user
    passes ``--unclassified`` while loading Kraken2 output. Classified reads
    with missing lineage ranks are counted under ``"unresolved"``.
    """
    if rank not in LINEAGE_RANKS:
        raise ValueError(f"Rank must be one of: {', '.join(LINEAGE_RANKS)}")
    if rank not in df.columns:
        raise ValueError(f"Rank '{rank}' not found in DataFrame columns.")

    counts = (
        df.assign(**{rank: df[rank].fillna(UNRESOLVED_LABEL)})
        .groupby(["sample", rank])
        .size()
        .reset_index(name="count")
        .sort_values(["sample", "count", rank], ascending=[True, False, True])
        .reset_index(drop=True)
    )
    return counts


def export_taxon_count_tables(
    df: pd.DataFrame,
    output: Path,
    prefix: str,
    ranks=LINEAGE_RANKS,
) -> List[Path]:
    """Write one TSV count table per requested lineage rank."""
    output.mkdir(parents=True, exist_ok=True)
    written: List[Path] = []
    for rank in ranks:
        out_path = output / f"{prefix}_{rank}_counts.tsv"
        taxon_count_table(df, rank).to_csv(out_path, sep="\t", index=False)
        written.append(out_path)
    return written
