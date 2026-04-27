"""
taxonomy.py — Load and query NCBI taxonomy data (nodes.dmp / names.dmp).
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, Optional, Tuple

logger = logging.getLogger(__name__)

# Ranks that are extracted into lineage columns
LINEAGE_RANKS = ("domain", "phylum", "family", "genus", "species")


def load_taxonomy(
    nodes_file: str | Path,
    names_file: str | Path,
) -> Tuple[Dict[int, int], Dict[int, str], Dict[int, str]]:
    """Parse ``nodes.dmp`` and ``names.dmp`` into in-memory dictionaries.

    Parameters
    ----------
    nodes_file:
        Path to the NCBI ``nodes.dmp`` file.
    names_file:
        Path to the NCBI ``names.dmp`` file.

    Returns
    -------
    parent : dict[int, int]
        Maps each taxid to its parent taxid.
    rank : dict[int, str]
        Maps each taxid to its rank string (e.g. ``"species"``).
    name : dict[int, str]
        Maps each taxid to its scientific name.

    Raises
    ------
    FileNotFoundError
        If either file does not exist.
    """
    nodes_file = Path(nodes_file)
    names_file = Path(names_file)

    for p in (nodes_file, names_file):
        if not p.exists():
            raise FileNotFoundError(f"Taxonomy file not found: {p}")

    parent: Dict[int, int] = {}
    rank: Dict[int, str] = {}
    name: Dict[int, str] = {}

    logger.info("Loading taxonomy nodes from %s …", nodes_file)
    with nodes_file.open() as fh:
        for line in fh:
            parts = line.strip().split("\t|\t")
            if len(parts) >= 3:
                taxid = int(parts[0])
                parent[taxid] = int(parts[1])
                rank[taxid] = parts[2]

    logger.info("Loading taxonomy names from %s …", names_file)
    with names_file.open() as fh:
        for line in fh:
            parts = line.strip().split("\t|\t")
            if len(parts) >= 2 and "scientific name" in line:
                taxid = int(parts[0])
                name[taxid] = parts[1]

    logger.info("Loaded %d taxa, %d scientific names.", len(parent), len(name))
    return parent, rank, name


def get_lineage(
    taxid: Optional[int],
    parent: Dict[int, int],
    rank: Dict[int, str],
    name: Dict[int, str],
) -> Dict[str, str]:
    """Walk the taxonomy tree upward from *taxid* and collect standard ranks.

    Parameters
    ----------
    taxid:
        NCBI taxonomy identifier of the query taxon.  ``None`` or missing
        taxids are handled gracefully and return all ranks as ``"unclassified"``.
    parent:
        Parent-taxid mapping (from :func:`load_taxonomy`).
    rank:
        Rank mapping (from :func:`load_taxonomy`).
    name:
        Scientific-name mapping (from :func:`load_taxonomy`).

    Returns
    -------
    dict
        Keys are the ranks defined in :data:`LINEAGE_RANKS`; values are the
        scientific name at that rank, or ``"unclassified"`` if not found.
    """
    lineage: Dict[str, str] = {r: "unclassified" for r in LINEAGE_RANKS}

    if taxid is None or taxid not in parent:
        return lineage

    t = int(taxid)
    visited: set[int] = set()

    while t in parent and t not in visited:
        visited.add(t)
        r = rank.get(t, "")
        if r in lineage and lineage[r] == "unclassified":
            lineage[r] = name.get(t, str(t))
        if t == parent[t]:
            break
        t = parent[t]

    return lineage
