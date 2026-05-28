"""
plots.py — Generate publication-quality taxonomy composition plots from a
classified-reads DataFrame.

All public functions accept a *pandas* DataFrame (as returned by
:func:`~krakenplot.loader.add_lineage`) and return a
:class:`matplotlib.figure.Figure` object.  Saving is handled by the CLI.
"""

from __future__ import annotations

import logging
from typing import List, Optional

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import numpy as np

matplotlib.use("Agg")  # non-interactive backend; overridden when --show is used

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

_CMAP = "tab20"
_ALPHA = 0.92
UNRESOLVED_LABEL = "unresolved"
UNCLASSIFIED_READ_LABEL = "Unclassified"
LEGACY_UNCLASSIFIED_LABEL = "unclassified"


def _top_n_and_other(
    series: pd.Series,
    n: int,
    other_label: str = "Other",
) -> pd.Series:
    """Keep the *n* largest values; fold the rest into *other_label*."""
    top = series.nlargest(n).index
    always_keep = {UNRESOLVED_LABEL, UNCLASSIFIED_READ_LABEL, LEGACY_UNCLASSIFIED_LABEL}
    return pd.Series(
        series.index.map(lambda x: x if x in top or x in always_keep else other_label),
        index=series.index,
    )


def _stacked_bar(
    pivot: pd.DataFrame,
    *,
    ylabel: str,
    title: str,
    figsize: tuple = (14, 7),
    legend_title: str = "",
    pct: bool = False,
    palette: str = "tab20",
) -> plt.Figure:
    """Shared helper that draws a stacked bar chart from a wide-format pivot."""
    if pct:
        data = pivot.div(pivot.sum(axis=1), axis=0).mul(100)
        ylabel = ylabel or "Percent of reads (%)"
    else:
        data = pivot

    n_cols = len(data.columns)
    cmap = plt.get_cmap(palette, max(n_cols, 1))
    colors = [cmap(i) for i in range(n_cols)]

    # Keep the two "not resolved to this rank" cases visually distinct.
    if UNCLASSIFIED_READ_LABEL in data.columns:
        idx = data.columns.get_loc(UNCLASSIFIED_READ_LABEL)
        colors[idx] = 'gray'
    if UNRESOLVED_LABEL in data.columns:
        idx = data.columns.get_loc(UNRESOLVED_LABEL)
        colors[idx] = 'lightgray'

    fig, ax = plt.subplots(figsize=figsize)
    data.plot(
        kind="bar",
        stacked=True,
        color=colors,
        width=0.75,
        edgecolor="none",
        alpha=_ALPHA,
        ax=ax,
    )
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_xlabel("Sample", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold", pad=12)
    ax.tick_params(axis="x", rotation=45)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(
        handles,
        labels,
        title=legend_title,
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
        fontsize=9,
    )
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Public plot functions
# ---------------------------------------------------------------------------


def plot_phylum_composition(
    df: pd.DataFrame,
    *,
    rank: str = "phylum",
    top_n: int = 10,
    pct: bool = True,
    figsize: tuple = (14, 7),
    title: str = None,
    palette: str = "tab20",
) -> plt.Figure:
    """Stacked bar chart of taxonomic composition across samples at specified rank.

    Parameters
    ----------
    df:
        DataFrame with lineage columns.
    rank:
        Taxonomic rank to plot ("domain", "phylum", "family", "genus", "species").
    top_n:
        Number of top taxa to show individually; the rest are merged into
        an *"Other"* category. Not used for "domain".
    pct:
        If ``True`` (default) plot percentages; otherwise raw read counts.
    figsize:
        Matplotlib figure size ``(width, height)`` in inches.
    title:
        Plot title. If None, auto-generated.
    palette:
        Color palette name.

    Returns
    -------
    matplotlib.figure.Figure
    """
    if title is None:
        title = f"{rank.capitalize()} composition across samples"

    counts = (
        df.assign(**{rank: df[rank].fillna(UNRESOLVED_LABEL)})
        .groupby(["sample", rank])
        .size()
        .reset_index(name="count")
    )

    if rank == "domain":
        # No top_n for domain
        pivot = counts.pivot(index="sample", columns=rank, values="count").fillna(0)
        pivot = pivot.loc[sorted(pivot.index)]
        pivot = pivot[pivot.sum().sort_values(ascending=False).index]
        legend_title = "Domain"
    else:
        global_totals = counts.groupby(rank)["count"].sum()
        grouped = _top_n_and_other(global_totals, top_n)
        counts[f"{rank}_grouped"] = counts[rank].map(grouped)

        agg = (
            counts.groupby(["sample", f"{rank}_grouped"])["count"]
            .sum()
            .reset_index()
        )
        pivot = agg.pivot(index="sample", columns=f"{rank}_grouped", values="count").fillna(0)
        pivot = pivot.loc[sorted(pivot.index)]
        pivot = pivot[pivot.sum().sort_values(ascending=False).index]
        legend_title = rank.capitalize()

    ylabel = "Percent of reads (%)" if pct else "Read count"
    return _stacked_bar(
        pivot, ylabel=ylabel, title=title, figsize=figsize,
        legend_title=legend_title, pct=pct, palette=palette,
    )


def plot_domain_composition(
    df: pd.DataFrame,
    *,
    pct: bool = True,
    figsize: tuple = (14, 7),
    title: str = None,
    palette: str = "tab20",
) -> plt.Figure:
    """Stacked bar chart of domain-level (superkingdom) composition.

    Parameters
    ----------
    df:
        DataFrame with lineage columns.
    pct:
        If ``True`` (default) plot percentages; otherwise raw counts.
    figsize:
        Matplotlib figure size.
    title:
        Plot title. If None, auto-generated.
    palette:
        Color palette name.

    Returns
    -------
    matplotlib.figure.Figure
    """
    return plot_phylum_composition(
        df, rank="domain", pct=pct, figsize=figsize, title=title, palette=palette
    )


def plot_phylum_breakdown(
    df: pd.DataFrame,
    phyla: List[str],
    *,
    top_n: int = 20,
    figsize_per_phylum: tuple = (14, 4),
) -> plt.Figure:
    """Per-phylum breakdown showing species / genus / family composition.

    For each phylum in *phyla*, a separate subplot is drawn.  The best
    available classification (species > genus > family) is used as the label.

    Parameters
    ----------
    df:
        DataFrame with lineage columns.
    phyla:
        List of phylum names to plot.
    top_n:
        Maximum number of distinct taxa to show per phylum; the rest are
        grouped into *"Other"*.
    figsize_per_phylum:
        ``(width, height)`` per subplot row.

    Returns
    -------
    matplotlib.figure.Figure
    """
    phyla_lower = [p.lower() for p in phyla]
    subset = df[df["phylum"].str.lower().isin(phyla_lower)].copy()

    if subset.empty:
        logger.warning("No reads found for the selected phyla: %s", phyla)
        fig, ax = plt.subplots(figsize=(8, 2))
        ax.text(0.5, 0.5, f"No data found for: {', '.join(phyla)}",
                ha="center", va="center", fontsize=12)
        ax.axis("off")
        return fig

    subset["best_label"] = subset.apply(
        lambda row: (
            row["species"] if row["species"] != UNRESOLVED_LABEL
            else (row["genus"] if row["genus"] != UNRESOLVED_LABEL else row["family"])
        ),
        axis=1,
    )

    # Build counts per (phylum, sample, best_label)
    counts = (
        subset.groupby(["phylum", "sample", "best_label"])
        .size()
        .reset_index(name="count")
    )

    # Keep top_n per phylum, group rest into 'Other'
    top_per_phylum = (
        counts.groupby(["phylum", "best_label"])["count"]
        .sum()
        .reset_index()
        .sort_values(["phylum", "count"], ascending=[True, False])
        .groupby("phylum")
        .head(top_n)
    )
    top_labels = set(top_per_phylum["best_label"].unique())
    counts["label"] = counts["best_label"].where(
        counts["best_label"].isin(top_labels), "Other"
    )
    counts = counts.groupby(["phylum", "sample", "label"])["count"].sum().reset_index()

    # Normalise to proportions per sample within each phylum
    sample_totals = counts.groupby(["phylum", "sample"])["count"].transform("sum")
    counts["proportion"] = counts["count"] / sample_totals

    # Only plot phyla that actually have data
    phyla_present = [p for p in phyla if p in counts["phylum"].unique()]
    if not phyla_present:
        logger.warning("None of the requested phyla were found after filtering.")
        fig, ax = plt.subplots(figsize=(8, 2))
        ax.text(0.5, 0.5, "No data found after filtering.",
                ha="center", va="center", fontsize=12)
        ax.axis("off")
        return fig

    n = len(phyla_present)
    fig, axes = plt.subplots(
        n, 1,
        figsize=(figsize_per_phylum[0], figsize_per_phylum[1] * n),
        sharex=True,
    )
    if n == 1:
        axes = [axes]

    for ax, phylum in zip(axes, phyla_present):
        data = (
            counts[counts["phylum"] == phylum]
            .pivot(index="sample", columns="label", values="proportion")
            .fillna(0)
        )
        data = data.loc[sorted(data.index)]

        n_cols = len(data.columns)
        cmap = plt.get_cmap(_CMAP, max(n_cols, 1))
        colors = [cmap(i) for i in range(n_cols)]

        data.plot(
            kind="bar",
            stacked=True,
            color=colors,
            width=0.75,
            edgecolor="none",
            alpha=_ALPHA,
            ax=ax,
        )
        ax.set_title(f"{phylum} — species / genus / family composition",
                     fontsize=12, fontweight="bold")
        ax.set_ylabel("Proportion within phylum", fontsize=10)
        ax.set_ylim(0, 1.05)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.legend(
            title="Top taxa",
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            frameon=False,
            fontsize=8,
        )

    axes[-1].set_xlabel("Sample", fontsize=12)
    for ax in axes:
        ax.tick_params(axis="x", rotation=45)

    fig.suptitle(
        "Intra-phylum composition across samples",
        fontsize=14, fontweight="bold", y=1.01,
    )
    fig.tight_layout()
    return fig


def plot_species_heatmap(
    df: pd.DataFrame,
    *,
    rank: str = "species",
    top_n: int = 30,
    figsize: tuple = (12, 10),
) -> plt.Figure:
    """Heatmap of relative abundance at the given *rank* across samples.

    Parameters
    ----------
    df:
        DataFrame with lineage columns.
    rank:
        Taxonomic rank to aggregate (default ``"species"``).
    top_n:
        Number of top taxa (by total reads) to include.
    figsize:
        Matplotlib figure size.
    include_unclassified:
        Include Kraken2 ``U`` reads labelled ``"Unclassified"`` when they are
        present. Classified reads with missing ranks remain labelled
        ``"unresolved"``.

    Returns
    -------
    matplotlib.figure.Figure
    """
    if rank not in df.columns:
        raise ValueError(f"Rank '{rank}' not found in DataFrame columns.")

    if include_unclassified:
        plot_df = df
    else:
        plot_df = df[~df[rank].isin([UNCLASSIFIED_READ_LABEL, LEGACY_UNCLASSIFIED_LABEL])]
    counts = (
        plot_df.groupby(["sample", rank])
        .size()
        .reset_index(name="count")
    )
    top_taxa = (
        counts.groupby(rank)["count"].sum().nlargest(top_n).index
    )
    counts = counts[counts[rank].isin(top_taxa)]

    pivot = counts.pivot(index=rank, columns="sample", values="count").fillna(0)
    # Normalise each sample to relative abundance
    pivot = pivot.div(pivot.sum(axis=0), axis=1)

    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(pivot.values, aspect="auto", cmap="YlOrRd", interpolation="nearest")
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=45, ha="right", fontsize=10)
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index, fontsize=8)
    ax.set_title(
        f"Top {top_n} {rank}s — relative abundance heatmap",
        fontsize=13, fontweight="bold", pad=10,
    )
    fig.colorbar(im, ax=ax, label="Relative abundance")
    fig.tight_layout()
    return fig


def plot_read_counts(
    df: pd.DataFrame,
    *,
    figsize: tuple = (10, 5),
    include_unclassified: bool = False,
) -> plt.Figure:
    """Simple bar chart of total loaded read counts per sample.

    Parameters
    ----------
    df:
        DataFrame with a ``sample`` column.
    figsize:
        Matplotlib figure size.
    Returns
    -------
    matplotlib.figure.Figure
    """
    counts = df.groupby("sample").size().sort_index()
    label = "Reads"
    fig, ax = plt.subplots(figsize=figsize)
    bars = ax.bar(counts.index, counts.values, color="#4C72B0", edgecolor="none",
                  alpha=0.85)
    ax.bar_label(bars, fmt="%d", padding=3, fontsize=9)
    ax.set_xlabel("Sample", fontsize=12)
    ax.set_ylabel(label, fontsize=12)
    ax.set_title(f"{label} per sample", fontsize=13, fontweight="bold")
    ax.tick_params(axis="x", rotation=45)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    return fig
