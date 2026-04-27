"""
cli.py — Command-line interface for KrakenPlot.

Usage examples
--------------
# Phylum composition across all *.txt files in a directory
kraken2plot phylum --nodes nodes.dmp --names names.dmp kraken2/

# Phylum-level breakdown for specific phyla
kraken2plot breakdown --nodes nodes.dmp --names names.dmp \\
    --phyla Pseudomonadota Bacteroidota \\
    kraken2/*.txt

# Heatmap of top-30 species
kraken2plot heatmap --nodes nodes.dmp --names names.dmp \\
    --rank species --top-n 30 kraken2/

# All plots at once
kraken2plot all --nodes nodes.dmp --names names.dmp kraken2/
"""

from __future__ import annotations

import glob
import logging
import sys
from pathlib import Path
from typing import List, Optional, Tuple

import click
import matplotlib
import matplotlib.pyplot as plt

from krakenplot import __version__
from krakenplot.loader import add_lineage, load_kraken_files
from krakenplot.plots import (
    plot_domain_composition,
    plot_phylum_breakdown,
    plot_phylum_composition,
    plot_read_counts,
    plot_species_heatmap,
)
from krakenplot.taxonomy import load_taxonomy

# ---------------------------------------------------------------------------
# Logging setup
# ---------------------------------------------------------------------------

logging.basicConfig(
    format="[%(levelname)s] %(message)s",
    level=logging.INFO,
    stream=sys.stderr,
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Shared options / helpers
# ---------------------------------------------------------------------------

def _resolve_inputs(inputs: Tuple[str, ...]) -> List[Path]:
    """Expand globs, directories (*.txt), and explicit files into Path list."""
    paths: List[Path] = []
    for inp in inputs:
        p = Path(inp)
        if p.is_dir():
            found = sorted(p.glob("*_output.txt"))
            if not found:
                # Fall back: any .txt in the directory
                found = sorted(p.glob("*.txt"))
            paths.extend(found)
        else:
            paths.extend(Path(g) for g in sorted(glob.glob(inp)) if Path(g).is_file())

    if not paths:
        click.echo("[ERROR] No input files found. Check your paths.", err=True)
        sys.exit(1)

    click.echo(f"[INFO] Found {len(paths)} input file(s):", err=True)
    for p in paths:
        click.echo(f"       {p}", err=True)
    return paths


def _save_or_show(fig: plt.Figure, output: Optional[Path], prefix: str, suffix: str, show: bool, dpi: int):
    """Save *fig* to *output* directory or display it."""
    if show:
        matplotlib.use("TkAgg")
        plt.show()
        return
    out_dir = Path(output) if output else Path(".")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{prefix}_{suffix}.png"
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    click.echo(f"Saved → {out_path}")
    plt.close(fig)


# Common options shared by all sub-commands
_common_options = [
    click.argument("inputs", nargs=-1, required=True, metavar="INPUT…"),
    click.option("--nodes", required=True, type=click.Path(exists=True),
                 help="Path to NCBI nodes.dmp"),
    click.option("--names", required=True, type=click.Path(exists=True),
                 help="Path to NCBI names.dmp"),
    click.option("--output", "-o", default=None, type=click.Path(),
                 help="Output directory (default: current directory)"),
    click.option("--prefix", default="kraken2plot",
                 help="File name prefix (default: kraken2plot)"),
    click.option("--dpi", default=150, show_default=True,
                 help="Figure resolution in dots per inch"),
    click.option("--show", is_flag=True,
                 help="Display interactive plot instead of saving"),
    click.option("--unclassified/--no-unclassified", default=False,
                 help="Include unclassified reads (default: exclude)"),
    click.option("--title", default=None,
                 help="Plot title (default: auto-generated)"),
    click.option("--palette", default="tab20", show_default=True,
                 type=click.Choice(["tab20", "viridis", "plasma", "Set1", "Set2", "tab10", "Paired"]),
                 help="Color palette for the plot"),
]


def add_common_options(func):
    for option in reversed(_common_options):
        func = option(func)
    return func


# ---------------------------------------------------------------------------
# CLI group
# ---------------------------------------------------------------------------

@click.group()
@click.version_option(__version__, prog_name="kraken2plot")
def cli():
    """KrakenPlot — visualise Kraken2 taxonomic classification output.

    Provide one or more Kraken2 *output* files (tab-separated per-read files)
    or a directory containing them, together with NCBI taxonomy dump files
    (nodes.dmp and names.dmp), and choose a plot type.

    \b
    Quick start
    -----------
    krakenplot phylum --nodes nodes.dmp --names names.dmp kraken2/
    krakenplot all    --nodes nodes.dmp --names names.dmp kraken2/
    """


# ---------------------------------------------------------------------------
# phylum sub-command
# ---------------------------------------------------------------------------

@cli.command()
@add_common_options
@click.option("--rank", default="phylum", show_default=True,
              type=click.Choice(["domain", "phylum", "family", "genus", "species", "all"]),
              help="Taxonomic rank to plot, or 'all' for all ranks")
@click.option("--top-n", default=10, show_default=True,
              help="Number of top taxa to show individually (not used for domain)")
@click.option("--counts", "show_counts", is_flag=True,
              help="Also save a read-count version of the plot")
def phylum(inputs, nodes, names, output, prefix, dpi, show, unclassified, title, palette, rank, top_n, show_counts):
    """Stacked bar chart of taxonomic composition across samples."""
    paths = _resolve_inputs(inputs)
    par, rnk, nam = load_taxonomy(nodes, names)
    df = load_kraken_files(paths, classified_only=not unclassified)
    df = add_lineage(df, par, rnk, nam)

    out = Path(output) if output else None

    ranks = ["domain", "phylum", "family", "genus", "species"] if rank == "all" else [rank]

    for r in ranks:
        plot_title = title or f"{r.capitalize()} composition across samples"
        fig = plot_phylum_composition(df, rank=r, top_n=top_n, pct=True, title=plot_title, palette=palette)
        suffix = f"{r}_pct"
        _save_or_show(fig, out, prefix, suffix, show, dpi)

        if show_counts:
            fig2 = plot_phylum_composition(df, rank=r, top_n=top_n, pct=False, title=plot_title, palette=palette)
            suffix2 = f"{r}_counts"
            _save_or_show(fig2, out, prefix, suffix2, show, dpi)


# ---------------------------------------------------------------------------
# domain sub-command
# ---------------------------------------------------------------------------

@cli.command()
@add_common_options
@click.option("--rank", default="domain", show_default=True,
              type=click.Choice(["domain", "phylum", "family", "genus", "species", "all"]),
              help="Taxonomic rank to plot, or 'all' for all ranks")
@click.option("--counts", "show_counts", is_flag=True,
              help="Also save a read-count version of the plot")
def domain(inputs, nodes, names, output, prefix, dpi, show, unclassified, title, palette, rank, show_counts):
    """Stacked bar chart of taxonomic composition across samples."""
    paths = _resolve_inputs(inputs)
    par, rnk, nam = load_taxonomy(nodes, names)
    df = load_kraken_files(paths, classified_only=not unclassified)
    df = add_lineage(df, par, rnk, nam)

    out = Path(output) if output else None

    ranks = ["domain", "phylum", "family", "genus", "species"] if rank == "all" else [rank]

    for r in ranks:
        plot_title = title or f"{r.capitalize()} composition across samples"
        fig = plot_domain_composition(df, rank=r, pct=True, title=plot_title, palette=palette)
        suffix = f"{r}_pct"
        _save_or_show(fig, out, prefix, suffix, show, dpi)

        if show_counts:
            fig2 = plot_domain_composition(df, rank=r, pct=False, title=plot_title, palette=palette)
            suffix2 = f"{r}_counts"
            _save_or_show(fig2, out, prefix, suffix2, show, dpi)


# ---------------------------------------------------------------------------
# breakdown sub-command
# ---------------------------------------------------------------------------

@cli.command()
@add_common_options
@click.option("--phyla", "-p", multiple=True,
              default=["Pseudomonadota", "Bacteroidota", "Bacillota",
                       "Campylobacterota"],
              show_default=True,
              help="Phylum name(s) to break down (repeatable)")
@click.option("--top-n", default=20, show_default=True,
              help="Max species/genus/family labels per phylum")
def breakdown(inputs, nodes, names, output, prefix, dpi, show, unclassified, title, palette, phyla, top_n):
    """Per-phylum species/genus/family breakdown across samples."""
    paths = _resolve_inputs(inputs)
    par, rnk, nam = load_taxonomy(nodes, names)
    df = load_kraken_files(paths, classified_only=not unclassified)
    df = add_lineage(df, par, rnk, nam)

    out = Path(output) if output else None
    fig = plot_phylum_breakdown(df, list(phyla), top_n=top_n)
    _save_or_show(fig, out, prefix, "phylum_breakdown", show, dpi)


# ---------------------------------------------------------------------------
# heatmap sub-command
# ---------------------------------------------------------------------------

@cli.command()
@add_common_options
@click.option("--rank", default="species", show_default=True,
              type=click.Choice(["domain", "phylum", "family", "genus", "species"]),
              help="Taxonomic rank for the heatmap rows")
@click.option("--top-n", default=30, show_default=True,
              help="Number of top taxa to display")
def heatmap(inputs, nodes, names, output, prefix, dpi, show, unclassified, title, palette, rank, top_n):
    """Relative-abundance heatmap at a chosen taxonomic rank."""
    paths = _resolve_inputs(inputs)
    par, rnk, nam = load_taxonomy(nodes, names)
    df = load_kraken_files(paths, classified_only=not unclassified)
    df = add_lineage(df, par, rnk, nam)

    out = Path(output) if output else None
    fig = plot_species_heatmap(df, rank=rank, top_n=top_n)
    _save_or_show(fig, out, prefix, f"heatmap_{rank}", show, dpi)


# ---------------------------------------------------------------------------
# readcounts sub-command
# ---------------------------------------------------------------------------

@cli.command()
@add_common_options
def readcounts(inputs, nodes, names, output, prefix, dpi, show, unclassified, title, palette):
    """Bar chart of total classified read counts per sample."""
    paths = _resolve_inputs(inputs)
    par, rnk, nam = load_taxonomy(nodes, names)
    df = load_kraken_files(paths, classified_only=not unclassified)
    df = add_lineage(df, par, rnk, nam)

    out = Path(output) if output else None
    fig = plot_read_counts(df)
    _save_or_show(fig, out, prefix, "readcounts", show, dpi)


# ---------------------------------------------------------------------------
# all sub-command
# ---------------------------------------------------------------------------

@cli.command(name="all")
@add_common_options
@click.option("--phyla", "-p", multiple=True,
              default=["Pseudomonadota", "Bacteroidota", "Bacillota",
                       "Campylobacterota"],
              show_default=True,
              help="Phyla for the breakdown plot (repeatable)")
@click.option("--top-n-phylum", default=10, show_default=True)
@click.option("--top-n-breakdown", default=20, show_default=True)
@click.option("--top-n-heatmap", default=30, show_default=True)
@click.option("--heatmap-rank", default="species", show_default=True,
              type=click.Choice(["domain", "phylum", "family", "genus", "species"]))
def all_plots(inputs, nodes, names, output, prefix, dpi, show, unclassified, title, palette,
              phyla, top_n_phylum, top_n_breakdown, top_n_heatmap, heatmap_rank):
    """Generate all available plots in one go."""
    paths = _resolve_inputs(inputs)
    par, rnk, nam = load_taxonomy(nodes, names)
    df = load_kraken_files(paths, classified_only=not unclassified)
    df = add_lineage(df, par, rnk, nam)

    out = Path(output) if output else None

    _save_or_show(plot_read_counts(df),              out, prefix, "readcounts",       show, dpi)
    _save_or_show(plot_domain_composition(df),        out, prefix, "domain_pct",       show, dpi)
    _save_or_show(plot_phylum_composition(df, top_n=top_n_phylum), out, prefix, "phylum_pct", show, dpi)
    _save_or_show(plot_phylum_composition(df, top_n=top_n_phylum, pct=False), out, prefix, "phylum_counts", show, dpi)
    _save_or_show(plot_phylum_breakdown(df, list(phyla), top_n=top_n_breakdown), out, prefix, "phylum_breakdown", show, dpi)
    _save_or_show(plot_species_heatmap(df, rank=heatmap_rank, top_n=top_n_heatmap), out, prefix, f"heatmap_{heatmap_rank}", show, dpi)

    click.echo("✓ All plots saved.")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    cli()


if __name__ == "__main__":
    main()
