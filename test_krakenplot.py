"""
Tests for krakenplot.

Run with:  pytest
"""

from __future__ import annotations

import io
import textwrap
from pathlib import Path

import pandas as pd
import pytest

from krakenplot.taxonomy import get_lineage, load_taxonomy
from krakenplot.loader import (
    add_lineage,
    export_taxon_count_tables,
    load_kraken_files,
    taxon_count_table,
)
from krakenplot.plots import (
    plot_domain_composition,
    plot_phylum_breakdown,
    plot_phylum_composition,
    plot_read_counts,
    plot_species_heatmap,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture()
def tmp_nodes(tmp_path):
    """Minimal nodes.dmp with root, domain (superkingdom), phylum, species."""
    lines = [
        "1\t|\t1\t|\tno rank\t|\t\t|",
        "2\t|\t1\t|\tsuperkingdom\t|\t\t|",
        "10\t|\t2\t|\tphylum\t|\t\t|",
        "100\t|\t10\t|\tfamily\t|\t\t|",
        "1000\t|\t100\t|\tgenus\t|\t\t|",
        "10000\t|\t1000\t|\tspecies\t|\t\t|",
        "3\t|\t1\t|\tsuperkingdom\t|\t\t|",
        "30\t|\t3\t|\tphylum\t|\t\t|",
        "300\t|\t30\t|\tspecies\t|\t\t|",
    ]
    p = tmp_path / "nodes.dmp"
    p.write_text("\n".join(lines))
    return p


@pytest.fixture()
def tmp_names(tmp_path):
    """Minimal names.dmp matching the nodes above."""
    lines = [
        "1\t|\troot\t|\t\t|\tscientific name\t|",
        "2\t|\tBacteria\t|\t\t|\tscientific name\t|",
        "10\t|\tFirmicutes\t|\t\t|\tscientific name\t|",
        "100\t|\tSomefamily\t|\t\t|\tscientific name\t|",
        "1000\t|\tSomegenus\t|\t\t|\tscientific name\t|",
        "10000\t|\tSomespecies\t|\t\t|\tscientific name\t|",
        "3\t|\tArchaea\t|\t\t|\tscientific name\t|",
        "30\t|\tEuryarchaeota\t|\t\t|\tscientific name\t|",
        "300\t|\tHalospecies\t|\t\t|\tscientific name\t|",
    ]
    p = tmp_path / "names.dmp"
    p.write_text("\n".join(lines))
    return p


@pytest.fixture()
def taxonomy(tmp_nodes, tmp_names):
    return load_taxonomy(tmp_nodes, tmp_names)


@pytest.fixture()
def kraken_file(tmp_path):
    """A small synthetic Kraken2 output file."""
    rows = "\n".join([
        "C\tread1\t10000\t150\tkmer_detail",
        "C\tread2\t300\t200\tkmer_detail",
        "C\tread3\t10000\t120\tkmer_detail",
        "U\tread4\t0\t180\t0",
        "C\tread5\t300\t90\tkmer_detail",
    ])
    p = tmp_path / "barcode01_output.txt"
    p.write_text(rows)
    return p


@pytest.fixture()
def kraken_file2(tmp_path):
    """A second synthetic Kraken2 output file."""
    rows = "\n".join([
        "C\tread1\t10000\t160\tkmer_detail",
        "C\tread2\t30\t130\tkmer_detail",
    ])
    p = tmp_path / "barcode02_output.txt"
    p.write_text(rows)
    return p


@pytest.fixture()
def df_with_lineage(kraken_file, kraken_file2, taxonomy):
    par, rnk, nam = taxonomy
    df = load_kraken_files([kraken_file, kraken_file2], classified_only=True)
    return add_lineage(df, par, rnk, nam)


# ---------------------------------------------------------------------------
# taxonomy tests
# ---------------------------------------------------------------------------

class TestLoadTaxonomy:
    def test_loads_parent_rank_name(self, taxonomy):
        par, rnk, nam = taxonomy
        assert par[10000] == 1000
        assert rnk[10000] == "species"
        assert nam[10000] == "Somespecies"

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            load_taxonomy(tmp_path / "nope.dmp", tmp_path / "nope2.dmp")


class TestGetLineage:
    def test_full_lineage(self, taxonomy):
        par, rnk, nam = taxonomy
        lin = get_lineage(10000, par, rnk, nam)
        assert lin["domain"] == "Bacteria"
        assert lin["phylum"] == "Firmicutes"
        assert lin["species"] == "Somespecies"

    def test_partial_lineage(self, taxonomy):
        par, rnk, nam = taxonomy
        # taxid 30 is a phylum with no genus/family in our minimal tree
        lin = get_lineage(300, par, rnk, nam)
        assert lin["phylum"] == "Euryarchaeota"
        assert lin["family"] == "unclassified"

    def test_none_taxid_returns_unassigned(self, taxonomy):
        par, rnk, nam = taxonomy
        lin = get_lineage(None, par, rnk, nam)
        assert all(v == "unclassified" for v in lin.values())

    def test_unknown_taxid_returns_unassigned(self, taxonomy):
        par, rnk, nam = taxonomy
        lin = get_lineage(999999, par, rnk, nam)
        assert all(v == "unclassified" for v in lin.values())


# ---------------------------------------------------------------------------
# loader tests
# ---------------------------------------------------------------------------

class TestLoadKrakenFiles:
    def test_loads_classified_only(self, kraken_file):
        df = load_kraken_files([kraken_file])
        assert (df["status"] == "C").all()
        assert len(df) == 4

    def test_includes_unclassified_when_requested(self, kraken_file):
        df = load_kraken_files([kraken_file], classified_only=False)
        assert "U" in df["status"].values

    def test_sample_name_extracted(self, kraken_file):
        df = load_kraken_files([kraken_file])
        assert "barcode01" in df["sample"].values

    def test_multiple_files(self, kraken_file, kraken_file2):
        df = load_kraken_files([kraken_file, kraken_file2])
        assert df["sample"].nunique() == 2

    def test_empty_list_raises(self):
        with pytest.raises(ValueError, match="No input files"):
            load_kraken_files([])


class TestAddLineage:
    def test_adds_lineage_columns(self, df_with_lineage):
        for col in ("domain", "phylum", "family", "genus", "species"):
            assert col in df_with_lineage.columns

    def test_resolved_correctly(self, df_with_lineage):
        row = df_with_lineage[df_with_lineage["taxid"] == 10000].iloc[0]
        assert row["species"] == "Somespecies"
        assert row["domain"] == "Bacteria"

    def test_missing_rank_on_classified_read_is_unresolved(self, df_with_lineage):
        row = df_with_lineage[df_with_lineage["taxid"] == 300].iloc[0]
        assert row["domain"] == "Archaea"
        assert row["family"] == "unresolved"

    def test_kraken_unclassified_read_is_unclassified(self, kraken_file, taxonomy):
        par, rnk, nam = taxonomy
        df = load_kraken_files([kraken_file], classified_only=False)
        df = add_lineage(df, par, rnk, nam)

        row = df[df["status"] == "U"].iloc[0]
        assert row["domain"] == "Unclassified"
        assert row["species"] == "Unclassified"


class TestTaxonCountTables:
    def test_counts_per_sample_and_taxon(self, df_with_lineage):
        counts = taxon_count_table(df_with_lineage, "species")

        row = counts[
            (counts["sample"] == "barcode01")
            & (counts["species"] == "Somespecies")
        ].iloc[0]
        assert row["count"] == 2

    def test_unclassified_rows_are_exported_when_loaded(self, kraken_file, taxonomy):
        par, rnk, nam = taxonomy
        df = load_kraken_files([kraken_file], classified_only=False)
        df = add_lineage(df, par, rnk, nam)

        counts = taxon_count_table(df, "domain")

        row = counts[
            (counts["sample"] == "barcode01")
            & (counts["domain"] == "Unclassified")
        ].iloc[0]
        assert row["count"] == 1

    def test_export_writes_prefixed_tsvs(self, df_with_lineage, tmp_path):
        paths = export_taxon_count_tables(
            df_with_lineage, tmp_path, "sample_prefix", ["domain", "species"]
        )

        assert paths == [
            tmp_path / "sample_prefix_domain_counts.tsv",
            tmp_path / "sample_prefix_species_counts.tsv",
        ]
        exported = pd.read_csv(paths[1], sep="\t")
        assert list(exported.columns) == ["sample", "species", "count"]


# ---------------------------------------------------------------------------
# plot tests (smoke tests — check figures are returned without error)
# ---------------------------------------------------------------------------

class TestPlots:
    def test_plot_read_counts(self, df_with_lineage):
        fig = plot_read_counts(df_with_lineage)
        assert fig is not None

    def test_plot_domain_composition_pct(self, df_with_lineage):
        fig = plot_domain_composition(df_with_lineage, pct=True)
        assert fig is not None

    def test_plot_domain_composition_counts(self, df_with_lineage):
        fig = plot_domain_composition(df_with_lineage, pct=False)
        assert fig is not None

    def test_plot_phylum_composition(self, df_with_lineage):
        fig = plot_phylum_composition(df_with_lineage, top_n=5)
        assert fig is not None

    def test_plot_phylum_breakdown_present(self, df_with_lineage):
        fig = plot_phylum_breakdown(df_with_lineage, ["Firmicutes"])
        assert fig is not None

    def test_plot_phylum_breakdown_missing(self, df_with_lineage):
        # Should not raise; should return a "no data" figure
        fig = plot_phylum_breakdown(df_with_lineage, ["NonexistentPhylum"])
        assert fig is not None

    def test_plot_species_heatmap(self, df_with_lineage):
        fig = plot_species_heatmap(df_with_lineage, rank="species", top_n=5)
        assert fig is not None

    def test_plot_species_heatmap_bad_rank(self, df_with_lineage):
        with pytest.raises(ValueError, match="Rank"):
            plot_species_heatmap(df_with_lineage, rank="bogus")
