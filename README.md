# Kraken2Plot

**A command-line tool for creating visualizations of Kraken2 taxonomic classification output.**

![Version](https://img.shields.io/badge/version-0.1.0-blue) ![License](https://img.shields.io/badge/license-MIT-green)

## Overview

Kraken2Plot takes the output from [Kraken2](https://ccb.jhu.edu/software/kraken2/) (a metagenomic sequence classifier) and transforms it into stacked bar charts and heatmaps that show the taxonomic composition of your samples.

Whether you're analyzing wastewater, environmental samples, or other microbial communities, Kraken2Plot makes it easy to visualize and compare the taxonomic profiles across multiple samples!

## Features

✨ **Key Capabilities:**

- 📊 **Stacked bar charts** showing taxonomic composition across samples
- 📈 **Multiple taxonomic ranks**: domain, phylum, family, genus, or species (or plot all at once!)
- 🎨 **Customizable color palettes**: Including color-blind friendly options (viridis, plasma)
- 🎯 **Flexible sample naming**: Automatically extracts barcode IDs (often used in Nanopore sequencing workflows) or uses full sample names
- 📂 **Organized output**: Save multiple plots to a directory with custom file prefix
- 🔧 **Consistent styling**: Professional publication-quality plots with consistent formatting
- 💾 **Count and percentage options**: View raw read counts or percentages
- 🧬 **Handles real-world data**: Automatically deals with unclassified reads and missing taxa

## Requirements

- **Python 3.9+**
- **NCBI Taxonomy files**: `nodes.dmp` and `names.dmp` from [NCBI's taxonomy dump](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/)

### Dependencies

The following Python packages are installed automatically:
- pandas
- matplotlib
- numpy
- click

## Installation

### 1. Clone or download this repository

```bash
git clone https://github.com/yourusername/kraken2plot.git
cd kraken2plot
```

### 2. Install the package in development mode

```bash
pip install -e .
```

This installs Kraken2Plot and makes the `kraken2plot` command available system-wide.

### 3. Download NCBI taxonomy files

Kraken2Plot needs the NCBI taxonomy database files to map taxon IDs to names:

```bash
# Create a directory to store taxonomy files
mkdir -p ~/databases/ncbi_taxonomy
cd ~/databases/ncbi_taxonomy

# Download and extract the taxonomy dump
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzf taxdump.tar.gz

# You now have nodes.dmp and names.dmp in this directory
```

## Quick Start

### Basic usage

```bash
# Set up paths to your files
data="/path/to/kraken2/output/directory"
nodes="~/databases/ncbi_taxonomy/nodes.dmp"
names="~/databases/ncbi_taxonomy/names.dmp"

# Generate a phylum-level composition plot
kraken2plot phylum \
  --nodes "$nodes" \
  --names "$names" \
  "$data"
```

This creates a plot file `kraken2plot_phylum_pct.png` in your current directory.

### Real-world example

For a complete analysis of wastewater samples with all taxonomic ranks:

```bash
data="/path/to/kraken/output"
nodes="~/databases/ncbi_taxonomy/nodes.dmp"
names="~/databases/ncbi_taxonomy/names.dmp"

kraken2plot phylum \
  --nodes "$nodes" \
  --names "$names" \
  --output ./results \
  --prefix composition \
  --title "Composition of wastewater DNA March 2026" \
  --palette tab10 \
  --rank all \
  "$data"
```

This generates all taxonomic rank plots in the `results/` directory with file names like:
- `composition_domain_pct.png`
- `composition_phylum_pct.png`
- `composition_family_pct.png`
- `composition_genus_pct.png`
- `composition_species_pct.png`

## Usage Guide

### Main Command: `kraken2plot phylum`

The primary command for generating taxonomic composition plots.

#### Required Arguments

- **INPUT**: Path(s) to Kraken2 output file(s) or a directory containing them
  - Can be a single file: `kraken_output.txt`
  - Can be a directory: `kraken_results/` (automatically finds `*_kraken_output.txt` files)
  - Can be a glob pattern: `kraken_results/*.txt`

#### Required Options

- `--nodes PATH`: Path to NCBI `nodes.dmp` file
- `--names PATH`: Path to NCBI `names.dmp` file

#### Output Options

- `--output OUTPUT_DIR` (or `-o`): Directory where to save plot files (default: current directory)
- `--prefix PREFIX`: File name prefix for output (default: `kraken2plot`)
  - Files will be saved as `{prefix}_{rank}_{type}.png`

#### Plot Customization

- `--rank {domain|phylum|family|genus|species|all}`: Which taxonomic rank to plot (default: `phylum`)
  - Use `all` to generate plots for every rank in one command
- `--title TITLE`: Custom plot title (default: auto-generated like "Phylum composition across samples")
- `--palette {tab20|viridis|plasma|Set1|Set2|tab10|Paired}`: Color palette (default: `tab20`)
  - `viridis` and `plasma` are color-blind friendly
- `--top-n N`: Number of top taxa to show individually; others grouped as "Other" (default: 10)
  - Not used for domain-level plots
- `--counts`: Also generate count-based versions of plots (default: off, shows percentages only)

#### Data Options

- `--unclassified`: Include unclassified reads in the plot (default: off, excludes them)
  - Unclassified reads appear in gray
- `--dpi DPI`: Resolution of output images in dots per inch (default: 150)
  - Use higher values (300+) for publication-quality figures

#### Interactive Options

- `--show`: Display plot in interactive viewer instead of saving to file

### Other Commands

#### `kraken2plot domain`
Similar to `phylum` but defaults to domain-level composition.

#### `kraken2plot breakdown`
Per-phylum species/genus/family breakdown. Useful for detailed examination of specific phyla.

```bash
kraken2plot breakdown \
  --nodes "$nodes" \
  --names "$names" \
  --phyla Pseudomonadota Bacteroidota \
  "$data"
```

#### `kraken2plot heatmap`
Relative abundance heatmap for a specific taxonomic rank.

```bash
kraken2plot heatmap \
  --nodes "$nodes" \
  --names "$names" \
  --rank species \
  --top-n 30 \
  "$data"
```

## Input File Format

Kraken2Plot expects standard Kraken2 output files (per-read output format):

```
C	read1	taxid	read_length	kmer_hits
C	read2	taxid	read_length	kmer_hits
U	read3	0	read_length	kmer_hits
```

Where:
- Column 1: Classification status (`C` = classified, `U` = unclassified)
- Column 2: Read/sequence identifier
- Column 3: NCBI taxonomy ID (0 for unclassified)
- Column 4: Sequence length
- Column 5: K-mer hits

## Sample Name Extraction

Kraken2Plot intelligently extracts sample names from file names:

| Filename | Sample Name |
|----------|------------|
| `barcode01_kraken_output.txt` | `barcode01` |
| `barcode22_kraken_output.txt` | `barcode22` |
| `unclassified_kraken_output.txt` | `unclassified` |
| `Antwerpen_Zuid_influent_kraken_output.txt` | `Antwerpen_Zuid_influent` |

## Output

Kraken2Plot generates PNG files with the following naming convention:

`{prefix}_{rank}_{type}.png`

Where:
- `prefix`: Custom prefix you specified (default: `kraken2plot`)
- `rank`: The taxonomic rank (domain, phylum, family, genus, species)
- `type`: Either `pct` (percentage) or `counts` (read counts)

Example outputs:
- `composition_phylum_pct.png` — Percentage-based phylum plot
- `composition_phylum_counts.png` — Count-based phylum plot (if `--counts` used)

## Color Handling

- **Unclassified reads**: Always shown in gray when included with `--unclassified`
- **"Other" category**: When taxa exceed `--top-n`, they're grouped into "Other" for clarity
- **Custom palettes**: Choose from matplotlib's colormaps; `viridis` and `plasma` are recommended for accessibility

## Examples

### Example 1: Simple phylum composition
```bash
kraken2plot phylum \
  --nodes nodes.dmp \
  --names names.dmp \
  kraken_results/
```

### Example 2: Compare genus composition across samples with counts
```bash
kraken2plot phylum \
  --nodes nodes.dmp \
  --names names.dmp \
  --rank genus \
  --counts \
  sample1.txt sample2.txt sample3.txt
```

### Example 3: Complete analysis with all ranks, custom styling
```bash
kraken2plot phylum \
  --nodes ~/db/nodes.dmp \
  --names ~/db/names.dmp \
  --output ./plots \
  --prefix wastewater_march2026 \
  --title "Wastewater Microbial Communities" \
  --palette viridis \
  --rank all \
  --top-n 15 \
  --unclassified \
  kraken_data/
```

### Example 4: Publication-quality figure (high DPI)
```bash
kraken2plot phylum \
  --nodes nodes.dmp \
  --names names.dmp \
  --dpi 300 \
  --palette plasma \
  --title "My Figure Title" \
  data.txt
```

## Troubleshooting

### Error: "No input files found"
- Check that your file paths are correct
- Ensure files end in `.txt`
- If using a directory, make sure it contains `*_kraken_output.txt` files

### Error: "ModuleNotFoundError: No module named 'kraken2plot'"
- Make sure you've installed the package: `pip install -e .`
- Check that you're in the correct directory

### Error: "No such file or directory: nodes.dmp"
- Download NCBI taxonomy files (see Installation section)
- Provide full path to the files, e.g., `~/databases/ncbi_taxonomy/nodes.dmp`

### Plots aren't saved to the output directory
- Use `--output ./your_directory` (with `./` prefix for relative paths)
- The directory will be created automatically if it doesn't exist

### Sample names look wrong
- For barcode files: Check that filenames follow the pattern `*_barcodeXX_kraken_output.txt`
- For other files: The full filename (minus `_kraken_output.txt`) is used as the sample name

## Tips for Lab Technicians

📝 **Simple workflow:**

1. **Run Kraken2** on your samples (produces `*.txt` output files)
2. **Download taxonomy files** once (reuse for all future analyses)
3. **Create a results directory**:
   ```bash
   mkdir results
   ```
4. **Run Kraken2Plot** with a descriptive prefix and title:
   ```bash
   kraken2plot phylum \
     --nodes ~/databases/ncbi_taxonomy/nodes.dmp \
     --names ~/databases/ncbi_taxonomy/names.dmp \
     --output results \
     --prefix my_project_name \
     --title "My Project - Sample Set" \
     ~/kraken_output/
   ```
5. **Open the generated PNG files** in your favorite image viewer or add to reports

🎨 **Color selection guide:**
- `tab20`: Default, distinct colors for up to 20 taxa
- `viridis`: Great for presentations, color-blind friendly
- `plasma`: Vibrant, also color-blind friendly
- `Set1`, `Set2`, `tab10`: More subtle, fewer colors

## Citation

```
Wim L. Cuypers. Kraken2Plot: Visualization tool for Kraken2 taxonomic output. 
Version 0.1.0. GitHub. 2026. Available at: https://github.com/Cuypers-Wim/kraken2plot

```

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## Authors

- Dr. Wim L. Cuypers

## Support

For questions or issues:
1. Check the [Troubleshooting](#troubleshooting) section above
2. Open an issue on GitHub
3. Contact the development team

---

**Happy plotting! 🧬📊**
