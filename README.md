<p align="center">
  <img src="ruta/a/tu/logo.png" alt="PopBugGenie logo" width="180"/>
  <img src="assets/PopBugGenie.png" alt="PopBugGenie logo" width="180"/>
</p>

<h1 align="center">🪄 PopBugGenie</h1>
<p align="center"><em>Unleashing the magic of microbial diversity</em></p>

---

## 🚀 Features

- ✅ Preprocessing of VCF files for compatibility with `scikit-allel`
- 🔍 Automatic population inference using PCA + KMeans
- 📊 Calculation of genetic diversity statistics globally and by population:
  - Nucleotide diversity (π)
  - Haplotype diversity (Hd)
  - Tajima’s D
  - Linkage disequilibrium (LD, measured as r²)
- 🌍 Window-based FST calculation between all population pairs
- 🔥 Z-score normalization and outlier detection (±2.326, i.e., 1st and 99th percentiles)
- 📈 High-quality visualizations:
  - PCA plots
  - Z-score tracks per population
  - Heatmaps combining multiple statistics
  - FST boxplots and heatmaps

---

## 🧬 Installation

PopBugGenie is written in **Python 3** and depends on:

- `scikit-allel`
- `numpy`
- `pandas`
- `matplotlib`
- `seaborn`
- `scikit-learn`

Install the dependencies:

```bash
conda create -n popbuggenie python=3.10
conda activate popbuggenie
pip install numpy pandas matplotlib seaborn scikit-allel scikit-learn
```

---

## 📂 Input

- A **multi-sample VCF file** (`.vcf` or `.vcf.gz`) with variant calls across bacterial samples.
- (Optional) Use `--replace-missing` to fill missing genotypes (`./.`), assuming clonality.

---

## 🧪 Usage

```bash
python PopBugGenie.py \
  --vcf path/to/your.vcf.gz \
  --output output_directory \
  --replace-missing \
  --window 50000 \
  --max_clusters 6 \
  --min_samples_per_pop 3
```

### Parameters

| Flag                   | Description                                                                 |
|------------------------|-----------------------------------------------------------------------------|
| `--vcf`                | Path to the input multi-sample VCF file                                     |
| `--output`             | Output directory for results                                                |
| `--replace-missing`    | (Optional) Replace missing genotypes (`./.`) with `0/0`                     |
| `--window`             | Non-overlapping genomic window size in base pairs (default: 50000)                          |
| `--max_clusters`       | Max number of KMeans clusters for PCA inference (default: 8)                |
| `--min_samples_per_pop`| Minimum samples per population to compute statistics (default: 3)           |

---

## 📊 Output Structure

```
output_directory/
├── auto_populations.csv                 # Cluster assignments from PCA + KMeans
├── global/                              # Global statistics (π, Tajima’s D, LD, Hd)
├── by_population/                       # Stats and plots per population
├── fst_pairwise/                        # Pairwise FST results and plots
├── comparison_plots/                    # Combined heatmaps and Z-score figures
├── zscore_summary_table.csv             # Summary table with combined statistics
└── pipeline_log_TIMESTAMP.log           # Log file for the pipeline run
```

---

## 📖 Citation

If you use PopBugGenie in your work, please cite:

> **Katalina Urbano Chaves (2025)**  
> *PopBugGenie: A pipeline for detecting genomic adaptation signatures in bacterial populations.*  
> Master’s Thesis, Universitat Autònoma de Barcelona.

---

## ❤️ Acknowledgements

PopBugGenie was developed as part of a Master's thesis in Bioinformatics at UAB.  
Special thanks to the open-source bioinformatics community and the authors of tools like:

- `scikit-allel`
- `bcftools`
- `vcftools`
- `snippy`
