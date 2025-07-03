<p align="center">
  <img src="ruta/a/tu/logo.png" alt="PopBugGenie logo" width="180"/>
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

Install dependencies:

```bash
conda create -n popbuggenie python=3.10
conda activate popbuggenie
pip install numpy pandas matplotlib seaborn scikit-allel scikit-learn
