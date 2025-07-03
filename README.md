# PopBugGenie
PopBugGenie
Unleashing the magic of microbial diversity

PopBugGenie is an automated and modular pipeline designed to analyze genetic diversity, population structure, selection, and differentiation in bacterial genomes from multi-sample VCF files. It provides a comprehensive and visual overview of genomic variation through Z-score normalized statistics and integrates population inference, sliding-window analyses, and graphical outputs.

🚀 Features
✅ Preprocessing of VCF files for compatibility with scikit-allel

🔍 Automatic population inference using PCA + KMeans

📊 Calculation of genetic diversity statistics globally and by population:

Nucleotide diversity (π)

Haplotype diversity (Hd)

Tajima’s D

Linkage disequilibrium (LD, measured as r²)

🌍 Window-based FST calculation between all population pairs

🔥 Z-score normalization and outlier detection (±2.326, i.e., 1st and 99th percentiles)

📈 High-quality visualizations:

PCA plots

Z-score tracks per population

Heatmaps combining multiple statistics

FST boxplots and heatmaps

🧬 Installation
PopBugGenie is written in Python 3 and depends on:

scikit-allel

numpy

pandas

matplotlib

seaborn

scikit-learn

You can install the dependencies via:

bash
Copiar
Editar
conda create -n popbuggenie python=3.10
conda activate popbuggenie
pip install numpy pandas matplotlib seaborn scikit-allel scikit-learn
📂 Input
A multi-sample VCF file (.vcf or .vcf.gz) with variant calls across bacterial samples

(Optional) Use the --replace-missing flag to replace missing genotypes with 0/0, assuming clonality

🧪 Usage
bash
Copiar
Editar
python PopBugGenie.py \
    --vcf path/to/your.vcf.gz \
    --output output_directory \
    --replace-missing \
    --window 50000 \
    --max_clusters 6 \
    --min_samples_per_pop 3
Arguments:
Flag	Description
--vcf	Path to the input multi-sample VCF file
--output	Output directory where results will be saved
--replace-missing	(Optional) Replace missing genotypes (./.) with 0/0
--window	Sliding window size in base pairs (default: 50000)
--max_clusters	Max number of clusters to test during population inference (default: 8)
--min_samples_per_pop	Minimum number of samples required per population (default: 3)

📊 Outputs
The pipeline produces:

bash
Copiar
Editar
output_directory/
├── auto_populations.csv                 # Clustered sample assignments
├── global/                              # Global statistics (π, Tajima’s D, LD, Hd)
├── by_population/
│   ├── population_0/
│   │   ├── pi.csv, tajima_d.csv, ld.csv, hd_by_window.csv
│   │   └── plots (.png)
│   ├── population_1/
│   └── ...
├── fst_pairwise/                        # FST matrices, boxplots, and comparisons
├── comparison_plots/                    # Z-score plots and combined heatmaps
├── zscore_summary_table.csv             # Final table with all Z-scores
└── pipeline_log_TIMESTAMP.log           # Log of the run
📖 Citation
If you use PopBugGenie in your work, please cite:

Katalina Urbano Chaves (2025). PopBugGenie: A pipeline for detecting genomic adaptation signatures in bacterial populations. Master’s Thesis. Universitat Autònoma de Barcelona.

❤️ Acknowledgements
PopBugGenie was developed as part of a Master’s Thesis in Bioinformatics at UAB. Special thanks to the bioinformatics community and open-source developers who maintain tools like scikit-allel, bcftools, and snippy.
