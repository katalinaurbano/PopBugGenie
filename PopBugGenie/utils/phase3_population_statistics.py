import allel
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import pickle

def calcular_estadisticos_por_poblacion(vcf_path, poblaciones_csv, output_dir, window_size=50000, min_samples=3):
    print("\n[FASE 3] Calculando estadísticos por población...")
    os.makedirs(output_dir, exist_ok=True)

    df_pobs = pd.read_csv(poblaciones_csv)
    muestras = df_pobs['Sample'].tolist()
    clusters = df_pobs['Cluster'].tolist()
    poblaciones = sorted(df_pobs['Cluster'].unique())
    print(f"Poblaciones encontradas: {poblaciones}")

    callset = allel.read_vcf(vcf_path, samples=muestras)
    genotypes = allel.GenotypeArray(callset['calldata/GT'])
    positions = callset['variants/POS']

    sorted_idx = np.argsort(positions)
    positions = positions[sorted_idx]
    genotypes = genotypes.take(sorted_idx, axis=0)

    sample_names = callset['samples']
    df_sample_cluster = pd.DataFrame({"sample": sample_names, "cluster": clusters})

    subpops = [df_sample_cluster[df_sample_cluster['cluster'] == pop].index.tolist() for pop in poblaciones]
    with open(os.path.join(output_dir, 'subpops_indices.pkl'), 'wb') as f:
        pickle.dump(subpops, f)

    window_starts = np.arange(positions.min(), positions.max(), window_size)
    windows = np.column_stack((window_starts, window_starts + window_size))

    for i, pop in enumerate(poblaciones):
        if len(subpops[i]) < min_samples:
            print(f"⚠️  Saltando población {pop} (solo {len(subpops[i])} muestra(s))")
            continue

        pop_dir = os.path.join(output_dir, f"poblacion_{pop}")
        os.makedirs(pop_dir, exist_ok=True)

        gen_pop = genotypes.take(subpops[i], axis=1)
        ac_pop = gen_pop.count_alleles()

        # π
        pi_total, _, _, n_sites = allel.windowed_diversity(pos=positions, ac=ac_pop, windows=windows)
        pi = pi_total / n_sites
        pi_z = (pi - np.nanmean(pi)) / np.nanstd(pi)
        pd.DataFrame({"start": windows[:, 0], "end": windows[:, 1], "pi": pi, "pi_zscore": pi_z}).to_csv(os.path.join(pop_dir, 'pi.csv'), index=False)
        plt.figure()
        plt.plot(windows[:, 0], pi_z)
        plt.axhline(y=2.326, color='red', linestyle='--', label='Z = +2.326')
        plt.xlabel('Posición (bp)')
        plt.ylabel('π (z-score)')
        plt.title(f'Diversidad nucleotídica - Población {pop}')
        plt.tight_layout()
        plt.savefig(os.path.join(pop_dir, 'pi.png'))
        plt.close()

        # Hd por ventana
        gn_pop = gen_pop.to_n_alt()
        polymorphic = np.any((gn_pop > 0) & (gn_pop < len(subpops[i])), axis=1)
        no_missing = ~np.any(gn_pop < 0, axis=1)
        mask_utiles = polymorphic & no_missing

        gen_filtrado = gen_pop.take(np.where(mask_utiles)[0], axis=0)
        pos_filtradas = positions[mask_utiles]

        # Hd global
        haps = gen_filtrado.to_haplotypes()
        with open(os.path.join(pop_dir, 'hd.txt'), 'w') as f:
            if haps.shape[0] == 0:
                f.write("No se pudo calcular Hd: sin SNPs polimórficos sin missing.\n")
            else:
                hd = allel.haplotype_diversity(haps.T)
                hd_mean = np.nanmean(hd)
                f.write(f"Hd promedio: {hd_mean:.4f}\n")

        # Hd por ventana
        hd_vals = []
        for start, end in windows:
            mask_win = (pos_filtradas >= start) & (pos_filtradas < end)
            snps_window = gen_filtrado.compress(mask_win, axis=0)
            if snps_window.shape[0] < 2:
                hd_vals.append(np.nan)
            else:
                haps_win = snps_window.to_haplotypes()
                hd = allel.haplotype_diversity(haps_win.T)
                hd_vals.append(np.nanmean(hd))
        hd_vals = np.array(hd_vals)
        hd_z = (hd_vals - np.nanmean(hd_vals)) / np.nanstd(hd_vals)
        pd.DataFrame({"start": windows[:, 0], "end": windows[:, 1], "hd": hd_vals, "hd_zscore": hd_z}).to_csv(os.path.join(pop_dir, 'hd_por_ventana.csv'), index=False)
        plt.figure()
        plt.plot(windows[:, 0], hd_z)
        plt.axhline(y=2.326, color='red', linestyle='--', label='Z = +2.326')
        plt.xlabel('Inicio ventana (bp)')
        plt.ylabel('Hd (z-score)')
        plt.title(f'Diversidad haplotípica - Población {pop}')
        plt.tight_layout()
        plt.savefig(os.path.join(pop_dir, 'hd.png'))
        plt.close()

        # Tajima's D
        if len(pos_filtradas) < 2:
            tajima_d = np.array([np.nan] * len(windows))
        else:
            ac_filtrado = gen_filtrado.count_alleles()
            tajima_d, _, _ = allel.windowed_tajima_d(pos=pos_filtradas, ac=ac_filtrado, windows=windows)
        tajima_z = (tajima_d - np.nanmean(tajima_d)) / np.nanstd(tajima_d)
        pd.DataFrame({"start": windows[:, 0], "end": windows[:, 1], "tajima_d": tajima_d, "tajima_d_zscore": tajima_z}).to_csv(os.path.join(pop_dir, 'tajima_d.csv'), index=False)
        plt.figure()
        plt.plot(windows[:, 0], tajima_z)
        plt.axhline(y=2.326, color='red', linestyle='--', label='Z = +2.326')
        plt.xlabel('Posición')
        plt.ylabel("Tajima's D (z-score)")
        plt.title(f"Tajima's D - Población {pop}")
        plt.tight_layout()
        plt.savefig(os.path.join(pop_dir, 'tajima.png'))
        plt.close()

        # LD
        gn = gn_pop[mask_utiles]
        pos_mask = positions[mask_utiles]
        ld_vals = []
        for start, end in windows:
            win_mask = (pos_mask >= start) & (pos_mask < end)
            gn_window = gn[win_mask, :]
            if gn_window.shape[0] < 2:
                ld_vals.append(np.nan)
            else:
                r = allel.rogers_huff_r(gn_window)
                r2 = r ** 2
                ld_vals.append(np.nanmean(r2))
        ld_z = (np.array(ld_vals) - np.nanmean(ld_vals)) / np.nanstd(ld_vals)
        pd.DataFrame({"start": windows[:, 0], "end": windows[:, 1], "ld_r2": ld_vals, "ld_r2_zscore": ld_z}).to_csv(os.path.join(pop_dir, 'ld.csv'), index=False)
        plt.figure()
        plt.plot(windows[:, 0], ld_z)
        plt.axhline(y=2.326, color='red', linestyle='--', label='Z = +2.326')
        plt.xlabel('Inicio ventana (bp)')
        plt.ylabel('LD (r² z-score)')
        plt.title(f'Desequilibrio ligamiento - Población {pop}')
        plt.tight_layout()
        plt.savefig(os.path.join(pop_dir, 'ld.png'))
        plt.close()

    # PCA coloreado
    gn_alt = genotypes.to_n_alt()
    gn_alt = gn_alt[~np.any(np.isnan(gn_alt), axis=1)]
    gn_scaled = StandardScaler().fit_transform(gn_alt.T)
    pca = PCA(n_components=2)
    coords = pca.fit_transform(gn_scaled)
    explained = pca.explained_variance_ratio_ * 100
    df_pca = pd.DataFrame(coords, columns=['PC1', 'PC2'])
    df_pca['cluster'] = clusters
    plt.figure()
    for pop in poblaciones:
        sub = df_pca[df_pca['cluster'] == pop]
        plt.scatter(sub['PC1'], sub['PC2'], label=f"Población {pop}", alpha=0.7)
    plt.legend()
    plt.xlabel(f'PC1 ({explained[0]:.1f}%)')
    plt.ylabel(f'PC2 ({explained[1]:.1f}%)')
    plt.title('PCA por población')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'pca_coloreado.png'))
    plt.close()

    print("✅ Estadísticos por población calculados y guardados en:", output_dir)

