import allel
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

def calcular_estadisticos_globales(vcf_path, output_dir, window_size=50000):
    print("\n[FASE 2] Calculando estadísticos globales...")
    os.makedirs(output_dir, exist_ok=True)

    callset = allel.read_vcf(vcf_path)
    genotypes = allel.GenotypeArray(callset['calldata/GT'])
    positions = callset['variants/POS']

    sorted_indices = np.argsort(positions)
    positions = positions[sorted_indices]
    genotypes = genotypes.take(sorted_indices, axis=0)

    window_starts = np.arange(positions.min(), positions.max(), window_size)
    windows = np.column_stack((window_starts, window_starts + window_size))

    ac = genotypes.count_alleles()

    # π
    pi_total, _, _, n_sites = allel.windowed_diversity(pos=positions, ac=ac, windows=windows)
    pi = pi_total / n_sites
    pi_z = (pi - np.nanmean(pi)) / np.nanstd(pi)
    pd.DataFrame({"start": windows[:,0], "end": windows[:,1], "pi": pi, "pi_zscore": pi_z}).to_csv(os.path.join(output_dir, 'pi_global.csv'), index=False)
    plt.figure()
    plt.plot(windows[:, 0], pi_z)
    plt.axhline(y=2.326, color='red', linestyle='--', label='Z = +2.326')
    plt.xlabel('Posición (bp)')
    plt.ylabel('π (z-score)')
    plt.title('π global (normalizado)')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'pi_global.png'))

    # Tajima's D
    tajima_d, _, _ = allel.windowed_tajima_d(pos=positions, ac=ac, windows=windows)
    tajima_z = (tajima_d - np.nanmean(tajima_d)) / np.nanstd(tajima_d)
    pd.DataFrame({"start": windows[:,0], "end": windows[:,1], "tajima_d": tajima_d, "tajima_d_zscore": tajima_z}).to_csv(os.path.join(output_dir, 'tajima_d_global.csv'), index=False)
    plt.figure()
    plt.plot(windows[:, 0], tajima_z)
    plt.axhline(y=2.326, color='red', linestyle='--', label='Z = +2.326')
    plt.xlabel('Posición (bp)')
    plt.ylabel("Tajima's D (z-score)")
    plt.title("Tajima's D global (normalizado)")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'tajima_global.png'))

    gn_alt = genotypes.to_n_alt()
    mask = ~np.any(np.isnan(gn_alt), axis=1)
    gn_alt = gn_alt[mask]
    pos_filtered = positions[mask]

    #LD
    ld_means = []
    for start, end in windows:
        mask_window = (pos_filtered >= start) & (pos_filtered < end)
        gn_window = gn_alt[mask_window, :]
        if gn_window.shape[0] < 2:
            ld_means.append(np.nan)
        else:
            r = allel.rogers_huff_r(gn_window)
            r2 = r ** 2
            ld_means.append(np.nanmean(r2))

    ld_z = (np.array(ld_means) - np.nanmean(ld_means)) / np.nanstd(ld_means)
    pd.DataFrame({"start": windows[:,0], "end": windows[:,1], "ld_r2": ld_means, "ld_r2_zscore": ld_z}).to_csv(os.path.join(output_dir, 'ld_global.csv'), index=False)
    plt.figure()
    plt.plot(windows[:, 0], ld_z)
    plt.axhline(y=2.326, color='red', linestyle='--', label='Z = +2.326')
    plt.xlabel('Posición (bp)')
    plt.ylabel('LD (r² z-score)')
    plt.title('LD global (normalizado)')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'ld_global.png'))


    #Hd
    haps = genotypes.to_haplotypes()
    haps = haps[~np.any(haps < 0, axis=1)]
    print("Número de SNPs sin missing para Hd:", haps.shape[0])
    with open(os.path.join(output_dir, 'hd_global.txt'), 'w') as f:
        if haps.shape[0] == 0:
            print("⚠️  No se pudo calcular la diversidad haplotípica global (Hd).")
            print("ℹ️  Esto ocurre cuando no hay suficientes SNPs polimórficos sin datos faltantes en TODAS las muestras.")
            print("ℹ️  Revisa si la variación genética está distribuida en subgrupos: el cálculo por población puede seguir siendo válido.")
            f.write("⚠️  No se pudo calcular la diversidad haplotípica global (Hd).\n")
            f.write("ℹ️  Esto puede deberse a que no hay suficientes SNPs polimórficos sin datos faltantes en todas las muestras.\n")
            f.write("ℹ️  Se recomienda revisar los resultados por población, donde puede existir diversidad haplotípica local.\n")
        else:
            num_haplotipos = np.unique(haps.T, axis=0).shape[0]
            print("Número de haplotipos únicos:", num_haplotipos)
            hd = allel.haplotype_diversity(haps.T)
            hd_mean = np.nanmean(hd)
            f.write(f"Diversidad haplotípica promedio: {hd_mean:.4f}\n")
            f.write(f"Número de haplotipos únicos: {num_haplotipos}\n")

    # Hd por ventanas (diversidad haplotípica)
    print("\n▶️ Calculando Hd por ventanas...")

    haps = genotypes.to_haplotypes()
    hd_vals = []

    for start, end in windows:
        # Filtrar posiciones en la ventana
        mask = (positions >= start) & (positions < end)
        haps_window = haps[mask]

        # Quitar SNPs con missing
        haps_window = haps_window[~np.any(haps_window < 0, axis=1)]

        # Necesitamos al menos 3 SNPs sin missing para que tenga sentido
        if haps_window.shape[0] < 3:
            hd_vals.append(np.nan)
            continue

        try:
            hd = allel.haplotype_diversity(haps_window.T)
            hd_mean = np.nanmean(hd)

            # Proteger contra valores inválidos o demasiado bajos
            if hd_mean is None or np.isnan(hd_mean) or hd_mean < 1e-6:
                hd_vals.append(np.nan)
            else:
                hd_vals.append(hd_mean)

        except Exception as e:
            print(f"❌ Error calculando Hd en ventana {start}-{end}: {e}")
            hd_vals.append(np.nan)

    hd_array = np.array(hd_vals)
    hd_z = (hd_array - np.nanmean(hd_array)) / np.nanstd(hd_array)

    # Guardar CSV
    df_hd = pd.DataFrame({
        "start": windows[:, 0],
        "end": windows[:, 1],
        "hd": hd_array,
        "hd_zscore": hd_z
    })
    df_hd.to_csv(os.path.join(output_dir, 'hd_global.csv'), index=False)

    # Graficar
    plt.figure()
    plt.plot(windows[:, 0], hd_z)
    plt.axhline(y=2.326, color='red', linestyle='--', label='Z = +2.326')
    plt.axhline(y=-2.326, color='red', linestyle='--')
    plt.xlabel('Posición (bp)')
    plt.ylabel('Hd (z-score)')
    plt.title('Hd global (normalizado)')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'hd_global.png'))

    print("✅ Hd por ventana calculado y guardado en:", os.path.join(output_dir, "hd_global.csv"))


    #PCA global

    gn_t = gn_alt.T
    gn_scaled = StandardScaler().fit_transform(gn_t)
    pca = PCA(n_components=2)
    coords = pca.fit_transform(gn_scaled)
    explained = pca.explained_variance_ratio_ * 100
    df_pca = pd.DataFrame(coords, columns=['PC1', 'PC2'])
    df_pca.to_csv(os.path.join(output_dir, 'pca_global.csv'), index=False)
    plt.figure()
    plt.scatter(df_pca['PC1'], df_pca['PC2'], alpha=0.7)
    plt.xlabel(f'PC1 ({explained[0]:.1f}%)')
    plt.ylabel(f'PC2 ({explained[1]:.1f}%)')
    plt.title('PCA global')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'pca_global.png'))

    print("✅ Estadísticos globales calculados y guardados en:", output_dir)

