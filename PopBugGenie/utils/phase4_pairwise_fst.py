import os
import numpy as np
import pandas as pd
import allel
import matplotlib.pyplot as plt
from itertools import combinations


def calcular_fst_por_pares(vcf_path, poblaciones_csv, output_dir, window_size=50000, min_samples=3):
    os.makedirs(output_dir, exist_ok=True)

    print("\n[FASE 4] Calculando FST por pares de poblaciones...")

    # Leer datos
    callset = allel.read_vcf(vcf_path)
    genotypes = allel.GenotypeArray(callset['calldata/GT'])
    positions = callset['variants/POS']
    muestras = callset['samples']

    df = pd.read_csv(poblaciones_csv)
    df.columns = df.columns.str.lower()

    # Agrupar muestras por cluster
    cluster_dict = {
        str(cluster): df[df['cluster'] == cluster]['sample'].tolist()
        for cluster in df['cluster'].unique()
    }

    # Filtrar clusters con suficientes muestras
    subpops = []
    nombres = []
    for cluster, muestras_cluster in cluster_dict.items():
        indices = [np.where(muestras == s)[0][0] for s in muestras_cluster if s in muestras]
        if len(indices) >= min_samples:
            subpops.append(indices)
            nombres.append(f"Pop{cluster}")
        else:
            print(f"⚠️  Saltando cluster {cluster}: solo {len(indices)} muestra(s)")

    if len(subpops) < 2:
        print("❌ No hay suficientes poblaciones válidas para comparar.")
        return

    window_starts = np.arange(positions.min(), positions.max(), window_size)
    windows = np.column_stack((window_starts, window_starts + window_size))

    all_fst_results = []
    labels = []

    for (i, j) in combinations(range(len(subpops)), 2):
        ac1 = genotypes[:, subpops[i]].count_alleles()
        ac2 = genotypes[:, subpops[j]].count_alleles()

        fst_vals = []
        for start, end in windows:
            mask = (positions >= start) & (positions < end)
            if np.sum(mask) < 2:
                fst_vals.append(np.nan)
                continue
            num, den = allel.hudson_fst(ac1[mask], ac2[mask])
            with np.errstate(invalid='ignore', divide='ignore'):
                fst_window = np.sum(num) / np.sum(den)
            fst_vals.append(fst_window)

        par_label = f"{nombres[i]} vs {nombres[j]}"
        labels.append(par_label)
        all_fst_results.append(fst_vals)

        df_result = pd.DataFrame({
            "start": windows[:, 0],
            "end": windows[:, 1],
            "fst": fst_vals
        })
        df_result.to_csv(os.path.join(output_dir, f"{par_label.replace(' ', '_')}_fst_ventanas.csv"), index=False)

        # Gráfico individual por par
        plt.figure(figsize=(10, 4))
        plt.plot(windows[:, 0], fst_vals, marker='o', linestyle='-', markersize=2)
        plt.title(f"FST por ventana: {par_label}")
        plt.xlabel("Posición genómica")
        plt.ylabel("FST (Hudson)")
        plt.ylim(0, 1)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{par_label.replace(' ', '_')}_fst_linea.png"))
        plt.close()

    # 🎨 Gráfico global comparativo - líneas
    fig, axs = plt.subplots(len(all_fst_results), 1, figsize=(10, 3.5 * len(all_fst_results)), sharex=True)
    if len(all_fst_results) == 1:
        axs = [axs]
    for i, (label, valores) in enumerate(zip(labels, all_fst_results)):
        axs[i].plot(windows[:, 0], valores, marker='o', linestyle='-', markersize=2)
        axs[i].set_title(label)
        axs[i].set_ylabel("FST")
        axs[i].set_ylim(0, 1)
    axs[-1].set_xlabel("Posición genómica")
    plt.suptitle("Comparación de FST por ventana entre pares de poblaciones", y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fst_lineas_comparativas.png"), bbox_inches="tight")
    plt.close()

    # 🎨 Gráfico global comparativo - boxplots
    plt.figure(figsize=(10, 5))
    plt.boxplot([pd.Series(vals).dropna() for vals in all_fst_results], labels=labels, patch_artist=True)
    plt.ylabel("FST")
    plt.title("Distribución de FST por pares de poblaciones")
    plt.ylim(0, 1)
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fst_boxplots_comparativos.png"))
    plt.close()

    # 🧾 Tabla resumen de FST promedio por comparación
    resumen = []
    for label, valores in zip(labels, all_fst_results):
        media = pd.Series(valores).mean(skipna=True)
        resumen.append({'Comparacion': label, 'FST_promedio': round(media, 5)})

    df_resumen = pd.DataFrame(resumen).sort_values(by='FST_promedio', ascending=False)
    df_resumen.to_csv(os.path.join(output_dir, "tabla_fst_promedio.csv"), index=False)
    print("📊 Tabla resumen guardada: tabla_fst_promedio.csv")

    # 🔥 Heatmap de FST promedio entre pares
    try:
        import seaborn as sns

        # Crear matriz cuadrada
        poblaciones_unicas = sorted(set(sum([label.split(" vs ") for label in labels], [])))
        matriz = pd.DataFrame(np.nan, index=poblaciones_unicas, columns=poblaciones_unicas)

        for label, valores in zip(labels, all_fst_results):
            pop1, pop2 = label.split(" vs ")
            media = pd.Series(valores).mean(skipna=True)
            matriz.loc[pop1, pop2] = media
            matriz.loc[pop2, pop1] = media  # Simétrico

        plt.figure(figsize=(8, 6))
        sns.heatmap(matriz, annot=True, fmt=".3f", cmap="YlGnBu", vmin=0, vmax=1)
        plt.title("Heatmap de FST promedio entre poblaciones")
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "fst_promedio_heatmap.png"))
        plt.close()

        print("🧊 Heatmap de FST promedio guardado: fst_promedio_heatmap.png")
    except ImportError:
        print("⚠️  No se pudo generar el heatmap (falta seaborn)")


    print("✅ FST por ventanas calculado y visualizado para todos los pares.")

if __name__ == "__main__":
    calcular_fst_por_pares("ruta/a/archivo.vcf", "ruta/a/poblaciones.csv", "salida/fst_por_pares")

