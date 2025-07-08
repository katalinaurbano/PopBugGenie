import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def graficar_comparativas_por_poblacion(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    estadisticos = {
        'pi': ('pi.csv', 'pi_zscore', 'Diversidad nucleotídica (π) - Z-score'),
        'tajima': ('tajima.csv', 'tajima_zscore', "Tajima's D - Z-score"),
        'ld': ('ld.csv', 'ld_zscore', 'Desequilibrio de ligamiento (r²) - Z-score')
    }

    percentil_superior = 2.326  # Z-score del percentil 99%

    poblaciones = [d for d in os.listdir(input_dir)
                   if os.path.isdir(os.path.join(input_dir, d))]
    poblaciones = sorted(poblaciones)

    for key, (filename, colname, titulo_general) in estadisticos.items():
        fig, axs = plt.subplots(len(poblaciones), 1, figsize=(10, 4 * len(poblaciones)), sharex=True)

        if len(poblaciones) == 1:
            axs = [axs]  # asegurar lista si solo una población

        for i, pop in enumerate(poblaciones):
            path = os.path.join(input_dir, pop, filename)
            if not os.path.exists(path):
                axs[i].set_title(f'{pop} (sin datos)')
                axs[i].axis('off')
                continue

            df = pd.read_csv(path)
            if 'start' not in df.columns or colname not in df.columns:
                axs[i].set_title(f'{pop} (columnas faltantes)')
                axs[i].axis('off')
                continue

            axs[i].plot(df['start'], df[colname], label=key.upper())
            axs[i].axhline(y=percentil_superior, color='red', linestyle='--', label='Z = +2.326')
            axs[i].set_title(pop)
            axs[i].set_ylabel('Z-score')
            axs[i].legend()

        plt.tight_layout()
        plt.suptitle(titulo_general, fontsize=16, y=1.02)
        plt.savefig(os.path.join(output_dir, f'{key}_comparativa.png'), bbox_inches='tight')
        plt.close()

