import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def graficar_heatmaps_combinados(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    # Ahora indicamos archivo y nombre real de columna
    archivo_csvs = {
        'pi_zscore': ('pi.csv', 'pi_zscore'),
        'tajima_zscore': ('tajima_d.csv', 'tajima_d_zscore'),
        'ld_zscore': ('ld.csv', 'ld_r2_zscore'),
        'hd_zscore': ('hd_por_ventana.csv', 'hd_zscore')

    }

    estadisticos = list(archivo_csvs.keys())
    poblaciones = [d for d in os.listdir(input_dir)
                   if os.path.isdir(os.path.join(input_dir, d))]
    poblaciones = sorted(poblaciones)

    for pop in poblaciones:
        datos = []
        posiciones = None

        for est, (archivo, columna) in archivo_csvs.items():
            path = os.path.join(input_dir, pop, archivo)
            if not os.path.exists(path):
                continue

            df = pd.read_csv(path)
            if 'start' not in df.columns or columna not in df.columns:
                continue

            if posiciones is None:
                posiciones = df['start']

            datos.append(df[columna].values)

        if datos and posiciones is not None:
            matriz = pd.DataFrame(datos, index=estadisticos, columns=posiciones)
            plt.figure(figsize=(14, 3))
            sns.heatmap(matriz, cmap="coolwarm", center=0, cbar_kws={'label': 'Z-score'})
            plt.title(f"Heatmap combinado de Z-scores - {pop}")
            plt.xlabel("Posición genómica (start)")
            plt.ylabel("Estadístico")
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f"{pop}_heatmap_combinado.png"))
            plt.close()

if __name__ == "__main__":
    graficar_heatmaps_combinados("resultados_por_poblacion", "heatmaps_combinados")

