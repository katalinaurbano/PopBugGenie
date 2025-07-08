import argparse
import os
import pandas as pd
import logging
from datetime import datetime

from utils.preprocesamiento import limpiar_vcf_para_scikit_allel
from utils.fase1_detectar_poblaciones import detectar_poblaciones_automaticamente
from utils.fase2_estadisticos_globales import calcular_estadisticos_globales
from utils.fase3_estadisticos_por_poblacion import calcular_estadisticos_por_poblacion
from utils.fase4_fst_hudson import calcular_fst_por_pares
from utils.graficas_comparativas import graficar_comparativas_por_poblacion
from utils.graficar_heatmaps_combinados import graficar_heatmaps_combinados


def main():
    import sys
    import shlex
    import time

    parser = argparse.ArgumentParser(description='Pipeline de análisis de diversidad genética en bacterias')
    parser.add_argument('--vcf', required=True, help='Ruta al archivo VCF de entrada')
    parser.add_argument('--output', required=True, help='Carpeta donde guardar los resultados')
    parser.add_argument('--max_clusters', type=int, default=8, help='Máximo número de clusters para evaluar')
    parser.add_argument('--window', type=int, default=50000, help='Tamaño de ventana para los cálculos por ventana')
    parser.add_argument('--replace-missing', action='store_true', help='Reemplazar genotipos ./., que son missing, por 0/0 (solo para bacterias clonales)')
    parser.add_argument('--min_samples_per_pop', type=int, default=3, help='Número mínimo de muestras por población para calcular los estadísticos')

    args = parser.parse_args()

    # Preparar log
    fecha = datetime.now().strftime('%Y%m%d_%H%M%S')
    os.makedirs(args.output, exist_ok=True)
    log_path = os.path.join(args.output, f"pipeline_log_{fecha}.log")

    logging.basicConfig(
        filename=log_path,
        filemode='w',
        level=logging.INFO,
        format='[%(asctime)s] %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    class DualLogger:
        def __init__(self):
            self.console = sys.stdout
        def write(self, message):
            self.console.write(message)
            logging.info(message.strip())
        def flush(self):
            self.console.flush()

    sys.stdout = DualLogger()
    sys.stderr = DualLogger()

    # Tiempo de inicio
    tiempo_inicio = time.time()
    comando = " ".join(shlex.quote(arg) for arg in sys.argv)

    print(f"\n🚀 INICIANDO PIPELINE - {fecha}")
    print(f"📝 Comando ejecutado: {comando}")
    print(f"VCF de entrada: {args.vcf}")
    print(f"Salida en: {args.output}")
    print(f"Ventana: {args.window} bp")
    print(f"Reemplazar missing: {'Sí' if args.replace_missing else 'No'}\n")


    print("[PREPROCESAMIENTO] Asegurando compatibilidad del VCF con scikit-allel...")
    if args.replace_missing:
        vcf_usable = limpiar_vcf_para_scikit_allel(args.vcf, reemplazar_missing=True)
    else:
        vcf_usable = args.vcf
    print(f"VCF que se usará en el análisis: {vcf_usable}")

    print("\n[FASE 1] Detección automática de poblaciones")
    detectar_poblaciones_automaticamente(
        vcf_path=vcf_usable,
        output_dir=args.output,
        max_clusters=args.max_clusters
    )

    print("\n[FASE 2] Estadísticos globales (normalizados)")
    calcular_estadisticos_globales(
        vcf_path=vcf_usable,
        output_dir=os.path.join(args.output, 'global'),
        window_size=args.window
    )

    print("\n[FASE 3] Estadísticos por población (normalizados)")
    calcular_estadisticos_por_poblacion(
        vcf_path=vcf_usable,
        poblaciones_csv=os.path.join(args.output, 'poblaciones_auto.csv'),
        output_dir=os.path.join(args.output, 'por_poblacion'),
        window_size=args.window,
        min_samples=args.min_samples_per_pop
    )

    print("\n[FASE 4] Comparación de FST por pares con método de Hudson")
    calcular_fst_por_pares(
        vcf_path=vcf_usable,
        poblaciones_csv=os.path.join(args.output, 'poblaciones_auto.csv'),
        output_dir=os.path.join(args.output, 'fst_por_pares'),
        window_size=args.window,
        min_samples=args.min_samples_per_pop
    )

    print("\n[VISUALIZACIÓN] Generando gráficas comparativas (normalizadas)")
    graficar_comparativas_por_poblacion(
        input_dir=os.path.join(args.output, 'por_poblacion'),
        output_dir=os.path.join(args.output, 'graficas_comparativas')
    )

    print("\n[RESUMEN] Generando tabla comparativa de Z-scores por población...")
    base_dir = os.path.join(args.output, 'por_poblacion')
    all_data = []
    for pop in sorted(os.listdir(base_dir)):
        pop_path = os.path.join(base_dir, pop)
        if not os.path.isdir(pop_path):
            continue
        try:
            df_pi = pd.read_csv(os.path.join(pop_path, 'pi.csv'))[['start', 'pi_zscore']]
            df_tajima = pd.read_csv(os.path.join(pop_path, 'tajima_d.csv'))[['tajima_d_zscore']]
            df_ld = pd.read_csv(os.path.join(pop_path, 'ld.csv'))[['ld_r2_zscore']]
            df_merged = pd.concat([df_pi, df_tajima, df_ld], axis=1)
            df_merged['Poblacion'] = pop
            all_data.append(df_merged)
        except Exception as e:
            print(f"⚠️  Error procesando {pop}: {e}")
            continue

    if all_data:
        df_final = pd.concat(all_data, ignore_index=True)
        df_final.to_csv(os.path.join(args.output, 'tabla_comparativa_zscores.csv'), index=False)
        print("✅ Tabla comparativa guardada en tabla_comparativa_zscores.csv")

        print("\n[VISUALIZACIÓN] Generando heatmaps combinados por población...")
        graficar_heatmaps_combinados(
            input_dir=os.path.join(args.output, 'por_poblacion'),
            output_dir=os.path.join(args.output, 'graficas_comparativas', 'heatmaps_combinados')
        )
    else:
        print("⚠️  No se pudo generar la tabla comparativa de Z-scores.")

    # Medir y mostrar tiempo total de ejecución
    tiempo_total = time.time() - tiempo_inicio
    minutos = int(tiempo_total // 60)
    segundos = int(tiempo_total % 60)
    print(f"\n⏱️ Duración total: {minutos} min {segundos} s")
    print("✅ Pipeline completado con éxito.")

if __name__ == '__main__':
    main()

