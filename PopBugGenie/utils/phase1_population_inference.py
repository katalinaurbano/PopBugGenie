import allel
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

def detectar_poblaciones_automaticamente(vcf_path, output_dir, max_clusters=8):
    os.makedirs(output_dir, exist_ok=True)

    print("\n[FASE 1] Cargando VCF...")
    callset = allel.read_vcf(vcf_path)
    genotypes = allel.GenotypeArray(callset['calldata/GT'])
    samples = callset['samples']

    print("-> Genotipos cargados:", genotypes.shape)

    print("-> Filtrando variantes con datos faltantes...")
    gn_alt = genotypes.to_n_alt()
    gn_alt = gn_alt[~np.any(np.isnan(gn_alt), axis=1)]
    print("-> Forma después del filtrado:", gn_alt.shape)

    print("-> Calculando PCA...")
    gn_transposed = gn_alt.T
    scaler = StandardScaler()
    gn_scaled = scaler.fit_transform(gn_transposed)

    pca = PCA(n_components=2)
    coords = pca.fit_transform(gn_scaled)
    explained = pca.explained_variance_ratio_ * 100
    print(f"-> Varianza explicada: {explained[:2]}")

    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score

    print("-> Calculando número óptimo de clusters...")
    inertia = []
    silhouettes = []
    K = range(2, max_clusters + 1)

    for k in K:
        km = KMeans(n_clusters=k, random_state=42)
        labels = km.fit_predict(coords)
        inertia.append(km.inertia_)
        silhouettes.append(silhouette_score(coords, labels))

    plt.figure()
    plt.plot(K, inertia, marker='o')
    plt.xlabel('Número de clusters (k)')
    plt.ylabel('Inercia total')
    plt.title('Método del codo')
    plt.grid(True)
    elbow_path = os.path.join(output_dir, 'codo_clusters.png')
    plt.savefig(elbow_path)
    print(f"Gráfico del codo guardado en {elbow_path}")

    best_k = K[np.argmax(silhouettes)]
    print(f"Número óptimo de clusters según silhouette: {best_k}")

    kmeans = KMeans(n_clusters=best_k, random_state=42)
    cluster_labels = kmeans.fit_predict(coords)

    df = pd.DataFrame({
        'Sample': samples,
        'PC1': coords[:, 0],
        'PC2': coords[:, 1],
        'Cluster': cluster_labels
    })

    output_csv = os.path.join(output_dir, 'poblaciones_auto.csv')
    df.to_csv(output_csv, index=False)
    print(f"✅ Poblaciones inferidas guardadas en {output_csv}")

    plt.figure(figsize=(10, 6))
    for cluster in sorted(df['Cluster'].unique()):
        sub = df[df['Cluster'] == cluster]
        plt.scatter(sub['PC1'], sub['PC2'], label=f'Cluster {cluster}', alpha=0.7)
    plt.title(f'PCA + Clustering automático (k={best_k})')
    plt.xlabel(f'PC1 ({explained[0]:.1f}%)')
    plt.ylabel(f'PC2 ({explained[1]:.1f}%)')
    plt.legend()
    plt.tight_layout()
    pca_path = os.path.join(output_dir, 'pca_clusters.png')
    plt.savefig(pca_path)
    print(f"Gráfico PCA guardado en {pca_path}")

    return df

