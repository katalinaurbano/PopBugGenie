import gzip
import os

def limpiar_vcf_para_scikit_allel(vcf_path, reemplazar_missing=False):
    """
    Limpia un archivo VCF para facilitar su lectura con scikit-allel.
    - Si reemplazar_missing=True, reemplaza campos ./., ./.:.:.:... por una cadena consistente con FORMAT.
    - También corrige casos raros como '0' → '0/0', '1' → '1/1'
    - Guarda el resultado como archivo comprimido .vcf.gz
    """
    vcf_out = vcf_path.replace('.vcf.gz', '_preprocesado.vcf.gz').replace('.vcf', '_preprocesado.vcf.gz')
    abrir = gzip.open if vcf_path.endswith('.gz') else open
    modo_lectura = 'rt' if vcf_path.endswith('.gz') else 'r'

    with abrir(vcf_path, modo_lectura) as infile, gzip.open(vcf_out, 'wt') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
            else:
                cols = line.strip().split('\t')
                format_fields = cols[8].split(':')

                # Construir el reemplazo dinámico
                reemplazo = '0/0'
                for campo in format_fields[1:]:
                    if campo in ['GQ', 'PL']:
                        reemplazo += ':30'
                    else:
                        reemplazo += ':1'

                fixed = cols[:9]
                for sample in cols[9:]:
                    # Casos ausentes: ./., ./.:.:.:.:...
                    if reemplazar_missing and (sample == './.' or sample.startswith('./.:')):
                        fixed.append(reemplazo)
                    # Casos raros: solo '0' o '1'
                    elif sample == '0':
                        fixed.append('0/0')
                    elif sample == '1':
                        fixed.append('1/1')
                    else:
                        fixed.append(sample)

                outfile.write('\t'.join(fixed) + '\n')

    print(f"✅ Archivo preprocesado guardado en: {vcf_out}")
    return vcf_out

