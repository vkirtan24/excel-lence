import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def extract_variant_metrics(vcf_file):
    variant_metrics = {
        'total_variant_count': 0,
        'passed_variant_count': 0,
        'snp_count': 0,
        'insertion_count': 0,
        'deletion_count': 0,
        'complex_indel_count': 0,
        'mixed_count': 0,
        'het_count': 0,
        'hom_var_count': 0,
        'singleton_count': 0,
        'called_genotype_count': 0
    }

    with open(vcf_file, 'r') as file:
        for line in file:
            if line.startswith('#'):  
                continue

            fields = line.strip().split('\t')

            variant_metrics['total_variant_count'] += 1

            if 'PASS' in fields[6]:
                variant_metrics['passed_variant_count'] += 1

            alt_alleles = fields[4].split(',')
            gt_field = fields[9].split(':')[0]

            if len(alt_alleles) == 1 and len(alt_alleles[0]) == 1:
                variant_metrics['snp_count'] += 1
            elif len(alt_alleles) == 1 and len(alt_alleles[0]) > 1:
                variant_metrics['insertion_count'] += 1
            elif len(alt_alleles) > 1 and any(len(alt) > 1 for alt in alt_alleles):
                variant_metrics['complex_indel_count'] += 1
            elif len(alt_alleles) == 1 and alt_alleles[0] == '*':
                variant_metrics['deletion_count'] += 1
            else:
                variant_metrics['mixed_count'] += 1

            if gt_field == '0/1':
                variant_metrics['het_count'] += 1
            elif gt_field == '1/1':
                variant_metrics['hom_var_count'] += 1

            if gt_field != './.':
                variant_metrics['called_genotype_count'] += 1

            if len(alt_alleles) == 1 and alt_alleles[0] != '*':
                if gt_field == '0/1' or gt_field == '1/1':
                    variant_metrics['singleton_count'] += 1

    return variant_metrics

directory = '' # Add the Directory path where you saved the multiple VCF files
vcf_dir = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith('.vcf')]

x_labels = list()
data = list()
colors = ['blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black', 'purple', 'orange', 'brown', 'pink']

for vcf_file, color in zip(vcf_dir, colors):
    sample_name = os.path.splitext(os.path.basename(vcf_file))[0]
    variant_metrics = extract_variant_metrics(vcf_file)
    x_labels.extend([f"{sample_name}_{variant_type}" for variant_type in variant_metrics.keys()])
    data.extend(list(variant_metrics.values()))

x_positions = np.arange(len(x_labels))

plt.figure(figsize=(12, 6))

plt.bar(x_positions, data, alpha=0.5, color=colors)

for x, y in zip(x_positions, data):
    plt.text(x, y, str(y), ha='center', va='bottom', rotation=90)

plt.xlabel('Sample_Variant Type')
plt.ylabel('Number of Variants')
plt.title('Number of each Variants for given Samples')
plt.xticks(x_positions, x_labels, rotation=90)

output_excel_path = os.path.join(directory, 'var.xlsx')
plt.savefig('var_plot.png')
plt.show()

print(f'Output saved in {output_excel_path}')
