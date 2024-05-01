import os
import pandas as pd
import matplotlib.pyplot as plt

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
            if line.startswith('#'):  # Skip header lines
                continue

            fields = line.strip().split('\t')

            # Update variant metrics
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

directory = 'E:/vcf_file/'
vcf_dir = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith('.vcf')]

results = []
for file in vcf_dir:
    results.append(extract_variant_metrics(file))

df = pd.DataFrame(results)

plt.figure(figsize=(12, 6))

variant_types = list(df.columns)

for i, variant_type in enumerate(variant_types):
    plt.bar(variant_type, df[variant_type], label=variant_type)
    for j, value in enumerate(df[variant_type]):
        plt.text(i, value, str(value), ha='center', va='bottom')

plt.xlabel('Variant Type')
plt.ylabel('Number of Variants')
plt.title('Number of each Variants for given Sample')
plt.xticks(rotation=90)

output_excel_path = os.path.join(directory, 'var.xlsx')
df.to_excel(output_excel_path, index=False)

print(f'Output saved in {output_excel_path}')
