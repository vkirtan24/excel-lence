import os
import pandas as pd

def count_variant_metrics(vcf_file):
    total_variant_count = 0
    passed_variant_count = 0
    snp_count = 0
    insertion_count = 0
    deletion_count = 0
    complex_indel_count = 0
    mixed_count = 0
    het_count = 0
    hom_var_count = 0
    singleton_count = 0
    called_genotype_count = 0

    with open(vcf_file, 'r') as file:
        for line in file:
            if line.startswith('##'):
                continue

            if line.startswith('##contig'):
                continue  

            fields = line.strip().split('\t')

            if len(fields) < 10:
                continue  

            alt_alleles = fields[4].split(',')

            total_variant_count += 1

            if 'PASS' in fields[6]:
                passed_variant_count += 1

            if len(alt_alleles) == 1 and len(alt_alleles[0]) == 1:
                snp_count += 1
            elif len(alt_alleles) == 1 and len(alt_alleles[0]) > 1:
                insertion_count += 1
            elif len(alt_alleles) > 1 and any(len(alt) > 1 for alt in alt_alleles):
                complex_indel_count += 1
            elif len(alt_alleles) == 1 and alt_alleles[0] == '*':
                deletion_count += 1
            else:
                mixed_count += 1

            genotype_info = fields[9].split(':')
            gt_field = genotype_info[0]

            if gt_field == '0/1':
                het_count += 1
            elif gt_field == '1/1':
                hom_var_count += 1

            if gt_field != './.':
                called_genotype_count += 1

            if len(alt_alleles) == 1 and alt_alleles[0] != '*':
                if gt_field == '0/1' or gt_field == '1/1':
                    singleton_count += 1

    return {
        'Sample Name': os.path.splitext(os.path.basename(vcf_file))[0],
        'Number of total variant loci': total_variant_count,
        'Number of passed variant loci': passed_variant_count,
        'Number of SNP loci': snp_count,
        'Number of insertions': insertion_count,
        'Number of deletions': deletion_count,
        'Number of complex indels': complex_indel_count,
        'Number of mixed loci': mixed_count,
        'Number of het loci': het_count,
        'Number of hom var loci': hom_var_count,
        'Number of singletons': singleton_count,
        'Number of called genotypes': called_genotype_count
    }

directory = '' #Here add the path of the directory in which you contain the Multiple VCF file
vcf_dir = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith('.vcf')]

results = []
for file in vcf_dir:
    results.append(count_variant_metrics(file))

df = pd.DataFrame(results)

output_excel_path = os.path.join(directory, 'var.xlsx')
df.to_excel(output_excel_path, index=False)

print(f'Output saved in {output_excel_path}')
