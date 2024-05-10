import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from tkinter import *
from tkinter import filedialog, messagebox

def extract_variant_metrics(vcf_file):
    variant_metrics = {
        'Sample Name': os.path.splitext(os.path.basename(vcf_file))[0],
        'Total Variant Count': 0,
        'Passed Variant Count': 0,
        'SNP Count': 0,
        'Insertion Count': 0,
        'Deletion Count': 0,
        'Complex Indel Count': 0,
        'Mixed Count': 0,
        'Het Count': 0,
        'Hom Var Count': 0,
        'Singleton Count': 0,
        'Called Genotype Count': 0
    }

    with open(vcf_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')

            variant_metrics['Total Variant Count'] += 1

            if 'PASS' in fields[6]:
                variant_metrics['Passed Variant Count'] += 1

            alt_alleles = fields[4].split(',')
            gt_field = fields[9].split(':')[0]

            if len(alt_alleles) == 1 and len(alt_alleles[0]) == 1:
                variant_metrics['SNP Count'] += 1
            elif len(alt_alleles) == 1 and len(alt_alleles[0]) > 1:
                variant_metrics['Insertion Count'] += 1
            elif len(alt_alleles) > 1 and any(len(alt) > 1 for alt in alt_alleles):
                variant_metrics['Complex Indel Count'] += 1
            elif len(alt_alleles) == 1 and alt_alleles[0] == '*':
                variant_metrics['Deletion Count'] += 1
            else:
                variant_metrics['Mixed Count'] += 1

            if gt_field == '0/1':
                variant_metrics['Het Count'] += 1
            elif gt_field == '1/1':
                variant_metrics['Hom Var Count'] += 1

            if gt_field != './.':
                variant_metrics['Called Genotype Count'] += 1

            if len(alt_alleles) == 1 and alt_alleles[0] != '*':
                if gt_field == '0/1' or gt_field == '1/1':
                    variant_metrics['Singleton Count'] += 1

    return variant_metrics

def browse_files():
    global directory
    directory = filedialog.askdirectory()
    execute_script()

def execute_script():
    global directory
    if directory:
        all_metrics = []

        for vcf_file in os.listdir(directory):
            if vcf_file.endswith('.vcf'):
                vcf_file_path = os.path.join(directory, vcf_file)
                variant_metrics = extract_variant_metrics(vcf_file_path)
                all_metrics.append(variant_metrics)

        output_excel_path = os.path.join(directory, 'var.xlsx')
        pd.DataFrame(all_metrics).to_excel(output_excel_path, index=False)

        x_labels = list()
        data = list()
        colors = ['blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black', 'purple', 'orange', 'brown', 'pink']

        for vcf_file, color in zip(all_metrics, colors):
            sample_name = vcf_file['Sample Name']
            x_labels.extend([f"{sample_name}_{variant_type}" for variant_type in vcf_file.keys() if variant_type != 'Sample Name'])
            data.extend(list(vcf_file.values())[1:])  
        x_positions = np.arange(len(x_labels))

        plt.figure(figsize=(18, 9)) 

        plt.bar(x_positions, data, alpha=0.5, color=colors)

        for x, y in zip(x_positions, data):
            plt.text(x, y, str(y), ha='center', va='bottom', rotation=90)

        plt.xlabel('Sample_Variant Type')
        plt.ylabel('Number of Variants')
        plt.title('Number of each Variants for given Samples')
        plt.xticks(x_positions, x_labels, rotation=90)

        plt.tight_layout()  
        plt.savefig(os.path.join(directory, 'var_plot.png'))
        plt.show()

        messagebox.showinfo("Success", f"Variant metrics calculated and saved as '{output_excel_path}'.")

root = Tk()
root.title("Excel-lence")

image = PhotoImage(file="excel-lence.png")  
resized_image = image.subsample(2, 2)  
image_label = Label(root, image=resized_image)
image_label.pack()

select_directory_button = Button(root, text="Select Directory", command=browse_files)
select_directory_button.pack(pady=10)

root.geometry("850x700") 
root.mainloop()
