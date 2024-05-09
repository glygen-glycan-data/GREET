from GTExData import *
from collections import defaultdict
import csv, gzip, urllib.request, json, io, random
import numpy as np
import os, sys, re


class GeneSet():
    sandbox_url = "https://edwardslab.bmcb.georgetown.edu/sandboxdev/api/getEnzymeMappings.php?limiter=no_filter&val="   
    
    def __init__(self):      
      self.gen_sdbox_data()

    def gen_sdbox_data(self):
      with urllib.request.urlopen(self.sandbox_url) as response:
        json_data = response.read()
        json_data = json_data.decode()
        self.sandbox_data = json.loads(json_data)
      
      self.total_enzymes = set()
      ext_data = self.sandbox_data["data"]
      self.enzymes_dict_set = defaultdict(set)
      for row in self.sandbox_data["data"]:
        if row['species'] == 'Homo sapiens' and row['gene_name'] != None:
          self.total_enzymes.add(row['gene_name'])
        for e in ext_data:
          if row['anomer'] == e['anomer'] and row['form_name'] == e['form_name'] and row['site'] == e['site']:
            if row['species'] == "Homo sapiens" and row['gene_name'] != None:
              anomer_site = row["anomer"] + row["site"]
              sd_group_name = f"{row['form_name']}-{anomer_site}-{row['parent_form_name']}"
              self.enzymes_dict_set[sd_group_name].add(row['gene_name'])
              
   
    def get_sdbox_data(self):
      sdbox_data = {}
      for k, gen_set in self.enzymes_dict_set.items():
        sdbox_data[tuple(sorted(gen_set))] = k 
      
      return sdbox_data

    def get_all_glyco_enz(self):   
      return self.total_enzymes
    
    def extract_glyco_set_at(self, group_name):   
      group_df = self.enzymes_dict_set.get(group_name)
      return  group_df
    
x = GeneSet()
x.get_sdbox_data()


def samples(url):
    samples_names = defaultdict(list)
    with urllib.request.urlopen(url) as response:
        decoded_response = response.read().decode('utf-8')
        reader = csv.DictReader(decoded_response.splitlines(), delimiter="\t")
        for row in reader:
            samples_names[row["SMTSD"]].append(row["SAMPID"])

    return samples_names




def extracting_data(extract_enzymes_tup, samples_names):
  dataset = ["prostate", "whole_blood", "vagina", "uterus", "thyroid", "testis", "stomach", "spleen", "small_intestine_terminal_ileum", "skin_sun_exposed_lower_leg", "skin_not_sun_exposed_suprapubic", "pancreas", "pituitary", "ovary", "nerve_tibial", "muscle_skeletal", "minor_salivary_gland",
           "lung", "kidney_cortex", "kidney_medulla", "heart_left_ventricle", "heart_atrial_appendage", "fallopian_tube", "esophagus_gastroesophageal_junction", "esophagus_mucosa", "esophagus_muscularis", "colon_sigmoid", "colon_transverse",
           "cervix_ectocervix", "cervix_endocervix", "cells_cultured_fibroblasts", "cells_ebv-transformed_lymphocytes", "breast_mammary_tissue", "brain_amygdala", "brain_anterior_cingulate_cortex_ba24", "brain_caudate_basal_ganglia", "brain_cerebellar_hemisphere", "brain_cerebellum", "brain_cortex",
           "brain_frontal_cortex_ba9", "brain_substantia_nigra", "brain_hippocampus", "brain_hypothalamus", "brain_nucleus_accumbens_basal_ganglia", "brain_putamen_basal_ganglia", "brain_spinal_cord_cervical_c-1", "bladder", "artery_aorta", "artery_coronary", "artery_tibial", "adipose_subcutaneous",
           "adipose_visceral_omentum", "adrenal_gland", "liver"]
  
  tissue_sets_data = defaultdict(list)
  for d in dataset:
    print(d + ".gct.gz")
    non_glyco_enzymes = set()
    gtex_url = "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/tpms-by-tissue/gene_tpm_2017-06-05_v8_" + d + ".gct.gz"
    with urllib.request.urlopen(gtex_url) as response:
      with gzip.open(io.BytesIO(response.read()), "rb") as file:
        header_lines = [next(file).decode().strip() for _ in range(3)]
        column_names = header_lines[2].split('\t')
        for line in file:
          line = line.decode()
          values = list(csv.reader([line], delimiter='\t'))[0]
          row_dict = dict(zip(column_names, values))
          gene_name = row_dict.get("Description")
          if gene_name in extract_enzymes_tup or len(non_glyco_enzymes) <= 5000:
            if gene_name not in extract_enzymes_tup:
              non_glyco_enzymes.add(gene_name)
            for tissue, samp in samples_names.items():
              tis_dict = defaultdict(list)
              for s in samp:
                s_get = row_dict.get(s)
                if s_get != None:
                  if tissue != None:
                    tis_dict[tissue].append(float(s_get))
              if tis_dict != {} :
                tissue_sets_data[gene_name].append(tis_dict)


  return tissue_sets_data, non_glyco_enzymes




def setting_values_per_tissue(dataset, total_sample=False): 
  tissue_set = {}
  tissue_type_dict = defaultdict(list)
  for k, v_dic in dataset.items():
    for v_lis in v_dic:
      if v_lis != None:
        if isinstance(v_lis, list):
          for val in v_lis:  
            for tis, values in val.items():
              tissue_set[tis] = len(values)
              for v in values:
                tissue_type_dict[k].append(v)
        else:
          for tis, values in v_lis.items():
            tissue_set[tis] = len(values)
            for v in values:
              tissue_type_dict[k].append(v)

  if total_sample == True:
    total_num_sample = 0
    for values in tissue_set.values():
      total_num_sample += values

  return tissue_type_dict, tissue_set



def adding_headers(tissue_dict_dataset, collapsed=False):
  headers = []
  collapsed_headers = []
  for key, val in tissue_dict_dataset.items():
    for i in range(val):
        headers.append(key)
    if collapsed == True:
      split_head = key.split('-')[0]
      if split_head != "Brain ":
        for i in range(val):
          collapsed_headers.append(split_head)
      else:
        for i in range(val):
          collapsed_headers.append(key)
  if collapsed == True:
    return collapsed_headers
  else:
    return headers
  
def assign_class(tissue_headers, tissue_name):
  clas_nm = []
  for l in tissue_headers:
    if l == tissue_name:
      ct = 1
      clas_nm.append(ct)

    elif l == tissue_name.split('-')[0]:
      ct = 1
      clas_nm.append(ct)

    else:
      ct = 0
      clas_nm.append(ct)

  clas_nm = np.array(clas_nm)

  return clas_nm


def reading_file_in(path):
  tis_names = set()
  temp_data = defaultdict(list)
  data = {}
  with open(path, "r") as file:
      dict_reader = csv.DictReader(file, delimiter="\t")
      for row in dict_reader:
        temp_data[row["Gene"]].append((row["Tissue"], float(row["Value"])))
  for gn, tis_val in temp_data.items():
    tissue_collection = defaultdict(list)
    for t_v in tis_val:
      tissue_collection[t_v[0]].append(t_v[1])
  
  
    total_num_sample = 0
    reduced_tissue_collection = {}
    for tissue, values in tissue_collection.items():
      if len(values) >= 100:
        reduced_tissue_collection[tissue] = values
        print(tissue, len(values))
        tis_names.add(tissue)
        total_num_sample += len(values)
      else:
        print(tissue, len(values))

    print(len(tis_names))
    print(total_num_sample)
    if total_num_sample == 17382:
      data[gn] = tissue_collection
  
  return data, sorted(tis_names)

      
def save_file_to(data, filename, delimiter='\t'):
    tis_names = set()
    full_path = os.path.join("data/" + filename)
    with open(full_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=delimiter)
        head = ["Gene", "Tissue", "Value"]
        writer.writerow(head)
        header = list(data.keys())
        for enz in header:
          for t in data.get(enz):
            for k_tis,val_lis in t.items():
              tis_names.add(k_tis)
              for v in val_lis:
                writer.writerow([enz, k_tis, v])
    
    return sorted(tis_names)



def check_for_file():
  file_path = "data/temp_data.csv"
  sam_url = "https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
  
  if os.path.exists(file_path):
    extracted_dataset, tissues_names = reading_file_in(file_path)

  else:
    samples_names = samples(sam_url)
    gn_set = GeneSet()
    total_enzymes = gn_set.get_all_glyco_enz()
    extracted_dataset, non_genes = extracting_data(total_enzymes, samples_names)
    tissues_names = save_file_to(extracted_dataset, "temp_data.csv")

  return extracted_dataset, tissues_names


extracted_dataset, tissues_names = check_for_file()


#genes = set(sys.argv[1:])
#gtd = GTexData()
#data = gtd.restrict(genes=genes,nottissue='Liver',includetissue=True)

#print(data)


