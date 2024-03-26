from collections import defaultdict
import csv, gzip, urllib.request, json, io, random
import numpy as np




def samples(url):
    samples_names = defaultdict(list)
    with urllib.request.urlopen(url) as response:
        decoded_response = response.read().decode('utf-8')
        reader = csv.DictReader(decoded_response.splitlines(), delimiter="\t")
        for row in reader:
            samples_names[row["SMTSD"]].append(row["SAMPID"])

    return samples_names



def sandbox_url_data(sandbox_url):
  with urllib.request.urlopen(sandbox_url) as response:
    json_data = response.read()
    json_data = json_data.decode()
    json_data = json.loads(json_data)


  total_enzymes = set() 
  sandbox_data = json_data["data"]
  ext_data = json_data["data"]
  enzymes_dict_set = defaultdict(set)
  for row in sandbox_data:
    if row['gene_name'] != None:
      total_enzymes.add(row['gene_name'])
    for e in ext_data:
      if row['anomer'] == e['anomer'] and row['form_name'] == e['form_name'] and row['site'] == e['site']:
        if row['species'] == "Homo sapiens" and row['gene_name'] != None:
          #key = re.split(r'\d+$', row['gene_name'])
          enzymes_dict_set[(row['form_name'], row['site'], row["anomer"], row['parent_form_name'])].add(row['gene_name'])
   
  return enzymes_dict_set, total_enzymes




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
          if gene_name in extract_enzymes_tup or len(non_glyco_enzymes) <= 100:
            if gene_name not in extract_enzymes_tup and '.' not in gene_name:
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
        for val in v_lis:
          for tis, v in val.items():
            tissue_set[tis] = len(v)
            for val in v:
              tissue_type_dict[k].append(val)

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
  print(tissue_name)
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


def creating_random_sets(gen_set, data, non_glyco):
  sample_enz_set = data.get(gen_set)
  enz_set = set()
  r_gen = defaultdict(set)
  for values in data.values():
    for v in values:
      if v not in sample_enz_set:
        enz_set.add(v)

  enz_len = len(sample_enz_set)
  for r_i in range(20):
    get_random_enz = random.sample(list(non_glyco), enz_len)
    r_gen[r_i + 1] = get_random_enz

  r_gen[21] = sample_enz_set

  return r_gen


