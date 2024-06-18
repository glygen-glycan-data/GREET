from collections import defaultdict
import os, csv, itertools
from Sets_of_Enzymes_exp import *
from preprocessing_setup import HPATissueClassifcation

def read_extracted_files(filepath):
  dir = "data/"
  path = os.path.join(dir,filepath)
  data = defaultdict(list)
  with open(path, "r") as file:
    for row in csv.DictReader(file, delimiter="\t"):
      data[(row["Group Name"], row["Tissue"])].append(row)
  return data 


def combine_files(d_lis, filepath):
  dir = "data/"
  path = os.path.join(dir,filepath)
  with open(path, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    head = ["Group Name", "Gene", "Tissue", "Z Score", "Real Set Recall", "Mean", "SD", "Old SD"] + [f"Set {i + 1}" for i in range(20)]  + ["Median"] + ["HPA Tissue", "Tissue Specificity"] #+ ["HPA Tissue", "Tissue Specificity"]
    writer.writerow(head)
    for d in d_lis: 
      for keys, values in d.items(): 
        for l in values:
          writer.writerow([v for v in l.values()])



def add_median(data):
  for values in data.values():
    srt_rows = sorted(values, key=lambda x:x["Z Score"])
    for i in range(len(srt_rows)):
      if i == 2:
        median_row = srt_rows[2]
        median_row["Median"] = "X"
      else:
        
        srt_rows[i]["Median"] = ""
    
  
  return data


class Groups:
  def __init__(self, data):
    self.data = data

  def singleton(self):
    single_sets = defaultdict(set)
    tract_genes = set()
    i = 1
    for key, values in self.data.items():
      for v in values:
        for enz in v.get("Gene").split(" "):
          if enz != "" and enz not in tract_genes:
            enz_dne = ["GGTA1", "GALNT17", "POGLUT2", "POGLUT3"]
            if  enz not in enz_dne:
              single_sets[i].add(enz)
              i += 1
              tract_genes.add(enz)    

    return single_sets
  
  def pairs(self):
    pair_set = {}
    i = 1
    for key, values in self.data.items():
      for v in values:
        zscore = float(v.get("Z Score"))
        recall = float(v.get("Real Set Recall"))
        median = v.get("Median")
        if zscore > 5 and recall  > 0.5:# and median == "X":
          lis = []
          track_sets = set()
          for enz in v.get("Gene").split(" "):
            if enz != "" and enz != "GALNT17":
              lis.append(enz)

          possible_pairs = itertools.combinations(lis, 2)
          for p in possible_pairs:
            if p not in track_sets:
              pair_set[i] = p
              i += 1
              track_sets.add(p)

    return pair_set

def single_pair_run(enz_sets, glyco_enz, extract_df, filename):
  total_gr_zscore = defaultdict(list)
  for group_number, gn_set in enz_sets.items():
   if len(gn_set) < 2:
      print(group_number, gn_set)
      l_g = list(gn_set)[0]
      total_gr_zscore[(group_number, l_g)] = set_enz_experiment(gn_set, glyco_enz, extract_df)
  
   else:
      print(group_number, gn_set)
      gn_group_set = ""
      for gn in gn_set:
        gn_group_set += f"{gn} "
      total_gr_zscore[(group_number, gn_group_set)] = set_enz_experiment(gn_set, glyco_enz, extract_df)

  save_zdata(total_gr_zscore, filename)



filename = "combine_U_HPA.tsv"
df = read_extracted_files(filename)
print(df)

"""
#Run Each experiment with 5 times with different file name, and read them in seperate
#d1 = read_extracted_files("z_table_csv_Logistic_HPA.tsv")
#d2 = read_extracted_files("z_table_csv_Logistic_HPA_2.tsv")
#d3 = read_extracted_files("z_table_csv_Logistic_HPA_3.tsv")
#d4 = read_extracted_files("z_table_csv_Logistic_HPA_4.tsv")
#d5 = read_extracted_files("z_table_csv_Logistic_HPA_5.tsv")
   
#df_lis = [d1,d2,d3,d4,d5]
filename = "combine_HPA.tsv"
#combine_files(df_lis, filename)



#u_filename = "z_table_csv_Singles_Logistics_U.tsv"
df = read_extracted_files(filename)
d = add_median(df)
#combine_files([d], "combine_U_HPA.tsv")



#g = Groups(df)
#singles = (g.singleton())
#pair_set = (g.pairs())
#print(singles)
#s = single_pair_run(singles, glyco_enzymes, extracted_dataset, "Logistic_HPA_5")
#pairs = single_pair_run(pair_set, glyco_enzymes, extracted_dataset, "Pair_set_L1")




#z_file = "z_table_csv_Pair_set.tsv"
#z_read = read_extracted_files(z_file)
#h = HPATissueClassifcation(HPATissue_file, Gspec_file)

data_Spec = h.extract_gene_specificity(d)
combine_files([data_Spec], "combine_U_HPA.tsv")

"""