from collections import defaultdict
import os, csv, statistics, itertools
from Sets_of_Enzymes_exp import *

def read_extracted_files(filepath):
  dir = "data/"
  path = os.path.join(dir,filepath)
  data = defaultdict(list)
  with open(path, "r") as file:
    for row in csv.DictReader(file, delimiter="\t"):
        row["Median"] = ""
        recall = float(row["Real Set Recall"])
        sd = float(row["SD"])
        mean = float(row["Mean"])
        if sd < 0.1:
          zscore = (recall - mean) / sd
          row["Z Score"] = zscore
        
        

        
        data[(row["Group Name"], row["Tissue"])].append(row)
        
  return data 


def combine_files(d_lis, filepath):
  dir = "data/"
  path = os.path.join(dir,filepath)
  with open(path, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    head = ["Group Name", "Gene", "Tissue", "Z Score", "Real Set Recall", "Mean", "SD", "Old SD"] + [f"Set {i + 1}" for i in range(20)] + ["Median"]
    writer.writerow(head)
    for d in d_lis: 
      for keys, values in d.items(): 
        for l in values:
          writer.writerow([v for v in l.values()])


#d1 = read_extracted_files("z_table_csv_3.tsv")
#d2 = read_extracted_files("z_table_csv_4.tsv")
#d3 = read_extracted_files("z_table_csv_5.tsv")
#d4 = read_extracted_files("z_table_csv_6.tsv")
#d5 = read_extracted_files("z_table_csv_7.tsv")


#df_lis = [d1,d2,d3,d4,d5]
filename = "combine_5_exp.tsv"
#combine_files(df_lis, filename)

def add_median(data):
  m_data = defaultdict(list)
  for key, values in data.items():
    track_score = []
    for v in values:
        track_score.append((float(v.get("Z Score"))))
    
    srt_lis = sorted(track_score)
    median = statistics.median(srt_lis)
    for v in values:
        if float(v.get("Z Score")) == median:
            v["Median"] = "X"

        m_data[key].append(v)      
  
  return m_data

u_filename = "Updated.tsv"
df = read_extracted_files(filename)
c_d = add_median(df)
c_d = [c_d]
combine_files(c_d, u_filename)

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
        if zscore > 5 and recall  > 0.5 and median == "X":
          lis = []
          for enz in v.get("Gene").split(" "):
            if enz != "":
              lis.append(enz)
          possible_pairs = itertools.combinations(lis, 2)
          for p in possible_pairs:
            pair_set[i] = p
            i += 1

    return pair_set


g = Groups(c_d)
singles = (g.singleton())
sin_pair_set = (g.pairs())

print(singles)

last_group_number = list(sin_pair_set)[-1]
for number, gn in singles.items():
  last_group_number += 1
  sin_pair_set.update({last_group_number: gn})
  new_n = list(sin_pair_set)[-1]


"""
total_gr_zscore = defaultdict(list)
g = GeneSet()
glyco_enzymes = g.get_all_glyco_enz()
for group_number, gn_set in sin_pair_set.items():
  #print(group_number, gn_set)
  if len(gn_set) < 2:
    if gn_set != {"GALNT17"}: 
      print(group_number, gn_set)
      l_g = list(gn_set)[0]
      total_gr_zscore[(group_number, l_g)] = set_enz_experiment(gn_set, glyco_enzymes, extracted_dataset)
  else:
    if "GALNT17" not in gn_set:
      print(group_number, gn_set)
      gn_group_set = ""
      for gn in gn_set:
        gn_group_set += f"{gn} "
      total_gr_zscore[(group_number, gn_group_set)] = set_enz_experiment(gn_set, glyco_enzymes, extracted_dataset)

save_zdata(total_gr_zscore, "Single_pair") 
"""