from collections import defaultdict
import os, csv, statistics

def read_extracted_files(filepath):
  dir = "data/"
  path = os.path.join(dir,filepath)
  data = defaultdict(list)
  with open(path, "r") as file:
    for row in csv.DictReader(file, delimiter="\t"):
      row["Median"] = ""
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

d1 = read_extracted_files("z_table_csv_3.tsv")
d2 = read_extracted_files("z_table_csv_4.tsv")
d3 = read_extracted_files("z_table_csv_5.tsv")
d4 = read_extracted_files("z_table_csv_6.tsv")
d5 = read_extracted_files("z_table_csv_7.tsv")


df_lis = [d1,d2,d3,d4,d5]
filename = "combine_5_exp.tsv"
combine_files(df_lis, filename)

def add_median(data):
  m_data = defaultdict(list)
  for key, values in data.items():
    track_score = []
    for v in values:
        track_score.append((v.get("Z Score")))
    
    srt_lis = sorted(track_score)
    median = statistics.median(srt_lis)
    for v in values:
        if v.get("Z Score") == median:
            v["Median"] = "X"

        m_data[key].append(v)      
  
  return m_data

u_filename = "Updated.tsv"
df = read_extracted_files(filename)
c_d = add_median(df)
c_d = [c_d]
combine_files(c_d, u_filename)
