from collections import defaultdict
import csv
import app, os

files = ["z_table_csv_Worker 1", "z_table_csv_Worker 2",
         "z_table_csv_Worker 3", "z_table_csv_Worker 4"]



ind_counts = {}
for row in csv.DictReader(open("data/ind_count.tsv"),delimiter="\t"):
    ind_counts[row["Gene"]] = row["Count"]


percent_detected = {}
for row in csv.DictReader(open("data/percent_detected_count.tsv"),delimiter="\t"):
    percent_detected[(row["Gene"], row["Tissue"])] = row["Percent Detected"]



cb = defaultdict()
gn_count = set()
for f in files:
    for row in csv.DictReader((open(f"data/{f}.tsv")), delimiter="\t"):
        row["Count"] = ind_counts.get(row["Tissue"])
        key = (row["Gene"], row["Tissue"])
        row["Percent Detected"] = percent_detected.get(key)
        cb[key] = row
        gn_count.add(row["Gene"])

        
print(len(gn_count))

filename = "Single_seq_comb_2.tsv"
full_path = os.path.join("data/" + filename)
with open(full_path, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter="\t")
    head = ["Group Name", "Gene", "Tissue", "Z Score", "Real Set Recall","Mean",  "SD", "Old SD"] + [f"Set {i + 1}" for i in range(20)]  + ["Counts", "Percent Detected"]
    writer.writerow(head)
    for values in cb.values():
        writer.writerow([v for v in values.values()])

