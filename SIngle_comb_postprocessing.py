from collections import defaultdict
import csv
import app, os, sys


files = ["z_table_csv_Worker 11", "z_table_csv_Worker 12",
        "z_table_csv_Worker 13", "z_table_csv_Worker 14"]


#files = ["z_table_csv_TEST_SIN"]

ind_counts = {}
#ct_convert = {}
for row in csv.DictReader(open("data/ind_count.tsv"),delimiter="\t"):
    ind_counts[row["New CT"]] = row["Count"]
    #ct_convert[row["Old CT"]] = row["New CT"]

percent_detected = {}
for row in csv.DictReader(open("data/percent_detected_count.tsv"),delimiter="\t"):
    percent_detected[(row["Gene"], row["Tissue"])] = [row["Percent Detected Class 0"], row["Percent Detected Class 1"]]



cb = defaultdict()
gn_count = set()
for f in files:

    for row in csv.DictReader((open(f"data/{f}.tsv")), delimiter="\t"):
        row["Count"] = ind_counts.get(row["Tissue"])

        key = (row["Gene"], row["Tissue"])
        p_d = percent_detected.get(key)
        
        row["Percent Detected Class 0"] = p_d[0]
        row["Percent Detected Class 1"] = p_d[1]
        
        cb[key] = row
        gn_count.add(row["Gene"])



filename = "T2.tsv"
full_path = os.path.join("data/" + filename)
with open(full_path, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter="\t")
    head = ["Group Name", "Gene", "Tissue", "Z Score", "Real Set Recall","Mean",  "SD", "Old SD"] + [f"Set {i + 1}" for i in range(20)]  + ["Counts", "Percent Detected Class 0", "Percent Detected Class 1"]
    writer.writerow(head)
    for values in cb.values():
        writer.writerow([v for v in values.values()])

