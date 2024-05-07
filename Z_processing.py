from plots import * 
from preprocessing_setup import GeneSet

def read_zfile(file, dir="data/"):
    z_data = defaultdict(list)
    file_path = os.path.join(dir, file)
    with open(file_path, "r") as f:
        dict_reader = csv.DictReader(f, delimiter="\t")
        #headers Group Name, Gene, Tissue, Z-Score, Recall-Score
        for row in dict_reader:
            z_data[(row["Group Name"], row["Gene"])].append({row["Tissue"]: (float(row["Z-Score"]), float(row["Recall-Score"]))})

    return z_data


#z_exp_1 = read_zfile("z_table_csv_1")
#z_exp_50_rand = read_zfile("z_table_csv_50")
z_exp_100_rand = read_zfile("z_table_csv_100")




#z_plot(z_exp_1, all_set=True, z_process=True, plt_show=True, plt_save=False)
#z_table(z_exp_1, z_processing=True)
z_plot(z_exp_100_rand, all_set=True, z_process=True, plt_show=True, plt_save=False)
z_table(z_exp_100_rand, z_processing=True)

gn = GeneSet()
print(len(gn.get_sdbox_data()))
print(len(gn.get_all_glyco_enz()))