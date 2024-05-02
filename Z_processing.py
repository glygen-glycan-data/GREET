from plots import * 


def read_zfile(file, dir="data/"):
    z_data = defaultdict(list)
    file_path = os.path.join(dir, file)
    with open(file_path, "r") as f:
        dict_reader = csv.DictReader(f, delimiter="\t")
        #headers Group Name, Gene, Tissue, Z-Score, Recall-Score
        for row in dict_reader:
            z_data[(row["Group Name"], row["Gene"])].append({row["Tissue"]: (float(row["Z-Score"]), float(row["Recall-Score"]))})

    return z_data


z_exp_1 = read_zfile("z_table_csv_1")
z_exp_50_rand = read_zfile("z_table_csv_2")


z_plot(z_exp_1, all_set=True, z_process=True, plt_show=True, plt_save=False)
#z_table(z_exp_1, z_processing=True)
z_plot(z_exp_50_rand, all_set=True, z_process=True, plt_show=True, plt_save=False)
z_table(z_exp_50_rand, z_processing=True)