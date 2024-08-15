
import scanpy as sc 
from app import GeneSet
import pandas as pd  
import numpy as np
from scipy.sparse import issparse, csr_matrix
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from app import Report, RecallScore, Scores, PRplots, create_random_sets
from ml_parameters import *
from plots import *
from sklearn.preprocessing import MaxAbsScaler
import time, multiprocessing
from queue import Empty




url = "https://storage.googleapis.com/adult-gtex/single-cell/v9/snrna-seq-data/GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad"
file = "data/GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad"
adata = sc.read_h5ad(file)

gn_sets = GeneSet()
all_sets = gn_sets.get_sdbox_data()
glyco_enzymes = gn_sets.get_all_glyco_enz()

non_glyco_set = adata.var.loc[~adata.var_names.isin(glyco_enzymes)]
non_glyco_genes = non_glyco_set.index
non_glyco_genes = non_glyco_genes.tolist()


names = ["Logistic Regression"]
classifiers = [LogisticRegression(max_iter=500, solver= "liblinear", penalty="l1", class_weight="balanced")]
ml_names_classi = {}
for i in range(len(names)):
  ml_names_classi[names[i]] = classifiers[i]


def print_process(process_name, start_time):
    end_time = time.time()
    elapsed_time = end_time - start_time
    timestamp = time.strftime("%H:%M:%S", time.localtime(end_time))
    print(f"{timestamp}: {process_name} - Elapsed time: {elapsed_time:.2f} seconds")
    return end_time


def reduce_df(data, enz_set):
    modified_reduced_df = {}
    for index, r in data.iterrows():
        ind = index.split("_")[0]
        tis = index.split("-")[1]
        for en in enz_set:
            unique_identifier = (ind + "." + tis + "." + r["Broad cell type"], en)
            unique_id_data = modified_reduced_df.get(unique_identifier)
            if unique_id_data:    
                unique_id_data += r[en]
                modified_reduced_df[unique_identifier] = unique_id_data
            else:
                modified_reduced_df[unique_identifier] = r[en]
        
    cell_tis_records = defaultdict(list)
    for keys, value in modified_reduced_df.items():
        tissue = keys[0].split(".")[1]
        cell_type = keys[0].split(".")[2]
        
        cell_tis_records[(tissue, cell_type)].append(value)
    
    record_lis = []
    df = pd.DataFrame(columns=["Cell Type", enz_set])
    
    for key, value in cell_tis_records.items():
        for v in value: 
            record_lis.append({"Tissue": key[0], "Broad cell type": key[1],  enz_set[0]: v})

    df = pd.DataFrame(record_lis)

    
    return df


def make_df(enz_set, reduce_func=reduce_df, modify_df=False,  raw_adata=adata):
    cm_gn = {}
    gn_exp_data = {}
    for gene in enz_set:
        if gene in raw_adata.var_names:
           gn_in = np.where(raw_adata.var_names == gene)[0]
           gn_expression = raw_adata.X[:, gn_in]
           gn_exp_data[gene] = gn_expression
           cm_gn[gene] = gn_in 
    
    gn_exp_df = pd.DataFrame(index=raw_adata.obs_names)
    for gn, sp_matx in gn_exp_data.items():
        gn_df = pd.DataFrame.sparse.from_spmatrix(sp_matx, index=raw_adata.obs_names, columns=[gn])
        gn_exp_df = pd.concat([gn_exp_df, gn_df], axis=1)

    obs_df = raw_adata.obs[['Broad cell type', 'Tissue']]
    gn_exp_df = pd.concat([obs_df, gn_exp_df], axis=1)
    

    if modify_df == True:
        gn_exp_df = reduce_func(gn_exp_df, enz_set)

    gn_exp_df = gn_exp_df[gn_exp_df["Broad cell type"] != "Unknown"]
    
    gn_exp_df = gn_exp_df.set_index(["Broad cell type", "Tissue"])
    return gn_exp_df


merged_data = make_df(glyco_enzymes)


Tissue_cells_types = merged_data.index.unique()
combined_tis_ct = set()
for ct in Tissue_cells_types:
    tis = ct[1]
    cell_type = ct[0]
    combined_tis_ct.add(ct)
    combined_tis_ct.add(tis)
    combined_tis_ct.add(cell_type)



gn_lis = []
gn_not_indata = []
for gn_set in all_sets:
  for gn in gn_set:
    temp_df = merged_data.get(gn)
    if temp_df is not None:
      gn_lis.append(gn)
    else:
      gn_not_indata.append(gn)



def find_ct_less_100(m_df, unique_ind):
    ct_less_100 = []
    updated_cell_types = []
    ct_types_counts = {}
    single_ct_track = []
    for ct in unique_ind:
        if type(ct) == tuple:
            m_df["Class"] = (m_df.index == ct).astype(int)
            
        else:
            index_flat = m_df.index.map(lambda x: " ".join(map(str,x)))
            m_df["Class"] = index_flat.str.contains(ct, regex=False)

        count_class_1 = (m_df["Class"] == 1).sum()
        ct_types_counts[ct] = count_class_1
        
        track_id = 0
        for check_type in unique_ind:
            if type(ct) != tuple:
                if ct == check_type or ct in check_type:
                    track_id += 1
    
        if track_id == 2:
            ct_less_100.append(ct)
            single_ct_track.append(ct)

        elif count_class_1 < 200:
            ct_less_100.append(ct)
        elif ct == "Unknown":
            ct_less_100.append(ct)
    
        else:
            updated_cell_types.append(ct)
   
    return ct_less_100, updated_cell_types, ct_types_counts, single_ct_track


drop_cts, reduced_cell_types, ind_counts, single_cts_track = find_ct_less_100(merged_data, combined_tis_ct)




def ct_name(cell_type_name, single_cts_track=single_cts_track):
    if type(cell_type_name) == tuple:
        if cell_type_name[0] in single_cts_track:
            n_ct = f"{cell_type_name[0]} in [{cell_type_name[1]}]"
            print(n_ct)
        else:
            n_ct = f"{cell_type_name[0]} in {cell_type_name[1]}"

    else:
        n_ct = cell_type_name

    return n_ct



print(len(reduced_cell_types))

print(len(reduced_cell_types))

### Run this script once to extract Percent Detected and Index Counts per CT type ####

"""

def percent_detected(df, unique_ind, gly_enz, gly_not_indata):
    p_d = defaultdict(list)
    for enz in gly_enz:
        if enz in gly_not_indata:
           continue

        print(enz)
        over_all_start = time.time()
        for ct in unique_ind:
            if type(ct) == tuple:
                df["Class"] = (df.index == ct).astype(int)
            
            else:
                index_flat = df.index.map(lambda x: " ".join(map(str,x)))
                df["Class"] = index_flat.str.contains(ct, regex=False)


            name_ct = ct_name(ct)
            for i in range(2):
                filtered_rows = df[df["Class"] == i]
                non_zero = len(filtered_rows[filtered_rows[enz] != 0])
                total_cells = len(filtered_rows[enz])
                detected = (non_zero/total_cells) * 100
                p_d[(enz, name_ct)].append(detected)

        print_process("Whole Set", over_all_start)
    
    return p_d




per_det = percent_detected(merged_data, combined_tis_ct, glyco_enzymes, gn_not_indata)

full_path = os.path.join("data/" + "percent_detected_count.tsv")
with open(full_path, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter="\t")
    head = ["Gene", "Tissue", "Percent Detected Class 0", "Percent Detected Class 1"]
    writer.writerow(head)
    for k,v in per_det.items():
        writer.writerow([k[0], k[1], v[0], v[1]]



full_path = os.path.join("data/" + "ind_count.tsv")
with open(full_path, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter="\t")
    head = ["Old CT", "New CT", "Count"]
    writer.writerow(head)
    for k,v in ind_counts.items():
        if k in reduced_cell_types:
            n_k = ct_name(k)
            writer.writerow([k, n_k, v])

"""

#print(gn_not_indata)


def single_enz_experiment(gene_set, all_enzymes, non_glyco_genes, single_cts_track):     
    over_all_start = time.time()
    cr = create_random_sets(gene_set, non_glyco_genes, all_enzymes, random_size=20)
    precision = 0.9
    col_zscore = {}
    for ct in reduced_cell_types:
        ml_score = RecallScore()
        for rand_num, ext_enz_set in cr.items():
            gnt = make_df(ext_enz_set)
            if type(ct) == tuple:
                gnt["Class"] = (gnt.index == ct).astype(int)

            
            else:
                index_flat = gnt.index.map(lambda x: " ".join(map(str,x)))
                gnt["Class"] = index_flat.str.contains(ct, regex=False)


            re = Report(gnt)
            pr_dic_scores, cdf_scores = re.execute_report(rand_num, ml_names_classi)
            cdf = ml_score.extract_cdf_scores(cdf_scores)
            pr = ml_score.extract_pr_scores(pr_dic_scores)



        name_ct = ct_name(ct)
        if cdf:
            sc = Scores(cdf)    
            col_zscore[name_ct] = sc.extract_zscore()
        

            #plots = PRplots(cdf, pr, precision, plt_show=True, plt_save=False)
            #plots.normalized_plt()
            #plots.box_plt()
            #plots.cdf_plt()
            #plots.histo_plt()
        

    #z_plot(col_zscore, plt_show=False, plt_save=False, title=gene_set)
    print_process("Whole Set", over_all_start)
    return col_zscore





class Workers:
    def __init__(self, glyco_enz, non_genes, data, genes_nt_data):
        self.glyco_enz = sorted(glyco_enz)
        self.non_gly_gn = non_genes
        self.data = data
        self.gn_nt_data = genes_nt_data

    def worker(self, queue, worker_id):
        total_gr_zscore = defaultdict(list)
        while True:
            try:
                gn = queue.get(timeout=2) 
            except Empty:
                break  

            if gn not in self.gn_nt_data:
                print(f"Worker {worker_id}", gn)
                test_gn = [gn]
                total_gr_zscore[(f"Worker {worker_id}", gn)] = single_enz_experiment(test_gn, self.glyco_enz, self.non_gly_gn, self.data)
                save_zdata(total_gr_zscore, f"Worker {worker_id}")

def worker_function(worker_id, glyco_enzymes, non_glyco_genes, merged_data, gns_nt_data, queue):
    w = Workers(glyco_enzymes, non_glyco_genes, merged_data, gns_nt_data)
    w.worker(queue, worker_id)

if __name__ == "__main__":
    glyco_enzymes = glyco_enzymes
    non_glyco_genes = non_glyco_genes
    merged_data = merged_data
    gn_not_indata = gn_not_indata


    queue = multiprocessing.Queue()
    for gn in glyco_enzymes:
        queue.put(gn)


    jobs = []
    for i in range(1, 5):
        p = multiprocessing.Process(target=worker_function, args=(i, glyco_enzymes, non_glyco_genes, merged_data, gn_not_indata, queue))
        jobs.append(p)
        p.start()


    for job in jobs:
        job.join()