from app import *
from ml_parameters import *
from plots import *
from sklearn.preprocessing import MaxAbsScaler
from collections import defaultdict
from queue import Empty
from scipy.sparse import issparse, csr_matrix
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
import os
import time, multiprocessing
import scanpy as sc 
import pandas as pd  
import numpy as np



#url = "https://storage.googleapis.com/adult-gtex/single-cell/v9/snrna-seq-data/GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad"
#file = "data/GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad"


adata = sc.read_h5ad(datafile)
# print(adata)

gn_sets = GeneSet()
all_sets = gn_sets.get_sdbox_data()
glyco_enzymes = gn_sets.get_all_glyco_enz()
# print(glyco_enzymes)
# print(all_sets)
    
non_glyco_set = adata.var.loc[~adata.var_names.isin(glyco_enzymes)]
non_glyco_genes = non_glyco_set.index
non_glyco_genes = non_glyco_genes.tolist()



names = ["Logistic Regression"]
classifiers = [LogisticRegression(max_iter=500, solver= "liblinear", penalty="l1", class_weight="balanced")]
ml_names_classi = {}
for i in range(len(names)):
  ml_names_classi[names[i]] = classifiers[i]



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
    # print(enz_set)
    cm_gn = {}
    gn_exp_data = {}
    for gene in enz_set:
        # print(gene,raw_adata.var_names)
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
# print(merged_data.info())
# print(merged_data)


Tissue_cells_types = merged_data.index.unique()
combined_tis_ct = set()
for ct in Tissue_cells_types:
    tis = ("",ct[1])
    cell_type = (ct[0],"")
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



def find_ct_less_100(m_df, unique_ind, ctc=threshold_count):
    ct_less_100 = []
    updated_cell_types = []
    ct_types_counts = {}
    single_ct_track = []
    for ct in sorted(unique_ind):
        cts = ":" + ((":".join(map(str,ct))).strip(":")) + ":"
        index_flat = m_df.index.map(lambda x: ":"+(":".join(map(str,x)))+":")
        m_df["Class"] = index_flat.str.contains(cts, regex=False)

        count_class_1 = (m_df["Class"] == 1).sum()
        ct_types_counts[ct] = count_class_1
        # print(ct,count_class_1)

        if count_class_1 < ctc: #check config file
            ct_less_100.append(ct)
        elif ct == "Unknown":
            ct_less_100.append(ct)
        else:
            updated_cell_types.append(ct)
    separate_counts = defaultdict(int)
    for ct in sorted(updated_cell_types):
        if ct[0] and ct[1]:
            separate_counts[ct[0]] += 1
            separate_counts[ct[1]] += 1
    for ct in sorted(list(updated_cell_types)):
        if not ct[0]:
            if separate_counts[ct[1]] <= 1:
                updated_cell_types.remove(ct)
                ct_less_100.append(ct)
                single_ct_track.append(ct[1])
        elif not ct[1]:
            if separate_counts[ct[0]] <= 1:
                updated_cell_types.remove(ct)
                ct_less_100.append(ct)
                single_ct_track.append(ct[0])

    return ct_less_100, sorted(updated_cell_types), ct_types_counts, single_ct_track


# print("combined_tis_ct",sorted(combined_tis_ct))
drop_cts, reduced_cell_types, ind_counts, single_cts_track = find_ct_less_100(merged_data, combined_tis_ct)
# print("drop_cts",sorted(drop_cts))
# print("reduced_cell_types",sorted(reduced_cell_types))
# print("single_cts_track",single_cts_track)

def ct_name(cell_type_name, single_cts_track=single_cts_track):
    if cell_type_name[0] in single_cts_track:
        n_ct = f"{cell_type_name[0]} [{cell_type_name[1]}]"
    elif cell_type_name[1] in single_cts_track:
        n_ct = f"{cell_type_name[0]} [{cell_type_name[1]}]"
    elif not cell_type_name[0]:
        n_ct = "Tissue: "+cell_type_name[1]
    elif not cell_type_name[1]:
        n_ct = "Cell-type: "+cell_type_name[0]
    else:
        n_ct = f"{cell_type_name[0]} in {cell_type_name[1]}"
    return n_ct


### Run this script once to extract Percent Detected and Index Counts per CT type ####

"""

def percent_detected(df, unique_ind, gly_enz, gly_not_indata):
    p_d = defaultdict(list)
    for enz in gly_enz:
        if enz in gly_not_indata:
           continue

        # print(enz)
        over_all_start = time.time()
        for ct in unique_ind:
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
# print(per_det)

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

def experiment(gene_set, all_enzymes, dataset, index_set, experiment_type = exp_type,  rts = random_test_size, rpa=recall_precision_threshold):     
    over_all_start = time.time()
    cr = create_random_sets(gene_set[2], dataset, all_enzymes, random_size=rts) #rts, check cofiguration file
    # print(cr)
    precision = rpa #rpa, check cofiguration file
    col_zscore = {}

    # index_set = index_set[:20]

    # i.e index are tissue_names in tissue experiment, and cell_types in cell experiemnt
    for j,ind in enumerate(index_set): 
        ml_score = RecallScore()
    
        for rand_num, ext_enz_set in cr.items():
            if experiment_type == "Tissue":
                glyco_enz_set_data = EnzymeData(dataset, ext_enz_set)        
                glyco_enz_set_data.add_parameters(ind)
                gnt = glyco_enz_set_data.get_gen_dataset() 
                glyco_enz_set_data.reset() 

            elif experiment_type == "Cell":
                gnt = make_df(ext_enz_set)
                
                inds = ":" + ((":".join(map(str,ind))).strip(":")) + ":"
                index_flat = gnt.index.map(lambda x: ":"+(":".join(map(str,x)))+":")
                gnt["Class"] = index_flat.str.contains(inds, regex=False)
                    
                re_ind = ct_name(ind) #renaming indices 

            re = Report(gnt)
            pr_dic_scores, cdf_scores = re.execute_report(rand_num, ml_names_classi)
            cdf = ml_score.extract_cdf_scores(cdf_scores)
            pr = ml_score.extract_pr_scores(pr_dic_scores)

        if cdf:
            sc = Scores(cdf)
            if experiment_type == "Cell":
                col_zscore[re_ind] = sc.extract_zscore()
            else:
                col_zscore[ind] = sc.extract_zscore()

        #plots = PRplots(cdf, pr, precision, plt_show=False, plt_save=False)
        #plots.normalized_plt()
        #plots.box_plt()
        #plots.cdf_plt()
        #plots.histo_plt()
        print("%d/%d"%(j+1,len(index_set)), "%d/%d"%(gene_set[0]+1,gene_set[1]), (" ".join(ind)).strip(), " ".join(sorted(gene_set[2])))
        sys.stdout.flush()
        
    
    #z_plot(col_zscore, plt_show=False, plt_save=False, title=gene_set)

    print_process("Whole Set", over_all_start)

    return col_zscore



def make_job_indices(total_workers,job_number, glyco_en):
    n_job_ind = []
    glyco_en = sorted(glyco_en)
    assert job_number <= total_workers, "Job number exceed the total number of workers"
    start_index = job_number - 1
    ln_glyco = len(glyco_en) - 1
    while start_index < ln_glyco or start_index <= ln_glyco:
        n_job_ind.append(start_index)
        start_index += total_workers

    return n_job_ind



class Workers:
    def __init__(self, glyco_enz, non_genes, data, genes_nt_data, ct_indices):
        self.glyco_enz = sorted(glyco_enz)
        self.non_gly_gn = non_genes
        self.data = data
        self.gn_nt_data = set(genes_nt_data)
        self.ct_names = ct_indices

    def worker(self, queue, worker_id):
        total_gr_zscore = defaultdict(list)
        while True:
            try:
                task = queue.get(timeout=2) 
            except Empty:
                break  

            if len(task[-1].intersection(self.gn_nt_data)) == 0:
                # print(f"Worker {worker_id}", gn)
                # test_gn = gn
                #total_gr_zscore[(f"Worker {worker_id}", gn)] = single_enz_experiment(test_gn, self.glyco_enz, self.non_gly_gn, self.data)
                total_gr_zscore[(f"Worker {worker_id}", ", ".join(sorted(task[-1])))] = experiment(task, self.glyco_enz, self.non_gly_gn, index_set = self.ct_names)

        return total_gr_zscore
    
def worker_function(worker_id, glyco_enzymes, non_glyco_genes, merged_data, gns_nt_data, ct_ind, queue, result_dict):
    w = Workers(glyco_enzymes, non_glyco_genes, merged_data, gns_nt_data, ct_ind)
    result_dict[worker_id] = w.worker(queue, worker_id)




if __name__ == "__main__":

    cpuCount = os.cpu_count()
    print("Number of CPUs in the system:", cpuCount, "\n")

    if not (args.workers and args.nworker):
        # parser.print_help()  # Print help message
        args.workers = 1
        args.nworker = 1

    if args.processors:
        p = args.processors
    else:
        p = cpuCount
        
    w = args.workers 
    j = args.nworker

    print(f"Number of Workers: {args.workers}")
    print(f"Job Number: {args.nworker}")
    print(f"Number of Processors: {p}")

    # glyco_enzyme_sets = [ (e,) for e in sorted(glyco_enzymes) ]
    # glyco_enzyme_sets = [ set((e1,e2)) for e1 in sorted(glyco_enzymes) for e2 in sorted(glyco_enzymes) if e1 < e2 ]

    glyco_enzyme_sets = set()
    allmonos = gn_sets.names()
    for m in sorted(allmonos):
        print(m)
    # allmonos = ['NeupNAc-a6-GalpNAc','GalpNAc-b4-GlcpNAc']
    allmonopairs = [('NeupNAc-a6-Galp','Galp-b4-GlcpNAc'),('NeupNAc-a6-GalpNAc','GalpNAc-b4-GlcpNAc'),
                    ('sulfate-n4-GalpNAc','GalpNAc-b4-GlcpNAc')]
    for pair in allmonopairs:
        for n1 in pair:
            for n2 in pair:
                if n1 >= n2:
                    continue
                for e1 in gn_sets.extract_glyco_set_at(n1):
                    for e2 in gn_sets.extract_glyco_set_at(n2):
                        if e1 == e2:
                            continue
                        glyco_enzyme_sets.add(tuple(sorted([e1,e2])))
    glyco_enzyme_sets = list(set(x) for x in glyco_enzyme_sets)
    # glyco_enzyme_sets = glyco_enzyme_sets[:8]

    print(f"Number of enzyme sets: {len(glyco_enzyme_sets)}")
    print(f"Number of cell-types: {len(reduced_cell_types)}")
    job_indices = make_job_indices(w, j, glyco_enzyme_sets)
    print(f"Number of enzymes sets to analyze: {len(job_indices)}")

    queue = multiprocessing.Queue()
    manager = multiprocessing.Manager()
    result_dict = manager.dict()

    for indx in job_indices:
        queue.put((indx,len(glyco_enzyme_sets),glyco_enzyme_sets[indx]))


    jobs = []
    for i in range(1, p + 1):
        p = multiprocessing.Process(target=worker_function, args=(i, glyco_enzymes, non_glyco_genes, merged_data, gn_not_indata, reduced_cell_types, queue, result_dict))
        jobs.append(p)
        p.start()


    for job in jobs:
        job.join()

    combined_results = {}
    for result in result_dict.values():
        for key, value in result.items():
            combined_results[key] = value 

    save_zdata(combined_results, f"Combined_{args.filename}")


