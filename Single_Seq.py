
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


url = "https://storage.googleapis.com/adult-gtex/single-cell/v9/snrna-seq-data/GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad"
file = "data/GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad"
adata = sc.read_h5ad(file)

gn_sets = GeneSet()
all_sets = gn_sets.get_sdbox_data()
glyco_enzymes = gn_sets.get_all_glyco_enz()

non_glyco_set = adata.var.loc[~adata.var_names.isin(glyco_enzymes)]
non_glyco_genes = non_glyco_set.index
non_glyco_genes = non_glyco_genes.tolist()


##Filter and Normalize ####
#sc.pp.filter_cells(adata, min_genes=200)
#sc.pp.filter_genes(adata, min_cells=3)
#upper_lim = np.quantile(adata.obs.n_genes.values, 0.98)
#lower_lim = np.quantile(adata.obs.n_genes.values, 0.02)
#print(f'{lower_lim} to {upper_lim}')
#adata = adata[(adata.obs.n_genes < upper_lim) & (adata.obs.n_genes > lower_lim)]
#sc.pp.normalize_total(adata, target_sum=1e4)
#sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

"""
"Scikit-Learn"
scaler = MaxAbsScaler()
adata.X = scaler.fit_transform(adata.X)
adata.X = csr_matrix(adata.X)
"""

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

    gn_exp_df[gn_exp_df["Broad cell type"].str.contains("Unknown")==False]
    
    gn_exp_df = gn_exp_df.set_index(["Broad cell type", "Tissue"])
    return gn_exp_df



def find_ct_less_100(m_df):
    cell_types = adata.obs["Broad cell type"].unique()
    ct_less_100 = []
    updated_cell_types = []
    for ct in cell_types:
        m_df["Class"] = (m_df.index == ct).astype(int)
        count_class_1 = (m_df["Class"] == 1).sum()
        if count_class_1 < 100:
            ct_less_100.append(ct)
        elif ct == "Unknown":
            ct_less_100.append(ct)
        else:
            updated_cell_types.append(ct)
   
    return ct_less_100, updated_cell_types


merged_data = make_df(glyco_enzymes)
#drop_cts, cell_types = find_ct_less_100(merged_data)
#print(drop_cts)
Tissue_cells_types = merged_data.index.unique()

gn_lis = []
gn_not_indata = []
for gn_set in all_sets:
  for gn in gn_set:
    temp_df = merged_data.get(gn)
    if temp_df is not None:
      gn_lis.append(gn)
    else:
      gn_not_indata.append(gn)


names = ["Logistic Regression"]
classifiers = [LogisticRegression(max_iter=500, solver= "liblinear", penalty="l1", class_weight="balanced")]
ml_names_classi = {}
for i in range(len(names)):
  ml_names_classi[names[i]] = classifiers[i]


def single_enz_experiment(gene_set, all_enzymes, non_glyco_genes, dataset, make_df):     
    cr = create_random_sets(gene_set, non_glyco_genes, all_enzymes, random_size=20)
    precision = 0.9
    col_zscore = {}
    #for ct in cell_types:
    cell_type_test = ('Endothelial cell (vascular)', 'Heart')
    for ct in Tissue_cells_types:
        #if ct != cell_type_test:
         # continue
        print(ct)
        ml_score = RecallScore()
        for rand_num, ext_enz_set in cr.items():
            gnt = make_df(ext_enz_set)
            #gnt = gnt[~gnt.index.isin(drop_cts)]
            gnt["Class"] = (gnt.index == ct).astype(int)
            filtered_rows = gnt[gnt["Class"] == 1]
            #for enz in ext_enz_set:
            #if any(len(filtered_rows[filtered_rows[enz] != 0]) > 100  for enz in ext_enz_set):
            re = Report(gnt)
            pr_dic_scores, cdf_scores = re.execute_report(rand_num, ml_names_classi)
            cdf = ml_score.extract_cdf_scores(cdf_scores)
            pr = ml_score.extract_pr_scores(pr_dic_scores)
            #else:
            #     cdf = []
        if cdf:
            sc = Scores(cdf)
            #print(ct, sc.extract_zscore())
            col_zscore[ct] = sc.extract_zscore()
        

            #plots = PRplots(cdf, pr, precision, plt_show=True, plt_save=False)
            #plots.normalized_plt()
            #plots.box_plt()
            #plots.cdf_plt()

            #plots.histo_plt()
        
    
    #z_plot(col_zscore, plt_show=False, plt_save=False, title=gene_set)
    return col_zscore


def run(enz_sets, glyco_enz, extracted_data, gene_nt_df, filename):
    total_gr_zscore = defaultdict(list)
    for i, (gn_group, group_name) in enumerate(enz_sets.items()):
        gn_group_set = ""
        gn_set = []
        for gn in gn_group:
            gn_group_set += f"{gn} " 
            if gn not in gene_nt_df:
                gn_set.append(gn)
        
        if list(gn_group)[0] != "GGTA1":
            #print(gn_set)
            print(gn_group_set)
            
            total_gr_zscore[(group_name, gn_group_set)] = single_enz_experiment(gn_set, glyco_enzymes, non_glyco_genes, merged_data, make_df)
        


    save_zdata(total_gr_zscore, filename)



#r = run(all_sets, glyco_enzymes, merged_data,gn_not_indata, "test")


test_set1 = ["ST6GALNAC1", "ST6GALNAC2", "ST6GALNAC3", "ST6GALNAC4", "ST6GALNAC5", "ST6GALNAC6"]
test_set = "ST6GALNAC1, ST6GALNAC2, ST6GALNAC3, ST6GALNAC4, ST6GALNAC5, ST6GALNAC6"

#test_set = ["CHST8", "CHST9"]
total_gr_zscore = defaultdict(list)
total_gr_zscore[("ST6", test_set)] = single_enz_experiment(test_set1, glyco_enzymes, non_glyco_genes, merged_data, make_df)
save_zdata(total_gr_zscore, "ST6GALNAC3")   