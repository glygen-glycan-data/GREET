from app import *
import time


#enz_set format or could be full sets of Enzymes:
#['form_name'], ['site'], ["anomer"], ['parent_form_name']

def print_process(process_name, start_time):
    end_time = time.time()
    elapsed_time = end_time - start_time
    timestamp = time.strftime("%H:%M:%S", time.localtime(end_time))
    print(f"{timestamp}: {process_name} - Elapsed time: {elapsed_time:.2f} seconds")
    return end_time



test = extracted_dataset.get("GGTA1")
#for t, d in test.items():
#    print(t, len(d))

def set_enz_experiment(gene_set, all_enzymes, dataset):     
    over_all_start = time.time()
    cr = create_random_sets(gene_set, dataset, all_enzymes)
    precision = 0.9
    col_zscore = {}
    for tissue in tissues_names:
        start_time = time.time()
        ml_score = RecallScore()
        for rand_num, ext_enz_set in cr.items():
            glyco_enz_set_data = EnzymeData(dataset, ext_enz_set)
            glyco_enz_set_data.add_parameters(tissue)
            gnt = glyco_enz_set_data.get_gen_dataset() 
            glyco_enz_set_data.reset() 
            re = Report(gnt)
            pr_dic_scores, cdf_scores = re.execute_report(rand_num, ml_names_classi)
            cdf = ml_score.extract_cdf_scores(cdf_scores)
            pr = ml_score.extract_pr_scores(pr_dic_scores)
            
        sc = Scores(cdf)
        col_zscore[tissue] = sc.extract_zscore()


        #plots = PRplots(cdf, pr, precision, plt_show=False, plt_save=False)
        #plots.normalized_plt()
        #plots.box_plt()
        #plots.cdf_plt()
        #plots.histo_plt()
        #print_process(tissue, start_time)
        
    
    #z_plot(col_zscore, plt_show=False, plt_save=False, title=gene_set)
    print_process("Whole Set", over_all_start)
    return col_zscore


#set_enz_experiment("GGTA1", extracted_dataset)



gn_sets = GeneSet()
#gn_set1 = gn_sets.extract_glyco_set_at('Fucp', '3', 'a', 'GlcpNAc')  
all_sets = gn_sets.get_sdbox_data()
glyco_enzymes = gn_sets.get_all_glyco_enz() 
total_gr_zscore = defaultdict(list)


for i, (gn_group, group_name) in enumerate(all_sets.items()):
    gn_group_set = ""
    for gn in gn_group:
        gn_group_set += f"{gn} "  

    if list(gn_group)[0] != "GGTA1":
        print(gn_group_set)
        total_gr_zscore[(group_name, gn_group_set)] = set_enz_experiment(gn_group, glyco_enzymes, extracted_dataset)



save_zdata(total_gr_zscore, 4)







