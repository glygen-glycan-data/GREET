from app import *

#enz_set format or could be full sets of Enzymes:
#['form_name'], ['site'], ["anomer"], ['parent_form_name']





"""
def enz_set_experiment(enz_set, violin=False, Normalized=True, Original=False, cdf=False, bx_plt=False):
    exp = Executing_Experiment()
    exp.add_enz_set(enz_set)
    exp.filtered_enz_dataset()
    tissues_names = exp.prep_dataframe()
    for t in tissues_names:
        enz_exp = Executing_Experiment()
        enz_exp.set_pr_cdf_()
        for rand_num, ext_enz_set in cr.items():
            enz_exp.add_enz_set(ext_enz_set)
            enz_exp.filtered_enz_dataset()
            enz_exp.prep_dataframe()
            enz_exp.setup_dataframe()

            if violin == True:
                enz_exp.exp_violin_plt()

            enz_exp.adding_parameters(t)
            enz_exp.execute_report(rand_num)
            enz_exp.extract_pr_scores()
            enz_exp.extract_cdf_scores()

            enz_exp.reset()

        if Normalized == True:
            enz_exp.exp_normalized_plt(0.9)
        if Original == True:
            enz_exp.exp_original_plt()
        if cdf == True:
            enz_exp.exp_cdf_plt()
        if bx_plt == True:
            enz_exp.exp_box_plt()

"""





def set_enz_experiment(gene_set):
    gen_set = GeneSet()
    all_enzymes = gen_set.get_all_glyco_enz()
    cr = create_random_sets(gene_set, extracted_dataset, all_enzymes)
    precision = 0.9
    col_zscore = {}
    for tissue in tissues_names:
        ml_score = RecallScore()
        for rand_num, ext_enz_set in cr.items():
            glyco_enz_set_data = EnzymeData(extracted_dataset, ext_enz_set)
            glyco_enz_set_data.add_parameters(tissue)
            gnt = glyco_enz_set_data.get_gen_dataset() 
            glyco_enz_set_data.reset()
            pr_dic_scores, cdf_scores = execute_report(gnt, rand_num)
            cdf = ml_score.extract_cdf_scores(cdf_scores)
            pr = ml_score.extract_pr_scores(pr_dic_scores)
        
        zsc = Z_Score(cdf)
        col_zscore[tissue] = zsc.extract_zscore()
    

        plots = PRplots(cdf, pr, precision, plt_show=False, plt_save=False)
        #plots.normalized_plt()
        #plots.box_plt()
        #plots.cdf_plt()
        #plots.histo_plt()

    z_plot(col_zscore, plt_show=True, plt_save=True, title=gene_set)

    return col_zscore




gn_sets = GeneSet()
gn_set1 = gn_sets.extract_glyco_set_at('Fucp', '3', 'a', 'GlcpNAc')  
all_sets = gn_sets.get_sdbox_data()


total_gr_zscore = defaultdict(list)
for gn_group in all_sets.values():
    print(gn_group)
    total_gr_zscore[(i, str(gn_group))] = set_enz_experiment(gn_group)
    



z_plot(total_gr_zscore, all_set= True, plt_show=False, plt_save=True)
