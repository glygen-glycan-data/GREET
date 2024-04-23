from app import *


def rm_per_gene(gn_group):
    ml_score = RecallScore()
    precision=0.9
    for tissue in tissues_names:
        print(tissue)
        glyco_enz_set1_data = EnzymeData(extracted_dataset, gn_group)
        exp_violin_plt(glyco_enz_set1_data.get_gen_dataset())
        glyco_enz_set1_data.add_parameters(tissue)
        gnt = glyco_enz_set1_data.get_gen_dataset() 
        re = Report(gnt)
        for i in range(len(gn_group)):
            print(i)
            pr_dic_scores, cdf_scores = re.execute_report(i, rm_exp=True)
            cdf = ml_score.extract_cdf_scores(cdf_scores)
            pr = ml_score.extract_pr_scores(pr_dic_scores)
    
        plots = PRplots(cdf, pr, precision, plt_show=True, plt_save=False)
        #plots.normalized_plt()  



gn_sets = GeneSet()
all_sets = gn_sets.get_sdbox_data()
gn_set1 = gn_sets.extract_glyco_set_at('Fucp', '3', 'a', 'GlcpNAc')
rm_per_gene(gn_set1)







"""

#enz_set format or could be full sets of Enzymes:
#['form_name'], ['site'], ["anomer"], ['parent_form_name']
gen_set = ('Fucp', '3', 'a', 'GlcpNAc')
ext_enz_set = enzymes_dict_set.get(gen_set)


def rm_per_gene(enz_set, violin=False, Normalized=True, Original=False, cdf=False, bx_plt=False):
    rm_exp = Executing_Experiment()
    rm_exp.add_enz_set(enz_set)
    rm_exp.filtered_enz_dataset()
    tissues_names = rm_exp.prep_dataframe()
    if violin == True:
        rm_exp.setup_dataframe() 
        rm_exp.exp_violin_plt()
    
    for t in tissues_names:
        rm_exp.set_pr_cdf_()
        rm_exp.setup_dataframe()
        rm_exp.adding_parameters(t)
        for i in range(len(enz_set)):
            rm_exp.execute_report(i,rm_exp=True)
            rm_exp.extract_pr_scores()
            rm_exp.extract_cdf_scores()

        rm_exp.extract_std()
        rm_exp.extract_zscore()

        if Normalized == True:
            rm_exp.exp_normalized_plt(0.9)
        if Original == True:
            rm_exp.exp_original_plt()
        if cdf == True:
            rm_exp.exp_cdf_plt()
        if bx_plt == True:
            rm_exp.exp_box_plt()
       


#rm_per_gene(ext_enz_set, Original=True)

"""

