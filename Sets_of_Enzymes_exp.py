from app import *

#enz_set format or could be full sets of Enzymes:
#['form_name'], ['site'], ["anomer"], ['parent_form_name']


gen_set = ('Fucp', '3', 'a', 'GlcpNAc')
p_enz_set = enzymes_dict_set.get(gen_set)
cr = creating_random_sets(gen_set, enzymes_dict_set, non_genes)

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


enz_set_experiment(p_enz_set)