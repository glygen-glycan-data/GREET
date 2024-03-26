from app import *


def enz_set_tissue(enz_set, violin=False, Normalized=True, Original=False, cdf=False, bx_plt=False):
    exp_enz_tissue = Executing_Experiment()
    exp_enz_tissue.add_enz_set(enz_set)
    exp_enz_tissue.filtered_enz_dataset()
    tissues_names = exp_enz_tissue.prep_dataframe()
    if violin == True:
        exp_enz_tissue.setup_dataframe()
        exp_enz_tissue.exp_violin_plt()

    exp_enz_tissue.set_pr_cdf_()
    for t in tissues_names:
        exp_enz_tissue.setup_dataframe()
        exp_enz_tissue.adding_parameters(t)
        exp_enz_tissue.execute_report(t)
        exp_enz_tissue.extract_cdf_scores()
        exp_enz_tissue.extract_pr_scores()
        
    
    if Normalized == True:
        exp_enz_tissue.exp_normalized_plt(0.9)
    if Original == True:
        exp_enz_tissue.exp_original_plt()
    if cdf == True:
        exp_enz_tissue.exp_cdf_plt()
    if bx_plt == True:
        exp_enz_tissue.exp_box_plt()



#enz_set could be set of Enzymes:
#enz_set format ::: ['form_name'], ['site'], ["anomer"], ['parent_form_name']
        
gen_set = ('Fucp', '3', 'a', 'GlcpNAc')
ext_enz_set = enzymes_dict_set.get(gen_set)

enz_set_tissue(ext_enz_set)
       
