from app import *


def rm_per_gene(gn_group):
    ml_score = RecallScore()
    precision=0.9
    for tissue in tissues_names:
        if tissue != "Cells - Cultured fibroblasts":
            continue
        print(tissue)
        glyco_enz_set1_data = EnzymeData(extracted_dataset, gn_group)
        exp_violin_plt(glyco_enz_set1_data.get_gen_dataset())
        glyco_enz_set1_data.add_parameters(tissue)
        gnt = glyco_enz_set1_data.get_gen_dataset() 
        re = Report(gnt)
        for i in range(len(gn_group)):
            i += 1
            if i != len(gn_group):
                pr_dic_scores, cdf_scores = re.execute_report(i, ml_names_classi, rm_exp=True)
                cdf = ml_score.extract_cdf_scores(cdf_scores)
                pr = ml_score.extract_pr_scores(pr_dic_scores)
    
        plots = PRplots(cdf, pr, precision, plt_show=True, plt_save=False)
        plots.normalized_plt()  



gn_sets = GeneSet()
all_sets = gn_sets.get_sdbox_data()
gn_set1 = gn_sets.extract_glyco_set_at("GlcpNAc-b3-Fucx")
rm_per_gene(gn_set1)



