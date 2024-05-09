from app import *

genes = set(sys.argv[1:])

precision = 0.9
gn_set = GeneSet()
gn_set1 = gn_set.extract_glyco_set_at('Fucp', '3', 'a', 'GlcpNAc')
ml_score = RecallScore()

for tissue in tissues_names:
    glyco_enz_set1_data = EnzymeData(extracted_dataset, gn_set1)
    glyco_enz_set1_data.add_parameters(tissue)
    gnt = glyco_enz_set1_data.get_gen_dataset() 
    re = Report(gnt)
    pr_dic_scores, cdf_scores = re.execute_report(tissue, ml_names_classi)
    cdf = ml_score.extract_cdf_scores(cdf_scores)
    pr = ml_score.extract_pr_scores(pr_dic_scores)


plots = PRplots(cdf, pr, precision, plt_show=True, plt_save=False)
plots.normalized_plt()
#plots.box_plt()
#plots.cdf_plt()
#plots.original_plt()



