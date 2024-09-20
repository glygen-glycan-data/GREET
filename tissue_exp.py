from app import *
from queue import Empty
import time, multiprocessing
import os



HPATissue_file = "HPA_tissue.txt"
Gspec_file = "proteinatlas.tsv"
h = HPATissueClassifcation(HPATissue_file, Gspec_file)
hpa_tis = h.get_HPA_Tissue()
extracted_dataset, tissues_names = check_for_file(hpa_tis, True)



def experiment(gene_set, all_enzymes, dataset, index_set, experiment_type = exp_type,  rts = random_test_size, rpa=recall_precision_threshold):     
    over_all_start = time.time()
    cr = create_random_sets(gene_set, dataset, all_enzymes, random_size=rts) #rts, check cofiguration file
    precision = rpa #rta, check cofiguration file
    col_zscore = {}

    # i.e index are tissue_names in tissue experiment, and cell_types in cell experiemnt
    for ind in index_set: 
        ml_score = RecallScore()
    
        for rand_num, ext_enz_set in cr.items():
            if experiment_type == "Tissue":
                glyco_enz_set_data = EnzymeData(dataset, ext_enz_set)        
                glyco_enz_set_data.add_parameters(ind)
                gnt = glyco_enz_set_data.get_gen_dataset() 
                glyco_enz_set_data.reset() 

            elif experiment_type == "Cell":
                print("This is cell experiment")
                gnt = make_df(ext_enz_set)
                
                if type(ind) == tuple:
                    gnt["Class"] = (gnt.index == ind).astype(int)
            
                else:
                    index_flat = gnt.index.map(lambda x: " ".join(map(str,x)))
                    gnt["Class"] = index_flat.str.contains(ind, regex=False)
                    
                
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
        
    
    #z_plot(col_zscore, plt_show=False, plt_save=False, title=gene_set)

    print_process("Whole Set", over_all_start)

    return col_zscore


def run(enz_sets, glyco_enz,extracted_data, filename):
    total_gr_zscore = defaultdict(list)
    for i, (gn_group, group_name) in enumerate(enz_sets.items()):
        gn_group_set = ""
        for gn in gn_group:
            gn_group_set += f"{gn} "  

        if list(gn_group)[0] != "GGTA1":
            print(gn_group_set)
            #total_gr_zscore[(group_name, gn_group_set)] = set_enz_experiment(gn_group, glyco_enz, extracted_data)
            total_gr_zscore[(group_name, gn_group_set)] = experiment(gn_group, glyco_enz, extracted_data, index_set=tissues_names)
            
        
    save_zdata(total_gr_zscore, filename)



gn_sets = GeneSet()
all_sets = gn_sets.get_sdbox_data()
glyco_enzymes = gn_sets.get_all_glyco_enz() 
r = run(all_sets, glyco_enzymes, extracted_dataset, tissue_filename)


### run the below commented line, if the hardware is able to support it. It runs, just comment out above line (r =)
## below commented runs the tissue experiment as single glyco enzyme based 
"""
gn_lis = []
gn_not_indata = []
for gn_set in all_sets:
  for gn in gn_set:
    temp_df = extracted_dataset.get(gn)
    if temp_df is not None:
      gn_lis.append(gn)
    else:
      gn_not_indata.append(gn)



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
    def __init__(self, glyco_enz, data, genes_nt_data, tis_indices):
        self.glyco_enz = sorted(glyco_enz)
        self.data = data
        self.gn_nt_data = genes_nt_data
        self.tis_names = tis_indices

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
                #total_gr_zscore[(f"Worker {worker_id}", gn)] = single_enz_experiment(test_gn, self.glyco_enz, self.non_gly_gn, self.data)
                total_gr_zscore[(f"Worker {worker_id}", gn)] = experiment(test_gn, self.glyco_enz, self.data, index_set=self.tis_names)

        return total_gr_zscore
    
def worker_function(worker_id, glyco_enzymes, data,gns_nt_data, t_ind, queue, result_dict):
    w = Workers(glyco_enzymes, data, gns_nt_data, t_ind)
    result_dict[worker_id] = w.worker(queue, worker_id)




if __name__ == "__main__":

    cpuCount = os.cpu_count()
    print("Number of CPUs in the system:", cpuCount, "\n")

    glyco_enzymes = sorted(glyco_enzymes)

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

    job_indices = make_job_indices(w, j, glyco_enzymes)
    print(job_indices)

    queue = multiprocessing.Queue()
    manager = multiprocessing.Manager()
    result_dict = manager.dict()


    for indx in job_indices:
        queue.put(glyco_enzymes[indx])


    jobs = []
    for i in range(1, p + 1):
        p = multiprocessing.Process(target=worker_function, args=(i, glyco_enzymes, extracted_dataset, gn_not_indata, tissues_names, queue, result_dict))
        jobs.append(p)
        p.start()


    for job in jobs:
        job.join()

    combined_results = {}
    for result in result_dict.values():
        for key, value in result.items():
            combined_results[key] = value 

    save_zdata(combined_results, f"Combined_{args.filename}")
"""