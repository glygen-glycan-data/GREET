from ml_parameters import *
from plots import *
from preprocessing_setup import *
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from scipy.stats import zscore
import warnings


# Ignore all warnings
warnings.filterwarnings("ignore")

sam_url = "https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
samples_names = samples(sam_url)

url = "https://edwardslab.bmcb.georgetown.edu/sandboxdev/api/getEnzymeMappings.php?limiter=no_filter&val="
enzymes_dict_set, total_enzymes = sandbox_url_data(url)
extracted_dataset, non_genes = extracting_data(total_enzymes, samples_names)





names = ["Logistic Regression", 'Modified LR L1', "KNN"] # 'QDA']
classifiers = [
    LogisticRegression(max_iter=50000),
    LogisticRegression(penalty='l1', solver='liblinear', max_iter=50000),
    KNeighborsClassifier()]
   # QuadraticDiscriminantAnalysis()]

ml_names_classi = {}
for i in range(len(names)):
  ml_names_classi[names[i]] = classifiers[i]




class Executing_Experiment:
  def __init__(self):
    self.filter_dataset = defaultdict(list)
    
  def reset(self):
    self.gen_set = {}
    self.filter_dataset = defaultdict(list)

  def add_enz_set(self, gen_set):
    self.gen_set = gen_set

  def filtered_enz_dataset(self):
    if self.gen_set != None:
      for enz in self.gen_set:
        self.filter_dataset[enz].append(extracted_dataset.get(enz))
    
  def prep_dataframe(self):
    self.enz_dict, tissue_set = setting_values_per_tissue(self.filter_dataset, total_sample=True)
    tissues_names = list(tissue_set.keys())
    self.head = adding_headers(tissue_set)
    
    return tissues_names

  def setup_dataframe(self):
    self.df = pd.DataFrame(self.enz_dict, index= self.head)

  def adding_parameters(self, t):
    tis_class = assign_class(self.head, t)
    self.df.insert(loc=0, column="Class", value=tis_class)
  
  def execute_report(self, enz_set_num, rm_exp=False):
    if rm_exp == True:
      self.pr_re_dic, self.collect_cdf, drop_key = gen_ml_report(self.df, ml_names_classi, enz_set_num, x_y_split, rm_exp)
      self.df = self.df.drop(columns=[drop_key])
      print(i, drop_key)

    else:
      self.pr_re_dic, self.collect_cdf = gen_ml_report(self.df, ml_names_classi, enz_set_num, x_y_split)

  def set_pr_cdf_(self):
    self.collect_cdf_scores = defaultdict(list)
    self.mon_pr = defaultdict(list)

  def extract_pr_scores(self):
    for m_name, t_val in self.pr_re_dic.items():
      self.mon_pr[m_name].append(t_val)

  def extract_cdf_scores(self):
    for m_name, v_lis in self.collect_cdf.items():
      for au in v_lis.values():
        self.collect_cdf_scores[m_name].append(au)

  def extract_std(self):
    self.all_values = np.concatenate(list(self.collect_cdf_scores.values())) 
    std_value = np.std(self.all_values)
    print("Standard Deviation:", std_value)

  def extract_zscore(self):
    z_score = zscore(self.all_values)
    print("Z-score:", z_score)

  def exp_violin_plt(self):
    v = Make_Violin(self.df)
    v.prepare_violin_data()
    v.violinplot
  
  def exp_box_plt(self):
    box_plot(self.collect_cdf_scores)
  
  def exp_normalized_plt(self, precision_at):
    monotonic_figure('Regularized', self.mon_pr, precision_at)

  def exp_original_plt(self):
    monotonic_figure('Original', self.mon_pr)

  def exp_cdf_plt(self):
    cdf_plot(self.collect_cdf_scores)
  
  def exp_histo_plt(self):
    histogram(self.collect_cdf_scores)




  