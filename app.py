from ml_parameters import *
from plots import *
from preprocessing_setup import *
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from scipy import stats
import warnings
import statistics, math


# Ignore all warnings
warnings.filterwarnings("ignore")




names = ["Logistic Regression", 'Modified LR L1', "KNN"] # 'QDA']
classifiers = [
    LogisticRegression(max_iter=50000),
    LogisticRegression(penalty='l1', solver='liblinear', max_iter=50000),
    KNeighborsClassifier()]
   # QuadraticDiscriminantAnalysis()]

ml_names_classi = {}
for i in range(len(names)):
  ml_names_classi[names[i]] = classifiers[i]





def create_random_sets(gene_set, dataset, glyco_set):
  r_gen = defaultdict(set)
  non_glyco_set = []
  for gn in dataset:
    if gn not in glyco_set:
      non_glyco_set.append(gn)
  
  for r_i in range(20):
    get_random_enz = random.sample(non_glyco_set, len(gene_set))
    r_gen[r_i + 1] = get_random_enz     
  
  r_gen[21] = gene_set

  return r_gen



class EnzymeData:
  def __init__(self, dataset, gen_set):
    self.dataset = dataset
    self.gen_set = gen_set
    self.enz_dict = []
    self.head = []
    self.generate_dataset()

  def generate_dataset(self):
    gen_dataset = defaultdict(list)
    for gn in self.gen_set:
      gen_dataset[gn].append(self.dataset.get(gn)) 
    self.enz_dict, tissue_set = setting_values_per_tissue(gen_dataset, total_sample=True)
    self.head = adding_headers(tissue_set)
    self.df = pd.DataFrame(self.enz_dict, index= self.head)

  def add_parameters(self, t):
    tis_class = assign_class(self.head, t)
    self.df.insert(loc=0, column="Class", value=tis_class)

  def get_gen_dataset(self):
    return self.df

  def reset(self):
    pass 



class Report:
  def __init__(self, gn_dataset):
    self.gn_dataset = gn_dataset
    
  def execute_report(self, enz_set_num, ml_names=ml_names_classi, rm_exp=False):
    if rm_exp == True:
      collect_pr_re_dic, collect_cdf, drop_key = gen_ml_report(self.gn_dataset, ml_names, enz_set_num, x_y_split, rm_exp)
      self.gn_dataset = self.gn_dataset.drop(columns=[drop_key])
      print(i, drop_key)
      return collect_pr_re_dic, collect_cdf
  
    else:
      collect_pr_re_dic, collect_cdf = gen_ml_report(self.gn_dataset, ml_names, enz_set_num, x_y_split, rm_exp)
      return collect_pr_re_dic, collect_cdf


class RecallScore:
  def __init__(self):
    self.collect_cdf_scores = defaultdict(list)
    self.mon_pr = defaultdict(list)

  def extract_pr_scores(self, pr_data):
    for m_name, t_val in pr_data.items():
      self.mon_pr[m_name].append(t_val)

    return self.mon_pr
  
  def extract_cdf_scores(self, cdf_data):
    for m_name, v_lis in cdf_data.items():
      for au in v_lis.values():
        self.collect_cdf_scores[m_name].append(au)

    return self.collect_cdf_scores
  

class Scores:
  def __init__(self, cdf_data):
    self.cdf_scores = cdf_data

  def extract_std(self):
    for k, values in self.cdf_scores.items():
      std_value = np.std(values)
      #print(f"Standard Deviation {k}:", std_value)
    
    return std_value

  def extract_zscore(self):
    z_data = []
    for k, values in self.cdf_scores.items():
      u = statistics.mean(values[:-1])
      std = np.std(values[:-1])
      z = (values[-1] - u)/std
      if std = 0:
        z = 5
      if k == "Logistic Regression":
        z_data.append((z, values[-1]))

    return z_data
  
  def extract_tscore(self):
    t_data = []
    for k, values in self.cdf_scores.items():
      u = statistics.mean(values[:-1])
      std = np.std(values[:-1])
      n = 20
      se = std/ math.sqrt(n)
      t = (values[-1] - u)/se
      if k == "Logistic Regression":
        t_data.append((t, values[-1]))


    return t_data


class PRplots:
  def __init__(self, cdf_data, pr_data, precision, plt_show, plt_save):
    self.cdf_data = cdf_data
    self.pr_data = pr_data
    self.precision = precision
    self.plt_show = plt_show
    self.plt_save = plt_save

  def box_plt(self):
    box_plot(self.cdf_data, self.plt_show, self.plt_save)
  
  def normalized_plt(self):
    monotonic_figure('Regularized', self.pr_data, self.plt_show, self.plt_save, pr_score=self.precision)

  def original_plt(self):
    monotonic_figure('Original', self.pr_data, self.plt_show, self.plt_save)

  def cdf_plt(self):
    cdf_plot(self.cdf_data, self.plt_show, self.plt_save)
  
  def histo_plt(self):
    histogram(self.cdf_data, self.plt_show, self.plt_save)


def exp_violin_plt(data):
    v = Make_Violin(data)
    v.prepare_violin_data()
    v.violinplot




"""
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
"""