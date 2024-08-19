from sklearn.metrics import classification_report
from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_recall_curve, roc_curve, auc
from sklearn.model_selection import StratifiedShuffleSplit, StratifiedKFold, cross_val_score
from imblearn.under_sampling import RandomUnderSampler
from sklearn.metrics import precision_recall_fscore_support
import numpy as np
from plots import * 
from app import *
import configparser


config = configparser.ConfigParser()
config.read('Sets_of_Enzymes_exp_config.ini')

recall_precision_at = config["ML Parameters"]["recall_precision_at"]
xy_test_resampling = config["ML Parameters"]["xy_test_resampling"]

  
def x_y_split(data):
    X = data.drop(['Class'], axis = 1)
    y = np.array(data["Class"])
    class_counts = np.bincount(y)

    count_class_0 = class_counts[0]
    count_class_1 = class_counts[1]
    strat_split = StratifiedShuffleSplit(n_splits=1, test_size=0.2)

    for train_index, test_index in strat_split.split(X, y):
      X_train, X_test = X.iloc[train_index], X.iloc[test_index]
      y_train, y_test = y[train_index], y[test_index]

    #Applied to both experiment
    if count_class_1 < 10000: 
      under = RandomUnderSampler(sampling_strategy=0.1)
      X_train, y_train = under.fit_resample(X_train, y_train)

    #Only for Single_Cell:
    #Downsample X_test for single cell
    
      under = RandomUnderSampler(sampling_strategy=1)
      X_test, y_test = under.fit_resample(X_test, y_test)
  
    return X_train, X_test, y_train, y_test


def cv_split(data):
    skf = StratifiedKFold(n_splits=5)
    X = data.drop(['Class'], axis = 1)
    y = np.array(data["Class"])
    for train_index, test_index in skf.split(X, y):
      X_train, X_test = X.iloc[train_index], X.iloc[test_index]
      y_train, y_test = y.iloc[train_index], y.iloc[test_index]
    
    return X_train, X_test, y_train, y_test


class ML_Parameters_Model:
  def __init__(self, data, model_name, model, X_train, X_test, y_train, y_test):
    self.X_train = X_train
    self.X_test = X_test
    self.y_train = y_train
    self.y_test = y_test
    self.data = data
    self.model_name = model_name
    self.model = model
    
  def set_up(self):
    self.model.fit(self.X_train, self.y_train)
    self.predictions = self.model.predict(self.X_test)
    return self.predictions
  
  def coef(self):
    coefficients = self.model.coef_
    feature_names = self.X_train.columns

    feat_coef = {}
    for feature, coef in zip(feature_names, coefficients[0]):
      feat_coef[feature] = coef

    return feat_coef
  
  def auc_val(self):
    roc_auc_val = roc_auc_score(self.y_test, self.predictions)
    return roc_auc_val
    
  def pr_curve(self):
    pred_proba = self.model.predict_proba(self.X_test)[:,1]

    precision, recall, _ = precision_recall_curve(self.y_test, pred_proba)
    p = Pre_Recall(recall, precision, self.model_name)
    interp_precision = p.interpolated_precision()

    return precision, recall, interp_precision

  def classi_map_data(self):
    pr,rec,f1_score, _ = precision_recall_fscore_support(self.y_test, self.predictions, average="macro")
    map_lis = [pr,rec,f1_score]
    return map_lis
       
  def validate(self):
    v_score = cross_val_score(self.model, self.X_train, self.y_train)
    v_result = (f"K-fold: %0.2f (+/- %0.2f)" % (v_score.mean(), v_score.std() * 2))
    return v_result

  

#def gen_ml_report(df, ml_names, tissue_name, x_y_split, rm_enz_exp=False):
def gen_ml_report(df, ml_names, tissue_name, cv_split, rm_enz_exp=False, recall_precision_at=recall_precision_at):
  
  X_train, X_test, y_train, y_test = x_y_split(df)
  #X_train, X_test, y_train, y_test = cv_split(df)

  pr_re_dic = {}
  ml_pred_lis = []
  classification_data = {}
  collect_auc_score = {}
  collect_cdf_score = {}
  tissue_pr_re = {}
  for model_name, model in ml_names.items():
    m = ML_Parameters_Model(df, model_name, model, X_train, X_test, y_train, y_test)
    pr = m.set_up()
    #print(classification_report(y_test, pr))
    auc_val = m.auc_val()
    precision, recall, interp_precision = m.pr_curve()
    clasi_m_data = m.classi_map_data()

    t = Pre_Recall(recall, interp_precision)
    re_at_90 = t.precision_at(0.9) #check configuration file
    collect_cdf_score[model_name] = {tissue_name: re_at_90}

    ml_pred_lis.append(pr)
    classification_data[model_name] = clasi_m_data
    collect_auc_score[model_name] = {tissue_name: auc_val}
    pr_re_dic[model_name] = {tissue_name: [recall, interp_precision]}
    tissue_pr_re[model_name] = [recall, precision]

    if rm_enz_exp == True:
      feat_weight = m.coef()
      srt = dict(sorted(feat_weight.items(), key=lambda x:abs(x[1])))
      drop_key = list(srt.keys())[0]
      print(srt.keys())
      print(srt.values())


  #heatmap_class_report(classification_data, tissue_name)
  #tissue_figure('Interp', tissue_pr_re, Pre_Recall,tissue_name)
  #confusion_table(ml_pred_lis, ml_names.keys(), y_test, tissue_name)

  if rm_enz_exp == False:
    return pr_re_dic, collect_cdf_score #collect_auc_score
  
  else:
    return pr_re_dic, collect_cdf_score, drop_key 