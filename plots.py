from collections import defaultdict
from sklearn.metrics import confusion_matrix
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
import numpy as np
import os

class Pre_Recall:
  def __init__(self, recall, precision, label=None):
    self.recall = recall
    self.precision = precision
    self.label = label

  def recall_at(self, recall=0.85):
    pre_reg = []
    rec_reg = []
    reg_set_label = ""
    for i in range(0, len(self.precision)):
      if self.precision[i] >= recall:
        pre_reg.append(self.precision[i])
        rec_reg.append(self.recall[i])
        if self.precision[i] != 1.0:
          reg_set_label = self.label
    
    self.precision = pre_reg
    self.recall = rec_reg
    self.label = reg_set_label

  def precision_at(self, precision=0.9):
    pre_reg = []
    rec_reg = []
    for i in range(0, len(self.precision)):
      if self.precision[i] >= precision:
        pre_reg.append(self.precision[i])
        rec_reg.append(self.recall[i])

    self.pre_target = precision
    last_pre = pre_reg[0]
    last_re = rec_reg[0]

    if sorted(rec_reg)[-1] > 0.00:
      self.last_pre = last_pre
      self.last_re = last_re


    return(last_re)



  def interpolated_precision(self):
    interp_precision = []
    max_precision = self.precision[0]
    for i in range(0, len(self.precision)):
        if self.precision[i] > max_precision:
            max_precision = self.precision[i]
        interp_precision.append(max_precision)

    self.interp_precision = interp_precision

    return interp_precision

  def plot(self, plot_type=["Original", "Interp", "Regularized"]):
    if plot_type == "Original":
      plt.plot(self.recall, self.precision,  label=self.label)
    elif plot_type == 'Interp':
      plt.plot(self.recall, self.interp_precision, label=self.label)
    elif plot_type == 'Regularized':
      if hasattr(self, 'last_pre') and hasattr(self, 'last_re'):
        plt.plot([self.last_re, self.last_re], [0, self.last_pre], color='lightsteelblue', linestyle='--',linewidth=1.1)
        plt.axhline(y=self.pre_target, color='lightgrey', linestyle='--', linewidth=1.1)
        plt.plot(self.recall, self.precision, markevery=25,markersize=9, label=self.label)
        plt.annotate(self.pre_target, fontweight='light',xy=(0, self.pre_target), xytext=(-50, self.pre_target),
                     arrowprops=dict(facecolor='black', arrowstyle='->'), textcoords="offset points", ha='left', fontsize=10)



def monotonic_figure(line, plot_data_dict, plt_show, plt_save, pr_score=1.0):
  for m_name, v_dic in plot_data_dict.items():
    plt.figure(figsize=(10, 8), layout='tight')
    for val in v_dic:
      for tis, v in val.items():
        if line == 'Regularized':
          p = Pre_Recall(v[0], v[1], tis)
          p.precision_at(pr_score)
        elif line == 'Original':
          p = Pre_Recall(v[0], v[1])
        p.interpolated_precision()
        p.plot(line)

    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title(f'Precision-Recall Curve {m_name}')
    plt.grid(True)
    plt.legend(loc='best', bbox_to_anchor=(1.0, 1.0))

    if plt_show == True:
      plt.show()
      
    if plt_save == True:
      file_name = f"Monethesitic_{line}_{m_name}"
      directory = f".\\plot_figures\\{line}_Pr_Re Plots\\"
      os.makedirs(directory, exist_ok=True)
      plt.savefig(os.path.join(directory, file_name))
    
    plt.close()

def tissue_figure(line, plot_data_dict, pre_recall_class,t_name):
  plt.figure(figsize=(12, 10))
  for m_name, v_dic in plot_data_dict.items():
    for val in v_dic:
      p = pre_recall_class(v_dic[0], v_dic[1], m_name)
    p.interpolated_precision()
    p.plot(line)
  plt.xlabel('Recall')
  plt.ylabel('Precision')
  plt.title(f'Precision-Recall Curve {t_name}')
  plt.grid(True)
  plt.legend(loc='best', bbox_to_anchor=(1.0, 1.0))
  #directory = f".\\plot_figures\\{gene_key[0]}\\{line}_Pr_Re Plots\\"
  #file_name = f"{line}_{t_name}"
  #os.makedirs(directory, exist_ok=True)
  #plt.legend(loc='best', bbox_to_anchor=(1.0, 1.0))
  #plt.savefig(os.path.join(directory, file_name))
  plt.close()



class Make_Violin:
  def __init__(self, data):
    self.data = data

  def prepare_violin_data(self):
    vio_data = {}
    for gn, values in self.data.items():
      if gn != 'Class':
        tem_data = defaultdict(list)
        for tissue_name, val in values.items():
          tem_data[tissue_name].append(val)

        srt = dict(sorted(tem_data.items()))
        vio_data[gn] = srt
    
    self.data = vio_data

  def violinplot(self):
    global gene_key
    sorted_dict = dict(sorted(self.data.items()))
    plt.figure(figsize=(15, 10))
    for i, (gn, val) in enumerate(sorted_dict.items(), start=1):
      gene_key = re.split(r'\d+$', gn)
      plt.subplot(len(self.data), 1, i)
      sns.violinplot(val, saturation=0.50, fill=False, inner="quart", density_norm='width')
      if i < len(self.data):
        plt.xticks([])
      else:
        plt.xticks(rotation=90)
        plt.xlabel('Tissue')

      plt.ylabel(gn)
    
    plt.show()
    #directory = ".\\plot_figures\\violin_plot\\"
    #file_name = f"{gene_key[0]}"
    #os.makedirs(directory, exist_ok=True)
    #plt.tight_layout()
    #plt.savefig(os.path.join(directory, file_name))
    plt.close()

def confusion_table(model_predictions, models_names,  data, tissue_name):
  m_dic = {}
  num = 1
  plt.figure(figsize=(10, 8))
  for l in model_predictions:
    CM = confusion_matrix(data, l)
    TN = CM[0][0]
    FN = CM[1][0]
    TP = CM[1][1]
    FP = CM[0][1]
    ACC = ((TP+TN)/(TP+FP+FN+TN) * 100)
    m_dic[num] = [TP, TN, FP, FN, ACC]
    num +=1

  m_dic = pd.DataFrame(m_dic)
  m_dic.columns = models_names
  m_dic = m_dic.transpose()
  xlabels = ["TP", "TN", "FP", "FN", "ACC"]
  graph = sns.heatmap(m_dic, xticklabels=xlabels ,cmap = "Blues", annot = True)
  plt.title(tissue_name)
  directory = f".\\plot_figures\\{gene_key[0]}\\Cofusion table\\"
  #file_name = f"{tissue_name}"
  #os.makedirs(directory, exist_ok=True)
  #plt.savefig(os.path.join(directory, file_name))
  plt.close()

def heatmap_class_report(compare_data, t_name):
  val_lis = []
  for ml_name, val in compare_data.items():
      val_lis.append(val)

  name = f"Classification Report {t_name}"
  Class_name = ["Precision", "Recall", "F1 Score"]

  ax = sns.heatmap(val_lis, xticklabels= Class_name,cmap = "Blues", annot = True)
  ax.set_yticklabels(compare_data.keys(), rotation='horizontal')
  plt.title(name)
  #directory = f".\\plot_figures\\{gene_key[0]}\\Classification Table\\"
  #file_name = f"heatmap_{t_name}"
  #os.makedirs(directory, exist_ok=True)
  #plt.savefig(os.path.join(directory, file_name))
  plt.close()


def histogram(data_dict, plt_show, plt_save):
    plt.figure(figsize=(10, 6))
    for model_name, val in data_dict.items():
        sns.histplot(val, bins=35, alpha=0.7, label=model_name, element='poly')
    plt.title('Histogram of Tissue AUC')
    plt.xlabel('Tissue AUC Score')
    plt.ylabel('Frequency')
    plt.legend()
    plt.tight_layout()
    if plt_show == True:
      plt.show()
      
    if plt_save == True:
      directory = f".\\plot_figures\\{gene_key[0]}\\Histo plots\\"
      file_name = f"histogram"
      os.makedirs(directory, exist_ok=True)
      plt.savefig(os.path.join(directory, file_name))
    
    plt.close()

def box_plot(data_dict,plt_show, plt_save):
    plt.figure(figsize=(10, 6))
    sns.boxenplot(data_dict, saturation=0.5)
    plt.title('Models Histogram  AUC')
    plt.xlabel('Tissue AUC Score')
    plt.ylabel('Frequency')
    plt.tight_layout()
    if plt_show == True:
      plt.show()
      
    if plt_save == True:
      directory = f".\\plot_figures\\{gene_key[0]}\\Histo plots\\"
      file_name = f"box_plot"
      os.makedirs(directory, exist_ok=True)
      plt.savefig(os.path.join(directory, file_name))
    
    plt.close()



def cdf_plot(data, plt_show, plt_save):
  plt.figure(figsize=(10, 8))
  for model_name, au_lis in data.items():
    au_lis_sorted = np.sort(au_lis)
    y = np.arange(len(au_lis_sorted)) / float(len(au_lis_sorted))
    plt.plot(y, au_lis_sorted, marker='o', label=model_name)

  plt.xlabel('Percent')
  plt.ylabel('Score')
  plt.title('(CDF) of AUC Scores')
  plt.grid(True)
  plt.legend()
  if plt_show == True:
    plt.show()
    plt.close()
  if plt_save == True:
    directory = f".\\plot_figures\\{gene_key[0]}\\Histo plots\\"
    file_name = "CDF"
    os.makedirs(directory, exist_ok=True)
    plt.savefig(os.path.join(directory, file_name))
  
  plt.close()


def z_plot(data, plt_show:bool, plt_save:bool, all_set = False, title = ""):
  z_sc = []
  pr_s = []
  t_label = {}
  
  for k_tis, z_pr in data.items():
    if all_set == True:
      for gn_set, z_pr_value in z_pr.items():
        z_sc.append(z_pr_value[0][0])
        pr_s.append(z_pr_value[0][1])
        t_label[z_pr_value[0]] = gn_set
    
    else:
      z_sc.append(z_pr[0][0])
      pr_s.append(z_pr[0][1])
      t_label[z_pr[0]] = k_tis

  plt.figure(figsize=(12, 10))
  plt.scatter(z_sc, pr_s)


  for x,y in zip(z_sc, pr_s):
    if x > 4:
      label = t_label.get((x,y))
      plt.annotate(label, (x,y), textcoords="offset points", xytext=(0,10), ha='center')
      #if all_set == True:
      #  print(k_tis[1])
      #  plt.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0))
        

  plt.xlabel("Z-Score")
  plt.ylabel("Recall")
  plt.grid()
  if title != "":
    plt.title(title)
  if plt_show == True:
    plt.show()
  if plt_save == True:
    directory = f".\\plot_figures\\Z_plots\\"
    if all_set == True:
      file_name = "Combine_zscore"
    else:
      file_name = f"{title}"

    os.makedirs(directory, exist_ok=True)
    plt.savefig(os.path.join(directory, file_name))
  
  plt.close()

def z_table(data, plt_show=False, plt_save=False):
  restricted_zlabels = {}
  for k_tis, z_pr in data.items():
    for gn_set, z_pr_value in z_pr.items():
      z_x = z_pr_value[0][0]
      recal_y = z_pr_value[0][1]
      if z_x > 0:#4 and recal_y > 0.5:
        restricted_zlabels[k_tis[0]] = [k_tis[1], gn_set, z_x, recal_y]

  restricted_zlabels = pd.DataFrame(restricted_zlabels)
  restricted_zlabels.columns = restricted_zlabels.keys()
  restricted_zlabels = restricted_zlabels.transpose()
  columns = ["Genes", "Tissue", "Z-Score", "Recall-Score"]
  restricted_zlabels.columns = columns
  restricted_zlabels = restricted_zlabels.sort_values(by=["Tissue"])
  print(restricted_zlabels)


