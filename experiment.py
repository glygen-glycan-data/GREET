
import numpy  as np
import pandas as pd

from sklearn.metrics import precision_recall_curve
from sklearn.model_selection import StratifiedShuffleSplit, StratifiedKFold
from sklearn.linear_model import LogisticRegression

from imblearn.under_sampling import RandomUnderSampler

import math
from scipy.stats import norm, beta

class MLModel(object):
    def __init__(self,name,model,train_downsample=None,test_downsample=None):
        self.model = model
        self.name = name
        self.undersample_train = None
        if train_downsample is not None and train_downsample > 0:
            self.undersample_train = RandomUnderSampler(sampling_strategy=train_downsample)
        self.undersample_test = None
        if test_downsample is not None and test_downsample > 0:
            self.undersample_test = RandomUnderSampler(sampling_strategy=test_downsample)
    
    def fit(self,X_train,y_train):
        self.model.fit(X_train,y_train)

    def proba(self,X_test):
        return self.model.predict_proba(X_test)[:,1]

    @staticmethod
    def interpolated_precision(precision,recall):
        # assumes precision is sorted for recall high to low...
        # assert len(precision) == len(recall)
        # assert sorted(recall,reverse=True) == list(recall)
        interp_precision = []
        max_precision = precision[0]
        for i in range(0, len(precision)):
            if precision[i] > max_precision:
                max_precision = precision[i]
            interp_precision.append(max_precision)
        return interp_precision,recall

    def prcurve(self,y_test,pred_proba):
        precision, recall, _ = precision_recall_curve(y_test, pred_proba)
        return self.interpolated_precision(precision,recall)

class KFoldCVMLModel(MLModel):
    def __init__(self,kfoldcv_folds=5,**kwargs):
        super().__init__(**kwargs)
        self.skf = StratifiedKFold(n_splits=kfoldcv_folds)

    def fit_and_proba(self,X,y):
        y_test_out = []
        proba_out = []
        for train_index, test_index in self.skf.split(X, y):
            X_train, X_test = X.iloc[train_index], X.iloc[test_index]
            y_train, y_test = y.iloc[train_index], y.iloc[test_index]
            if self.undersample_train is not None:
                try:
                    X_train, y_train = self.undersample_train.fit_resample(X_train, y_train)
                except ValueError:
                    pass
            self.fit(X_train,y_train)
            if self.undersample_test is not None:
                try:
                    X_test, y_test = self.undersample_test.fit_resample(X_test, y_test)
                except ValueError:
                    pass
            proba = self.proba(X_test)
            y_test_out.append(y_test)
            proba_out.append(proba)
        y_test_out = pd.concat(y_test_out).reset_index(drop=True)
        proba_out = np.hstack(proba_out)
        return y_test_out, proba_out
            
class KSplitMLModel(MLModel):
    def __init__(self,ksplit_splits=5,test_split=0.2,**kwargs):
        super().__init__(**kwargs)
        self.strat_split = StratifiedShuffleSplit(n_splits=ksplit_splits, test_size=test_split)

    def fit_and_proba(self,X,y):
        y_test_out = []
        proba_out = []
        for train_index, test_index in self.strat_split.split(X, y):
            X_train, X_test = X.iloc[train_index], X.iloc[test_index]
            y_train, y_test = y.iloc[train_index], y.iloc[test_index]
            if self.undersample_train is not None:
                try:
                    X_train, y_train = self.undersample_train.fit_resample(X_train, y_train)
                except ValueError:
                    pass
            self.fit(X_train,y_train)
            if self.undersample_test is not None:
                try:
                    X_test, y_test = self.undersample_test.fit_resample(X_test, y_test)
                except ValueError:
                    pass
            proba = self.proba(X_test)
            y_test_out.append(y_test)
            proba_out.append(proba)
        y_test_out = pd.concat(y_test_out).reset_index(drop=True)
        proba_out = np.hstack(proba_out)
        return y_test_out,proba_out

class KFoldCVLRModel(KFoldCVMLModel):
    params = ['kfoldcv_folds','train_downsample','test_downsample']
    def __init__(self,model_params={},**kwargs):
        model_params.update(dict(max_iter=50000))
        super().__init__(name="K-Fold CV LogisticRegression",model=LogisticRegression(**model_params),**kwargs)

class KFoldCVLassoLRModel(KFoldCVMLModel):
    params = ['kfoldcv_folds','train_downsample','test_downsample']
    def __init__(self,model_params={},**kwargs):
        model_params.update(dict(penalty='l1', solver='liblinear', max_iter=50000))
        super().__init__(name="K-Fold CV Lasso LogisticRegression",model=LogisticRegression(**model_params),**kwargs)

class KSplitLRModel(KSplitMLModel):
    params = ['ksplit_splits','test_split','train_downsample','test_downsample']
    def __init__(self,model_params={},**kwargs):
        model_params.update(dict(max_iter=50000))
        super().__init__(name="K-Split LogisticRegression",model=LogisticRegression(**model_params),**kwargs)

class MLScore(object):
    def __init__(self,model):
        self.mlmodel = model

    def get_prcurve(self,X,y):
        y_test,proba = self.mlmodel.fit_and_proba(X,y)
        return self.mlmodel.prcurve(y_test,proba)

class RecallAtPrecision(MLScore):
    params = ['precision']
    def __init__(self,precision,**kwargs):
        self.precision = precision
        super().__init__(**kwargs)

    def score(self,X,y):
        precision,recall = self.get_prcurve(X,y)
        # assert sorted(precision) == list(precision)
        # assert sorted(recall,reverse=True) == list(recall)
        for p,r in zip(precision,recall):
            if p >= self.precision:
                return r
        return 0

class Experiment(object):
    params = ['pval_floor','nonzero_score_count','nonzero_stdev_floor','min_non_glycoenzyme_genesets','max_non_glycoenzyme_genesets','replicates']
    def __init__(self, data, score, **kwargs):
        self.data = data
        self.scorer = score
        self.minpv = kwargs.get('pval_floor')
        self.minnzsd = kwargs.get('nonzero_stdev_floor')
        self.nzsccnt = kwargs.get('nonzero_score_count')
        self.maxnonglycosets = kwargs.get('max_non_glycoenzyme_genesets')
        self.minnonglycosets = kwargs.get('min_non_glycoenzyme_genesets')
        self.replicates = kwargs.get('replicates')

    def do_nullmodel(self,task,**kwargs):
        ngenes = task['ngenes']
        samples = task['samples']
        replicate = task['replicate']

        nzobs = []
        obs = []
        for df in self.data.create_random_dataframes(ngenes,samples,self.maxnonglycosets):
            X = df.drop(['Class'], axis = 1)
            y = df["Class"]
            score = self.scorer.score(X,y)
            obs.append(score)
            if score > 0:
                nzobs.append(score)
                if len(obs) >= self.minnonglycosets and len(nzobs) >= self.nzsccnt:
                    break

        robust_mean,robust_stdev = self.robust(nzobs)

        if len(obs) == self.maxnonglycosets:
            pnz = len(nzobs)/len(obs)
        else:
            pnz = (len(nzobs)-1)/(len(obs)-1)
        # bparams = beta.fit(nzobs)
        # criticalpvals = [ (x,beta.sf(x,*bparams)) for x in (0.5,0.6,0.7,0.8,0.9,0.99) ]
        nzobs1 = [ min(x,0.99) for x in nzobs ]
        bparams1 = beta.fit(nzobs1,floc=0,fscale=1)
        criticalpvals1 = [ (x,beta.sf(x,*bparams1)) for x in (0.5,0.6,0.7,0.8,0.9,0.99) ]
        row = dict(key=(samples,ngenes),replicate=replicate,
                   n=len(obs),nnz=len(nzobs),pnz=pnz,mean=np.mean(obs),stdev=np.std(obs),
                   nzmean=robust_mean,nzstdev=robust_stdev,nzobs=sorted(nzobs),
                   betaparams=bparams1,cpv=criticalpvals1)
        return row

    def significance(self,score,ns):
        if score == 0.0:
            pval = (1-ns['pnz'])
            zscore = 0.0
        else:
            if score <= ns['nzmean']:
                zscore = 0
            else:
                zscore = (score-ns['nzmean'])/max(ns['nzstdev'],self.minnzsd)
            pval = ns['pnz']*beta.sf(min(score,0.99),*ns['betaparams'])
            # pval = ns['pnz']*norm.sf(zscore)
        val = dict(pval=pval,zscore=zscore)
        if pval < 1e-10:
            pval = 1e-10
        val['neg10logpval'] = abs(-10*math.log(pval,10))
        return val

    def robust(self,obs):
        med = np.median(obs)
        mad = np.median([ abs(x - med) for x in obs ])
        return med,1.482*mad

    def do_analysis(self,task,**kwargs):
        genes = task['geneset']
        samples = task['samples']
        nullstats = task['nullstats']

        rows = []
        for i in range(self.replicates):
            df = self.data.create_genebased_dataframe(genes,samples)
            X = df.drop(['Class'], axis = 1)
            y = df["Class"]
            score = self.scorer.score(X,y)
            rows.append(dict(ngenes=len(genes),genes=",".join(sorted(genes)),topredict=task['topredict'],
                             replicate=i+1,task_id=kwargs['task_index'],score=score))
        
        medianrow = dict(sorted(rows,key=lambda d: d['score'])[len(rows)//2].items())
        del medianrow['replicate']

        sigs = []
        for ns in nullstats:
            sigs.append(self.significance(medianrow['score'],ns))
        mediansig = sorted(sigs,key=lambda d: (d['pval'],-d['zscore']))[len(sigs)//2]
        medianrow.update(mediansig)
        for key in ('score','neg10logpval','zscore'):
            medianrow[key] = round(medianrow[key],3)
        for key in ('pval',):
            medianrow[key] = "%.3g"%(medianrow[key],)

        rows = []
        for gsn in task['geneset_names']:
            rows.append(dict(geneset_name=gsn,**medianrow))
        return rows
