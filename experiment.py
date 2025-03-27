
import numpy  as np
import pandas as pd

from sklearn.metrics import precision_recall_curve
from sklearn.model_selection import StratifiedShuffleSplit, StratifiedKFold
from sklearn.linear_model import LogisticRegression

from imblearn.under_sampling import RandomUnderSampler

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
    params = ['stdev_floor','non_glycoenzyme_genesets','replicates']
    def __init__(self, data, score, **kwargs):
        self.data = data
        self.scorer = score
        self.minsd = kwargs.get('stdev_floor',0.1)
        self.nonglycosets = kwargs.get('non_glycoenzyme_genesets',20)
        self.replicates = kwargs.get('replicates',5)

    def zscore(self,x,randobs):
        std = np.std(randobs)
        if std < self.minsd:
            std = self.minsd
        mean = np.mean(randobs)
        return max(x-mean,0)/std

    def do_analysis(self,task,**kwargs):
        genes = task['geneset']
        samples = task['samples']

        rows = []
        for i in range(self.replicates):
            obs = []
            for df in self.data.create_dataframes(genes,samples,self.nonglycosets):
                X = df.drop(['Class'], axis = 1)
                y = df["Class"]
                score = self.scorer.score(X,y)
                obs.append(score)
            zscore = self.zscore(obs[0],obs[1:])
            rows.append(dict(geneset_name=task['geneset_name'],genes=",".join(sorted(genes)),topredict=task['topredict'],
                             replicate=i+1,task_id=kwargs['task_index'],
                             score=obs[0],zscore=zscore,randscores=obs[1:]))
        
        medianrow = dict(sorted(rows,key=lambda d: d['zscore'])[len(rows)//2].items())
        del medianrow['replicate']
        medianrow['score'] = round(medianrow['score'],3)
        medianrow['zscore'] = round(medianrow['zscore'],3)
        return medianrow
