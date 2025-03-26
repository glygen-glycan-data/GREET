import configparser, argparse

from glycoenzymes import GlycoEnzymes
from dataio import *
from experiment import *

class GREETConfig(object):
    GTEX_TISSUE_RNASEQ = 1
    GTEX_CELLTYPE_SCRNASEQ = 2

    def __init__(self,config_file):
        self.config = configparser.ConfigParser()
        self.config.read(config_file)

    @staticmethod
    def convert_exp_type(s):
        assert s in ("GTEX_TISSUE_RNASEQ","GTEX_CELLTYPE_SCRNASEQ"), "Experiment type %s not valid"%(s,)
        return getattr(GREETConfig,s)

    def get_param(self,section,param,conversion=str):
        return conversion(self.config[section][param])
        
    def sandbox_groupby(self):
        return self.get_param("Enzymes","groupby")

    def sandbox_datafile(self):
        return self.get_param("Enzymes","datafile")

    def experiment_type(self):
        return self.get_param("Data","experiment_type",GREETConfig.convert_exp_type)

    def min_samples_per_sampletype(self):
        return self.get_param("Data","min_samples_per_sampletype",float)

    def non_glycoenzyme_genesets(self):
        return self.get_param("Parameters","non_glycoenzyme_genesets",int)
    
    def stdev_floor(self):
        return self.get_param("Parameters","stdev_floor",float)

    def replicates(self):
        return self.get_param("Parameters","replicates",int)

    def model_class(self):
        return self.get_param("MachineLearning","model_class",str)

    def kfoldcv_folds(self):
        return self.get_param("MachineLearning","kfoldcv_folds",int)

    def ksplit_splits(self):
        return self.get_param("MachineLearning","ksplit_splits",int)

    def test_split(self):
        return self.get_param("MachineLearning","test_split",float)

    def train_downsample(self):
        return self.get_param("MachineLearning","train_downsample",float)

    def test_downsample(self):
        return self.get_param("MachineLearning","test_downsample",float)

    def score_class(self):
        return self.get_param("MachineLearning","score_class",str)

    def precision(self):
        return self.get_param("MachineLearning","precision",float)

    def get_enzymes(self):
        return GlycoEnzymes(datafile=self.sandbox_datafile(),
                            groupby=self.sandbox_groupby())

    def get_data(self,datafile):
        if self.experiment_type() == self.GTEX_TISSUE_RNASEQ:
            dataio = GTExTissueRNASeq(self)
        elif self.experiment_type() == self.GTEX_CELLTYPE_SCRNASEQ:
            dataio = GTExCelltypeSCRNASeq(self)
        else:
            raise LookupError("Bad experiment type")
        dataio.read(datafile)
        return dataio

    def get_instance(self,clsname,*args,**kwargs):
        cls = eval(clsname)
        for p in cls.params:
            kwargs[p] = getattr(self,p)()
        return cls(*args,**kwargs)

    def get_model_instance(self):
        return self.get_instance(self.model_class())

    def get_score_instance(self,model):
        return self.get_instance(self.score_class(),model=model)

    def get_experiment_instance(self,data,score):
        return self.get_instance('Experiment',data,score=score)

    def get_experiment(self,data):
        model = self.get_model_instance()
        score = self.get_score_instance(model)
        return self.get_experiment_instance(data,score)
