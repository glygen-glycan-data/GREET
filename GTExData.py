#!.venv/bin/python

import os, sys, gzip, csv
import pandas as pd
import numpy as np
import time

starts = []
def tic():
    starts.append(time.time())

def toc(msg="Elapsed"):
    elapsed = time.time()-starts[-1]
    starts.pop(-1)
    print(f"{msg}: {elapsed:.2f}",file=sys.stderr)

class GTExData:
    tpmdataext='_gene_tpm.gct.gz'
    sampledataext='_SampleAttributesDS.txt'
    bindataext='.bin'

    def __init__(self,filename=None,ensgid=False,force=False):
        self.ensgid = ensgid
        if filename.endswith(self.bindataext) and os.path.exists(filename):
            self.readbin(filename)
        elif not force and filename.endswith(self.tpmdataext) and os.path.exists(filename.replace(self.tpmdataext,self.bindataext)):
            self.readbin(filename.replace(self.tpmdataext,self.bindataext))
        else:
            assert filename.endswith(self.tpmdataext)
            self.tmpdata = filename
            self.sampledata = filename.replace(self.tpmdataext,self.sampledataext)
            self.bindata = filename.replace(self.tpmdataext,self.bindataext)
            datadir = os.path.split(filename)[0]
            self.read(datadir)
            self.writebin(datadir)

    def read(self,datadir='data'):
        self.tpmfilename = os.path.join(datadir,self.tpmdata)
        self.samplefilename = os.path.join(datadir,self.sampledata)
        assert(os.path.exists(self.tpmfilename))
        assert(os.path.exists(self.samplefilename))
        tic()
        self.readsamplemetadata()
        toc("Read sample metatdata")
        self.df = self.getdataframe()

    def readsamplemetadata(self):
        self.samplemd = dict()
        for row in csv.DictReader(open(self.samplefilename), delimiter="\t"):
            sid = row["SAMPID"]
            tissue = row["SMTSD"]
            self.samplemd[sid] = dict(tissue=tissue)

    def get_sample_md(self,sid,mdkey):
        return self.samplemd[sid].get(mdkey,None)

    def get_sample_tissue(self,sid):
        return self.get_sample_md(sid,'tissue')

    def get_all_tissues(self):
        return self.tissues

    def getdataframe(self):
        tic()
        data = dict()
        h = gzip.open(self.tpmfilename,'rt')
        # Burn top two lines
        next(h);next(h)       
        samples = None
        for line in h:
            if not samples:
                samples = line.split()[2:]
                tissues = [ self.get_sample_tissue(s) for s in samples ]
                data['Sample'] = samples
                data['Tissue'] = tissues
                continue
            ensgid,genename,rest = line.split(None,2)
            if self.ensgid:
                data[ensgid] = np.array(rest.split(), dtype=np.float32)
            else:
                data[genename] = np.array(rest.split(), dtype=np.float32)
        h.close()
        toc("Read data I/O")
        tic()
        df = pd.DataFrame(data=data)
        toc("make data-frame")
        tic()
        df.set_index('Sample',inplace=True)
        toc("set index")
        return df

    def writebin(self,datadir='data'):
        self.df.reset_index(inplace=True)
        self.df.to_feather(os.path.join(datadir,self.bindata))
        self.df.set_index('Sample',inplace=True)

    def readbin(self,filename):
        self.df = pd.read_feather(filename)
        self.df.set_index('Sample',inplace=True)
        self.tissues = set(self.df['Tissue'])
        self.samplemd = dict(map(lambda t: (t[0],dict(tisssue=t[1])),zip(self.df.index,self.df['Tissue'])))

    def genes(self):
        return set(self.df.columns)

    def restrict(self,genes=None,tissue=None,nottissue=None,includetissue=False):
        if includetissue:
            genes = ['Tissue'] + list(genes)
        else:
            genes = list(genes)
        if tissue is not None:
            tissuesamp = (self.df['Tissue'] == tissue)
        elif nottissue is not None:
            tissuesamp = (self.df['Tissue'] != nottissue)
        else:
            tissuesamp = (self.df['Tissue'] != "")
        return self.df[tissuesamp][genes]

if __name__ == '__main__':

    genes = set(sys.argv[1:])
    tic()
    gtd = GTexData()
    toc("Constructor")
    print(gtd.restrict(genes=genes,nottissue='Liver',includetissue=True))
