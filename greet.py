#!.venv/bin/python

import sys, os, copy
from collections import defaultdict

os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'

import warnings
warnings.simplefilter("ignore", RuntimeWarning)

import argparse
from distproc import DistributedProcessing as dp

parser = argparse.ArgumentParser(description='Process output filename')

parser.add_argument('-c', "--config", metavar='config_file', type=str, required=True,
                    help='Path to the configuration file (.ini)')

parser.add_argument('-d', "--datafile", metavar='data_file', type=str, required=True,
                    help='Path to the data file')

parser.add_argument('-o', "--outfile", metavar='output_file', type=str, default=None,
                    help='Path to the output file')

parser.add_argument('-v', '--verbose', action = 'store_true', default = False,
                    help = 'Verbose logging.')

dp.add_arguments(parser,'workers')

args = parser.parse_args()
workers = dp.parse_args(parser)

from config import GREETConfig
config = GREETConfig(args.config)

data = config.get_data(args.datafile)
experiment = config.get_experiment(data)

proc = dp.stage_process_init(workers,args.verbose,[experiment.do_nullmodel,experiment.do_analysis])

gslens = [ len(geneset) for gsnames,geneset in data.glyco_enzyme_genesets() ]
mingslen = min(gslens)
maxgslen = max(gslens)

tasks = []
for name,samples in sorted(data.topredict(),key=lambda t: t[1]):
    for i in range(mingslen,maxgslen+1):
        for j in range(config.nullmodel_replicates()):
            tasks.append(dict(samples=samples,ngenes=i,replicate=j+1))

nullstats = defaultdict(list)
for row in proc.stage_process(0,tasks):
    print(row)
    sys.stdout.flush()
    nullstats[row['key']].append(copy.copy(row))

tasks = []
for gsnames,geneset in data.glyco_enzyme_genesets():
    for name,samples in data.topredict():
        tasks.append(dict(geneset_names=gsnames,geneset=geneset,
                          topredict=name,samples=samples,
                          nullstats=nullstats[samples,len(geneset)]))

output = sys.stdout
if args.outfile is not None:
    output = open(args.outfile,'w')

headers = """
  GeneSet Genes SampleType Recall ZScore p-Value Neg10LogPVal
""".split()
print("\t".join(headers),file=output)

for rows in proc.stage_process(1,tasks):
    pval = None
    for row in rows:
        print("\t".join([ str(row.get(h,"")) for h in "geneset_name genes topredict score zscore pval neg10logpval".split()]),file=output)
    output.flush()
if args.outfile is not None:
    output.close()

proc.stage_process_finish()
