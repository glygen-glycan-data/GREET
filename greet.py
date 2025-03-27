#!.venv/bin/python

import sys, os

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

dp.start_if_worker(workers,experiment.do_analysis)

tasks = []
for gsname,geneset in data.glyco_enzyme_genesets():
    for name,samples in data.topredict():
        tasks.append(dict(geneset_name=gsname,geneset=geneset,
                          topredict=name,samples=samples))
        # print(tasks[-1])

output = sys.stdout
if args.outfile is not None:
    output = open(args.outfile,'w')
logtempl = "worker_id: %(hostname)s:%(worker_index)s task_id: %(task_index)s runtime: %(runtime)s progress: %(progress)s remaining: %(remaining)s hrs"
headers = """
  GeneSet Genes SampleType Recall ZScore
""".split()
print("\t".join(headers),file=output)
for row in dp.process(workers,experiment.do_analysis,tasks,args.verbose,logtempl=logtempl):
    print("\t".join([ str(row.get(h,"")) for h in "geneset_name genes topredict score zscore".split()]),file=output)
    output.flush()
if args.outfile is not None:
    output.close()
