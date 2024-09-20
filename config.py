import configparser, argparse



parser = argparse.ArgumentParser(description='Process output filename')

parser.add_argument('-f', "--filename", metavar='output filename', type=str, required=True,
                    help='name of the output file')

parser.add_argument('-c', "--config", metavar='config_file', type=str, required=True,
                    help='Path to the configuration file (.ini)')

parser.add_argument('-d', "--datafile", metavar='data_file', type=str, required=True,
                    help='Path to the data file')

parser.add_argument('-w', "--workers", metavar='Number of Workers', type=int,
                    help='an integer for the number of workers')

parser.add_argument('-n', "--nworker", metavar='nworkers', type=int,
                    help='an integer for the job number')

parser.add_argument('-p', "--processors", metavar='Number of processor', type=int,
                    help='an integer for the number of processors')


args = parser.parse_args()


w = args.workers
j = args.nworker


tissue_filename = args.filename
config_file = args.config
datafile = args.datafile

config = configparser.ConfigParser()
config.read(config_file)


threshold_count = float(config["Preprocessing"]["threshold_count"])
ngg_temp_file = int(config["Preprocessing"]["random_gene_set_size"])

  

random_test_size = int(config["Parameters"]["random_test_sample_size"])
exp_type = config["Parameters"]["experiment_type"]

recall_precision_threshold = float(config["ML Parameters"]["recall_precision_at"])
train_under_sample = float(config["ML Parameters"]["train_under_sample"])
test_under_sample = float(config["ML Parameters"]["test_under_sample"])

sd_threshold = float(config["Z Score"]["sd_less_than_threshold"])