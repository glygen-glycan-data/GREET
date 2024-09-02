# GREET
Glycosylation-Related Enzyme Expression Toolkit

# Scripts
preprocessing_setup.py
app.py
ml_parameters.py
plots.py


# PREPROCESSING
1. Classes:
    ## GeneSet:
        Fetches enzyme mapping data from a remote API [SandBox].
        Organizes and filters the data, particularly focusing on Glyco Enzymes.
        Provides methods to retrieve specific enzyme sets or all glyco-enzymes.

    ## FileStatus:
        Handles reading and writing to files.
        Processes gene expression data by tissues.
        Can filter and save the processed data.

    ## HPATissueClassifcation:
        Reads tissue classification data and gene-specific classification.
        Allows retrieval of tissue data and applies tissue-specific classifications.

2. Usage:
    ## GeneSet Class: Use to retrieve enzyme data filtered by specific criteria, such as species and gene name.
    ## FileStatus Class: Use for file management related to gene expression data, including reading, writing, and extracting tissue-specific data.
    ## HPATissueClassifcation Class: Use to apply tissue-specific classifications to gene expression data.
    ## Helper Functions: These are utility functions that assist with fetching, processing, and classifying gene expression data.

3. Functions:
    ## samples(url): Fetches sample names from a URL and organizes them by tissue.
    ## extracting_data(extract_enzymes_tup, samples_names): Extracts and organizes gene expression data based on specific tissues and enzymes.
    ## setting_values_per_tissue(dataset, hpa_tissues={}, total_sample=False): Aggregates and processes tissue data, possibly summing up samples.
    ## adding_headers(tissue_dict_dataset, collapsed=False): Generates headers based on tissue data, with an option to collapse similar tissues.
    ## assign_class(tissue_headers, tissue_name): Assigns a classification (binary) to tissue samples based on their names.
    ## check_for_file(hpa_tissue=None, return_hpa=False): Checks if a data file exists, processes it if found, or fetches and processes data otherwise.


4. Key URLs:
    ## Enzyme Mapping API: https://edwardslab.bmcb.georgetown.edu/sandboxdev/api/getEnzymeMappings.php
    ## Sample Attributes URL: https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt

5. Example Workflow:
    ## Initialize GeneSet: Retrieve enzyme mappings and filter by Homo sapiens.
    ## Process File: Use FileStatus to read gene expression data and filter by tissues.
    ## Classify Data: Use HPATissueClassifcation to apply tissue-specific classifications to the processed data.
    ## Output: Save the processed data or use it for further analysis.