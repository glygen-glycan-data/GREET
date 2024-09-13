# GREET
Glycosylation-Related Enzyme Expression Toolkit

# Scripts
* preprocessing_setup.py
* app.py
* ml_parameters.py
* plots.py


# PREPROCESSING
1. Classes:
    ###  GeneSet:
        Fetches enzyme mapping data from a remote API [SandBox].
        Organizes and filters the data, particularly focusing on Homo Sapiens.
        Provides methods to retrieve specific enzyme sets or all glyco-enzymes.
        Use to retrieve enzyme data filtered by specific criteria, such as species and gene name.

    ### FileStatus:
        Handles reading and writing to files.
        Processes gene expression data by tissues.
        Can filter and save the processed data.
        Use for file management related to gene expression data, including reading, writing, and extracting tissue-specific data.

    ### HPATissueClassifcation:
        Reads tissue classification data and gene-specific classification.
        Allows retrieval of tissue data and applies tissue-specific classifications.
        Use to apply tissue-specific classifications to gene expression data.

2. Functions:
    * samples(url): 
        Fetches sample names from a URL and organizes them by tissue.
    * extracting_data(extract_enzymes_tup, samples_names): 
        Extracts and organizes gene expression data based on specific tissues and enzymes.
    * setting_values_per_tissue(dataset, hpa_tissues={}, total_sample=False): 
        Aggregates and processes tissue data, possibly summing up samples.
    * adding_headers(tissue_dict_dataset, collapsed=False): 
        Generates headers based on tissue data, with an option to collapse similar tissues.
    * assign_class(tissue_headers, tissue_name): 
        Assigns a classification (binary) to tissue samples based on their names.
    * check_for_file(hpa_tissue=None, return_hpa=False): 
        Checks if a data file exists, processes it if found, or fetches and processes data otherwise.

4. Key URLs:
    * Enzyme Mapping API: https://edwardslab.bmcb.georgetown.edu/sandboxdev/api/getEnzymeMappings.php
    * Sample Attributes URL: https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt

5. Usage:
    ### Example Workflow
        Initialize GeneSet: Retrieve enzyme mappings and filter by Homo sapiens.
        Process File: Use FileStatus to read gene expression data and filter by tissues.
        Classify Data: Use HPATissueClassifcation to apply tissue-specific classifications to the processed data.
        Output: Save the processed data or use it for further analysis.

# Plots.py
1. Precision-Recall Curve Analysis:
    ### Pre_Recall: 
    This class handles precision-recall data and provides methods for calculating interpolated precision, plotting different types of precision-recall curves, and filtering the data based on recall or precision thresholds.

2. Visualization Functions:
    * monotonic_figure: Creates precision-recall curve plots for multiple models.
    * tissue_figure: Generates precision-recall curve plots for specific tissue types.
    * Make_Violin: A class to generate violin plots for gene expression data across tissues.
    * confusion_table: Generates confusion matrices for different models and visualizes them as heatmaps.
    * heatmap_class_report: Creates a heatmap for classification reports comparing multiple models.
    * histogram: Plots histograms of tissue AUC scores.
    * box_plot: Generates box plots for models' AUC scores.
    * cdf_plot: Plots cumulative distribution functions (CDF) of AUC scores.
    * z_plot: Plots Z-scores against recall scores for various models.
    * z_table: Generates and prints a table of Z-scores and recall scores filtered by specific criteria.


# ML_Parameters.py
1. Data Splitting Functions
    ### x_y_split(data)
        Purpose: Splits the dataset into training and testing sets while handling class imbalance.
        Process:
            * Separates features (X) and target variable (y).
            * Calculates the distribution of classes.
            * Uses StratifiedShuffleSplit to maintain class proportions in splits.
            * Applies RandomUnderSampler to balance the classes if necessary:
            * Training Set: Under-samples majority class to achieve a 10:1 ratio.
            * Testing Set: Under-samples to achieve a 1:1 ratio (specific to single-cell experiments).

        Returns: X_train, X_test, y_train, y_test

    ### cv_split(data)
        Purpose: Performs cross-validation splitting of the dataset.
        Process:
            * Uses StratifiedKFold to split the data into 5 folds while preserving class distribution.
        Returns: X_train, X_test, y_train, y_test for each fold.

2. Machine Learning Model Class
    ### ML_Parameters_Model
        Description: A class that encapsulates the training and evaluation process of a machine learning model.
        Initialization Parameters:
            * data: The complete dataset.
            * model_name: A string identifier for the model.
            * model: An instance of a scikit-learn estimator.
            * X_train, X_test, y_train, y_test: Training and testing datasets.
        
        Methods:
            * set_up():
                # Trains the model using X_train and y_train.
                # Predicts labels for X_test.
                # Returns: Predicted labels for the test set. 

            * coef():
                # Retrieves model coefficients (for linear models).
                # Maps feature names to their corresponding coefficients.
                # Returns: Dictionary of feature coefficients.

            * auc_val():
                # Calculates the ROC AUC score based on predictions.
                # Returns: ROC AUC score.

            * pr_curve():
                # Computes precision-recall curve data.
                # Utilizes predicted probabilities for positive class.
                # Returns: Precision, recall, and interpolated precision arrays.

            * classi_map_data():
                # Calculates precision, recall, and F1-score using macro averaging. 
                # Returns: List containing precision, recall, and F1-score.

            * validate():
                # Performs cross-validation and computes validation scores.
                # Returns: Formatted string summarizing cross-validation results.

3. Reporting Function
    ### gen_ml_report(df, ml_names, tissue_name, rm_enz_exp=False, recall_precision_at=recall_precision_at)
        Purpose: Automates the process of training multiple models and collecting their performance metrics.
        Parameters:
            * df: The dataset to be used.
            * ml_names: Dictionary mapping model names to their instances.
            * tissue_name: Identifier for the specific tissue or condition being analyzed.s
            * rm_enz_exp (optional): Boolean flag indicating whether to remove the least significant feature based on coefficients.
            * recall_precision_at (optional): Threshold for precision or recall, possibly loaded from a configuration file.

        Process:
            * Data splitting: Function used to split the data (e.g., x_y_split or cv_split).
            * Model Training and Evaluation:
            * Iterates over each model in ml_names.
            * Trains and evaluates each model using ML_Parameters_Model.
            * Collects various metrics including:
                # Precision, recall, and interpolated precision from precision-recall curves.
                # ROC AUC scores.
                # Classification metrics (precision, recall, F1-score).
                # Optionally identifies and records the least significant feature to be removed.
                
            * Result Compilation: Stores collected metrics in dictionaries for further analysis or visualization.

        Returns:
            * If rm_enz_exp is False: Dictionaries containing precision-recall data and recall at specified precision thresholds.
            * If rm_enz_exp is True: Same as above, plus the name of the least significant feature to be potentially removed.


# app.py
1. Defining Classifiers
    *   Purpose: This snippet defines a list of machine learning model names and corresponding classifier objects. It then maps each model name to its classifier using a dictionary ml_names_classi.
    *   Usage: This is useful for initializing and organizing multiple models that will be used in the analysis pipeline.

