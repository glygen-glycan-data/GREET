#!/bin/sh
DATE=`date "+%Y-%m-%d"`
wget -O sandbox-"$DATE".json 'https://sandbox.glyomics.org/api/getEnzymeMappings.php?limiter=no_filter&val='
rm -f sandbox.json
ln -s sandbox-"$DATE".json sandbox.json
wget -O GTEx_gene_tpm-"$DATE".gct.gz 'https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz'
rm -f GTEx_gene_tpm.gct.gz
ln -s GTEx_gene_tpm-"$DATE".gct.gz GTEx_gene_tpm.gct.gz
wget -O GTEx_SampleAttributesDS-"$DATE".txt 'https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt'
rm -f GTEx_SampleAttributesDS.txt
ln -s GTEx_SampleAttributesDS-"$DATE".txt GTEx_SampleAttributesDS.txt

wget -O GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad 'https://storage.googleapis.com/adult-gtex/single-cell/v9/snrna-seq-data/GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad'



