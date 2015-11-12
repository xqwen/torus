# Tutorial: Analysis of GEUVADIS Data 

In this tutorial, we go through the process of eGene discovery using the real GEUVADIS data. The data set contains multi-population eQTL data and is originally distributed from [this website](http://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/analysis_results/). We performed additional pre-processing steps that are documented in [Wen et al, 2015](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005176). 



## Step 1: Prepare Input files

Due to the size limitation, all the gzipped files are placed on an external server.

* Pre-computed single SNP Bayes factors by sbams_sslr: [geuv.summary.bf.gz](http://www-personal.umich.edu/~xwen/dap/data/geuv/geuv.summary.bf.gz)

* SNP position file: [geuv.snp.map.gz](http://www-personal.umich.edu/~xwen/dap/data/geuv/geuv.snp.map.gz)

* Gene TSS position file: [geuv.gene.map.gz](http://www-personal.umich.edu/~xwen/dap/data/geuv/geuv.gene.map.gz)

* Binding variants annotation file: [geuv.annot.gz](http://www-personal.umich.edu/~xwen/dap/data/geuv/geuv.annot.gz)


## Step 2: Enrichment Analysis

If the goal is solely for enrichment analysis and not for QTL discovery, one can run the following command to obtain point and uncertainty estimates for the enrichment parameters
```
torus -d geuv.summary.bf.gz --load_bf -smap geuv.snp.map.gz -gmap geuv.gene.map.gz -annot geuv.annot.gz  -est
```

## Step 3: QTL discovery

For QTL discovery accounting for annotations of 1) SNP distance to TSS and 2) Binding variants annotations, issue the following command 
```
torus -d geuv.summary.bf.gz --load_bf -smap geuv.snp.map.gz -gmap geuv.gene.map.gz -annot geuv.annot.gz  -qtl
```

## Step 4: Prepare Bayesian Prior for Fine-mapping Analysis
