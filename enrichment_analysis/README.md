# Using GlycoEnzOnto in a Custom Enrichment Analysis

`GlycoEnzOnto.gmt` contains the glycoenzyme pathway memberships in a tab-delimited format.  This file can be directly used in the [GSEA](https://www.gsea-msigdb.org/gsea/index.jsp), or can be processed in R to use other enrichment tools.  Brief instructions on how to use `GlycoEnzOnto.gmt` in various contexts is described below:

## GSEA GUI Application:

The [GSEA desktop application](https://www.gsea-msigdb.org/gsea/index.jsp) runs a [gene set enrichment analysis](https://www.pnas.org/doi/10.1073/pnas.0506580102) on a provided OMICs dataset.  Users can provide the `GlycoEnzOnto.gmt` file to GSEA to perform glycosylation pathway enrichment analyses on their own omics datasets.  Instructions follow below:

### 1. Download GSEA:

The GSEA program is available for [Mac](https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/gsea/software/desktop/4.3/GSEA_MacApp_4.3.2.app.zip), [Windows](https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/gsea/software/desktop/4.3/GSEA_Win_4.3.2-installer.exe), and [Linux](https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/gsea/software/desktop/4.3/GSEA_Linux_4.3.2.zip) platforms.  Follow the installation instructions, then start GSEA.

### 2. Load `GlycoEnzOnto.gmt` gene sets into GSEA:

After starting GSEA, navigate to the "Load data" tab on the left-hand side of the application.  Upload the `GlycoEnzOnto.gmt` file by using "Method 1: Browse for files" in the GUI.  Confirm that GlycoEnzOnto pathways have been successfully loaded into GSEA by checking the Object cache at the bottom right-hand part of the GUI.  A branch in the object tree called "Gene set databases" should have a child entity called "GlycoEnzOnto.gmt" present underneath.  Double-click the GlycoEnzOnto gene set name, and select "GeneSetMatrixViewer2" to view the GlycoEnzOnto pathways and their gene memberships.

### 3. Load omics datasets into GSEA: 

Users can upload their omics datasets using the same "Load data" interface to read `GlycoEnzOnto.gmt`.  While it is possible to provide a gene-by-sample matrix and sample phenotype labels to GSEA for computing gene-to-phenotype correlations, it is recommended that users provide more rigorous ranking criteria from other bioinformatics tools (limma or DESeq2 differential expression analyses, WGCNA for gene correlation network analysis, etc).  Once statistics are computed from these other methods, they can be provided to GSEA for a "GSEAPreranked" analysis in the [RNK](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Using_RNA-seq_Datasets_with_GSEA) format.

For instance, if one has RNA-Seq data from a cell line and would like to find differentially expressed glycosyation pathways in a drug-treated condition with respect to an undrugged condition, one could use DESeq2 to compute log2 fold changes in the drug-treated condition with respect to the untreated condition.  The log2 fold changes could then be converted to RNK format and loaded into GSEA for a GSEAPreranked Analysis.

* It is important that omics data either have gene names in the HGNC format, or that users provide an identifier-to-gene mapping file (called a chip file in GSEA) to ensure enrichment proceeds correctly.  For instance, if one is performing enrichment analyses on proteomics data where proteins are labeled with UniProt accessions, one should provide the corresponding chip file [Human\_UniProt\_Gene\_ID\_MSigDB.v7.2.chip](https://data.broadinstitute.org/gsea-msigdb/msigdb/annotations/human/Human_UniProt_Gene_ID_MSigDB.v7.2.chip).

### 4. Run the GSEA Preranked Analysis using GlycoEnzOnto Pathways:

Click on the "Run GSEAPreranked" option under "Tools" in the GSEA GUI.  Fill the following fields listed below to run the analysis:

* Required Fields:
  + Gene sets database: load GlycoEnzOnto.gmt by navigating to the "Local GMT/GMX" tab
  + Number of permutations: 1000
  + Ranked List: specify the RNK file to process for GSEA enrichment.
  + Collapse/Remap to gene symbols: 
    - Select "No_Collapase" if the RNK file contains HGNC symbols.  This will most likely be the case for any genomics datasets (RNA-Seq)
    - Select "Collapse" if the RNK file contains non-HGNC symbols.  Then, find the right chip mapping file from GSEA's website. Human mapping files can be found [here](https://data.broadinstitute.org/gsea-msigdb/msigdb/annotations/human/).
  + Chip platform:  If gene IDs need to be collapsed, upload the corresponding chip file here.

* Basic fields:
  + Analysis Name: Provide the name of the comparison being made in the RNK file being processed (drug vs no drug, for example)
  + Enrichment statistic: "weighted"
  + Max size: exclude larger sets = 500
  + Min size: exclude smaller sets = 2

### 5. Downloading results:

Once GSEA preranked is complete, an entry will appear in the GSEA reports window at the bottom left-hand corner of the GSEA GUI.  Click on the "Success" button to open an HTML page describing the enrichment results. 

## `fgsea` R Package

Another way to run GSEA analyses programmatically on GlycoEnzOnto pathways is to use the [`fgsea`](https://bioconductor.org/packages/release/bioc/html/fgsea.html) package.  It comes with a way to load pathways from a GMT file, and comes with plotting functionalities to visualize GSEA results.  A template for running `fgsea` on RNA-Seq differential expression results is shown below:

``` R
library(fgsea)
###############################
# Data loading and preparation:
###############################
#GlycoEnzOnto Pathways:
glycoenzonto_gmt_path<-'path/to/GlycoEnzOnto.gmt' ### Replace with path to GlycoEnzOnto.gmt
GlycoEnzOnto_pathways <- gmtPathways(glycoenzonto_gmt_path)
# Load some differential expression analysis data:
### For example, let's assume DESeq2 was used to produce log2 fold changes, and 
### you would like to use log2 fold change as the ranking statistic:
### NOTE: make sure HGNC gene names are the feature names
deseq2_results<-readRDS('./deseq2result.rds')
lfc_dta<-results(deseq2_results,contrast='statusdrugged')
###############
# Running GSEA
###############
glycoenzonto_gsea_results <- fgsea(
    GlycoEnzOnto_pathways, lfc_dta[,'log2FoldChange'], minSize=2, maxSize=500
    )
```

Please refer to the package's [vignette](https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html) for more specific use cases.
