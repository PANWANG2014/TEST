

## miRACLe

miRACLe: improving the prediction of miRNA-mRNA interactions by a random contact model

## Introduction
'miRACLe' is an R package developed to infer individual-specific miRNA-mRNA interactions (MMIs) from a single paired miRNA-mRNA expression profile by applying a random contact model. Evaluation by a variety of measures shows that fitting a sequence-based algorithm into the framework of miRACLe can significantly improve its prediction power, and the combination of miRACLe and the cumulative weighted context++ scores (CWCSs) from TargetScan consistently outperforms state-of-the-art methods in prediction accuracy, regulatory potential and biological relevance. The package is easy to apply and fast in computation. It typically requires less than 10 seconds of CPU time to complete the prediction for an individual sample on a laptop.

## Installation
```r
install.packages('devtools') #skip this step if it is already installed
library('devtools')
install_github('lookgene/miRACLe')
library('miRACLe')
```

## Quick start with an example
1. Let’s generate individual-specific miRNA-mRNA interactions using a single pair of miRNA-mRNA expression profile of HeLa cells.<br>

```r
 data(seqScore.Rdata) # to load sequence-based interaction score, default is 'TargetScan7_CWCS_cons'
 data(Test_data.Rdata) # to load test datasets
 mirExpr <- Test_HeLa_miRNA
 tarExpr <- Test_HeLa_mRNA
 final_output_ind <- miracle_ind(seqScore, mirExpr, tarExpr, OutputSelect = TRUE)
```

* The essential inputs for 'miracle_ind()' are `seqScore` (sequence-based interaction scores), `mirExpr`(the expression profile of miRNA), `tarExpr`(the expression profile of mRNA). Another input is optional: `OutputSelect`(logical variable, select “TRUE” to return the top 10 percent-ranked predictions by scores, and “FALSE” to return the whole prediction result. Default is TRUE.)<br>




2. If the expression data of multiple samples are provided, we can generate miRNA-mRNA interactions at both individual and population levels. To do this, type following lines.<br>

```r
 data(seqScore.Rdata) # to load sequence-based interaction score, default is 'TargetScan7_CWCS_cons'
 data(Test_data.Rdata) # to load test datasets
 mirExpr <- Test_DLBC_miRNA
 tarExpr <- Test_DLBC_mRNA
 sampleMatch <- Test_DLBC_sampleMatch
 sampleSelect = c("TCGA-FA-A4BB-01A-11R-A31S-13", "TCGA-FA-A4XK-01A-11R-A31S-13", "TCGA-FA-A6HN-01A-11R-A31S-13") # samples selected from the test dataset to analyze
 final_output <- miracle(seqScore, sampleMatch, mirExpr, tarExpr, samSelect = sampleSelect, exprFilter = 1, OutputSelect = TRUE)
final_output$Ind	#Individual-level result
final_output$Pop	#Population-level result
```

* The output of ‘miracle()’ is a list containing two elements. The first "Ind" contains the individual-level prediction results ranked by miracle scores, and the second "Pop" contains the population-level prediction result ranked by miracle scores.

* The essential inputs for 'miracle()' are `seqScore` (sequence-based interaction scores), `sampleMatch`(corresponding relationships between samples from miRNA expression data and mRNA expression data), `mirExpr`(the expression profile of miRNA), `tarExpr`(the expression profile of mRNA). Another three inputs are optional: `samSelect`(sample selection, users can select a subset of all samples to analyze, default is NULL, no selection applied), `exprFilter`(filter of expression profile, miRNAs/mRNAs that are not expressed in more than a given percentage of samples will be removed, default is 1), `OutputSelect`(logical variable, select “TRUE” to return the top 10 percent-ranked predictions by scores, and “FALSE” to return the whole prediction result. Default is TRUE.)<br>

## About seqScore

* The `seqScore` is a Z x 3 data.frame that contains sequence-based interaction scores for putative miRNA-mRNA pairs. These scores are originally obtained from [TargetSan v7.2](http://www.targetscan.org/vert_72/) (TargetScan7\_CWCS\_cons and TargetScan7\_CWCS), [DIANA-microT-CDS](http://diana.imis.athena-innovation.gr/DianaTools/index.php?r=microT_CDS/index) (DIANA\_microT\_CDS), [MirTarget v4]( http://mirdb.org/) (MirTarget4), [miRanda-mirSVR](http://www.microrna.org/microrna/getDownloads.do) (miRanda\_mirSVR) and compiled by the developers to fit the model. Default is **TargetScan7\_CWCS\_cons** and is provided in the package. The other scores can be downloaded [here](https://figshare.com/s/0b7c68cd5152da27a191).<br>

* User can also provide their own sequence-based interaction score, as long as the format of input file meets the requirements. Specifically, the first line must contain the label Names for mRNAs, miRNAs and their associated interaction scores. The remainder of the file contains RNA identifiers corresponding to those used in the expression files and the scores for each miRNA-mRNA pair. Note that the first column must contain identifiers for mRNAs, the second column must contain identifiers for miRNAs with the third column containing the associated scores.<br>

