

## miRACLe

miRACLe: improving the prediction of miRNA-mRNA interactions by a random contact model

## Introduction
'miRACLe' is an R package developed to infer individual-specific miRNA-mRNA interactions (MMIs) from a single paired miRNA-mRNA expression profile by applying a random contact model. Evaluation by a variety of measures shows that fitting a sequence-based algorithm into the framework of miRACLe can significantly improve its prediction power, and the combination of miRACLe and the cumulative weighted context++ scores from TargetScan consistently outperforms state-of-the-art methods in prediction accuracy, regulatory potential and biological relevance. Empirical test suggests that on a laptop Intel Core i7-4712HQ personal computer with a 2.30 GHz CPU and 16 GB of RAM, our source code implementation requires less than 10 seconds of CPU time to complete the prediction for an individual sample.

## Installation
```
**install.packages('devtools') (skip this step if it is already installed.)**
**library('devtools')**
**install_github('lookgene/miRACLe')**
**library('miRACLe')**
```

## Quick start with an example
Let’s generate individual-specific miRNA-mRNA interactions using a single pair of miRNA-mRNA expression profile of HeLa cells.<br>

```
**load(test_data.Rdata) # to load test datasets**
**mirExpr<- Test_HeLa_miRNA**
**tarExpr<- Test_HeLa_mRNA**
**final_output_ind <- miracle_ind(seqScore, mirExpr, tarExpr, OutputSelect = TRUE)**
**write.table(final_output_ind, file = "HeLa prediction result.txt", quote = FALSE, sep = "\t", row.names = FALSE)**
```

\#the essential inputs for 'miracle_ind()' are `seqScore` (sequence-based interaction scores), `mirExpr`(the expression profile of miRNA), `tarExpr`(the expression profile of mRNA). Another input is optional: `OutputSelect`(logical variable, select “TRUE” to return the top 10 percent-ranked predictions by scores, and “FALSE” to return the whole prediction result. Default is TRUE.)<br>

If the expression data of multiple samples are provided, we can generate miRNA-mRNA interactions at both individual and population levels. To do this, type following lines.<br>

```
**load(test_data.Rdata) # to load test datasets**
**mirExpr <- Test_DLBC_miRNA**
**tarExpr <- Test_DLBC_mRNA**
**sampleMatch <- Test_DLBC_sampleMatch**
**final_output <- miracle(seqScore, sampleMatch, mirExpr, tarExpr, samSelect = NULL, exprFilter = 1, OutputSelect = TRUE)**
**write.table(final_output$Ind, file = "Individual-level predictions.txt", quote = FALSE, sep = "\t", row.names = FALSE)	#Individual-level result**
**write.table(final_output$Pop, file = "Population-level predictions.txt", quote = FALSE, sep = "\t", row.names = FALSE)	#Population-level result**
```

\# the essential inputs for 'miracle()' are `seqScore` (sequence-based interaction scores), `sampleMatch`(corresponding relationships between samples from miRNA expression data and mRNA expression data), `mirExpr`(the expression profile of miRNA), `tarExpr`(the expression profile of mRNA). Another three inputs are optional: `samSelect`(sample selection, users can select a subset of all samples to analyze, default is no selection applied), `exprFilter`(filter of expression profile, miRNAs/mRNAs that are not expressed in more than a given percentage of samples will be removed, default is 1), `OutputSelect`(logical variable, select “TRUE” to return the top 10 percent-ranked predictions by scores, and “FALSE” to return the whole prediction result. Default is TRUE.)<br>


