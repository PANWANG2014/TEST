# miRACLe: inference of individual-specific miRNA-mRNA interactions
<br>

## Table of contents
1. <a href="#1">Introduction</a><br>
2. <a href="#2">Executing miRACLe</a><br>
	2.1 <a href="#3">Files required</a><br>
	2.2 <a href="#4">Script Execution</a><br>
3. <a href="#5">Reproduction of the miRACLe paper`s analyses</a><br>
4. <a href="#6">Examples of predicted results</a>



---
### <a name="1">1. Introduction</a>
miRACLe (<u>miR</u>NA <u>A</u>nalysis by a <u>C</u>ontact mode<u>L</u>) is a newly developed miRNA target prediction tool. It integrates both local and system-wide factors that affect the targeting process into a random contact model (Figure below). Assuming that the different types of RNAs in an individual sample are identical particles that move freely and contact each other randomly, miRACLe identifies target preferences by computing the relative probability of an effective contact occurring for each potential miRNA-mRNA combination. Unlike most of the existing methods, miRACLe does not require pre-selection of candidates (e.g. differentially expressed molecules), nor any learning process on sample data. It infers miRNA-mRNA interactions (MMIs) directly from high-dimensional data provided by a single sample. Therefore, it supports prediction of individual-specific miRNA targets.Our test suggests that on a laptop PC as an Intel Core i7-4712HQ with a 2.30 GHz CPU and 16 GB RAM, the source code implementation requires less than 30 seconds of CPU time to complete the prediction for one sample.<br>

![test](https://github.com/PANWANG2014/miRACLe/blob/master/Work%20flow%20of%20miRACLe.png)<br>

This figure shows the workflow of the miRACLe algorithm. (A) miRACLe is based on a random contact model in which the RNAs in a cell are assumed to be indifferent particles moving freely and randomly contacting each other. Different types of miRNAs compete for a common target mRNA, while different mRNAs compete for common miRNAs. Once the particles contact each other, only a certain percentage of the contacts lead to potent interactions under the condition that these hits are transformed into effective contacts marked by the formation of stable RISCs (RNA-induced silencing complex). (B) The process of targeting is determined by contact count and binding efficiency. The former correlates with the expression levels of the molecules, while the latter is quantified by the CWCS (cumulative weighted context score) of [TargetScan](https://elifesciences.org/articles/05005). With these parameters, miRACLe computes a score that measures the relative activity of each MMI.

---
### <a name="2">2. Executing miRACLe</a>
####  <a name="3">2.1 Files required</a>
In order to run the current version of miRACLe, the users should provide two data files describing the expression levels of each miRNA and mRNA for the same sample. And one additional file describing the correspondence of samples between the miRNA and mRNA data files. All files are tab-delimited ASCII text files and must comply with the following specifications:

1. **Input miRNA file** is organized as follows:<br>

	| miRNA | TCGA-05-4384-01A-01T-1754-13 | TCGA-05-4390-01A-02T-1754-13 | TCGA-05-4396-01A-21H-1857-13 | TCGA-50-5066-01A-01T-1627-13|
	| :-------------: |:-------------:| :-----:| :-----:|:-----:|
	| hsa-let-7a-5p | 19.0144019352802 | 16.2421134450244 | 19.2817036827951 | 18.072179718816 |
	| hsa-let-7a-3p | 7.31288295528436 | 6.20945336562895 | 7.83920378809694 | 6.2667865406949 |
	| hsa-let-7a-2-3p | 6.52356195605701 | 5.4594316186373 | 3.70043971814109 | 7.38370429247405 |
	| hsa-let-7b-5p | 16.9613140283806 | 15.5496944365338 | 17.8444238072634 | 16.9950715116303 |
	| hsa-let-7b-3p | 7.92481250360578 | 5.20945336562895 | 7.66533591718518 | 6.82017896241519 |

	The first line contains the labels Name followed by the identifiers for each sample in the dataset. <br>
	>Line format: `Name(tab)(sample 1 name)(tab)(sample 2 name) (tab) ... (sample N name)`<br>
	>Example: `miRNAName	sample_1	sample_2	...	sample_n`<br>

	The remainder of the file contains data for each of the miRNAs. There is one line for each miRNA. Each line contains the miRNA name and a value for each sample in the dataset.<br>

2. **Input mRNA file** is organized as follows:<br>

	| Gene | TCGA-05-4384-01 | TCGA-05-4390-01 | TCGA-05-4396-01 | TCGA-50-5066-01|
	| :-------------: |:-------------:| :-----:| :-----:|:-----:|
	| AARS | 10.70943229 | 11.69327441 | 12.42829508 | 11.04643008 |
	| AASDHPPT | 9.908138588 | 9.671621204 | 10.11131958 | 9.98327168 |
	| AASDH | 7.94715708 | 7.289783756 | 8.321654408 | 7.627425906 |
	| AASS | 9.964902673 | 7.775210583 | 9.172307988 | 5.950624713 |
	| AATF | 9.952503787 | 9.538021465 | 9.367084237 | 8.437501997 |

	The first line contains the labels Name followed by the identifiers for each sample in the dataset. <br>
	>Line format: `Name(tab)(sample 1 name)(tab)(sample 2 name) (tab) ... (sample N name)`<br>
	>Example: `GeneName	sample_a	sample_b	...	sample_m`<br>

	The remainder of the file contains data for each of the mRNAs. There is one line for each mRNA. Each line contains the mRNA name and a value for each sample in the dataset.<br>

3. **Sample-to-sample file** generally contains two columns, which shows the corresponding relationship of the sample identifiers in miRNA expression file and mRNA expression file. It also serves as a index to denote which samples we choose to analyze. It is organized as follows:<br>

	| miRNA | Gene | 
	| :-------------: |:-------------:| 
	| TCGA-50-5066-01A-01T-1627-13 | TCGA-50-5066-01 |
	| TCGA-05-4384-01A-01T-1754-13 | TCGA-05-4384-01 |
	| TCGA-05-4390-01A-02T-1754-13 | TCGA-05-4390-01 |
	| TCGA-05-4396-01A-21H-1857-13 | TCGA-05-4396-01 |

	The first line must contain the label Names for samples in each expression dataset with the first column for miRNA and second column for mRNA. <br>
	>Line format: `(sample name in miRNA file)(tab)(sample name in mRNA file)`<br>
	>Example: `sample_1	sample_a`<br>

	The remainder of the file contains sample identifiers used in the miRNA and mRNA expression files. There is one line for each sample. Each line contains the identifiers for that sample.<br>

#### <a name="4">2.2 Script Execution</a><br>
miRACLe is written in R. Thus the users first need to download and install the R software on the platform (refer to [R-project](https://www.r-project.org/) for details). [The code](https://github.com/PANWANG2014/miRTIGO/blob/master/miRTIGO_algorithm/miRTIGO.R) of miRACLe consists of three parts, namely, 'FUNCTIONS', 'DATA INPUT' and 'MAIN PROGRAM'. The users only need to focus on the 'DATA INPUT' part and fill in the relevant files described as follows:<br>
	
>178 mrna = as.matrix(read.table("`mrna_list.txt`", head = TRUE, sep = "\t"))<br>
>179 mirna = as.matrix(read.table("`mirna_list.txt`", head = TRUE, sep = "\t"))<br>
>180 CWCS_matrix1 = as.matrix(read.table("`wMRE_all.txt`", head = TRUE, sep = "\t"))<br>

Those three files serve as the basis to define the sequence matching scores between miRNAs and mRNAs, which are compiled from TargetScan and provided [here](https://github.com/PANWANG2014/miRTIGO/blob/master/miRTIGO_algorithm/Sequence%20matching%20scores.7z). 
	
>188 name\_cancer = as.matrix(read.table("`Sample-to-sample.txt`", head = TRUE, sep = "\t"))<br>
>189 mirna\_cancer = as.matrix(read.table("`Input miRNA expression.txt`", head = FALSE, sep = "\t"))<br>
>190 mrna\_cancer = as.matrix(read.table("`Input mRNA expression.txt`", head = FALSE, sep = "\t"))<br>

Those three files should be provided by the users, which contains the paired expression profiles of miRNAs and mRNAs of the samples that they are interested in. Note that the input miRNA/mRNA file should be transformed into a non-negative matrix, in order for the main program to execute correctly. To achieve optimal prediction on the RNA sequencing data, we strongly recommend that users provide log2 transformed normalized counts (e.g. RSEM or RPM) as the input for our program.

---
### <a name="5">3. Reproduction of the miRACLe paper`s analyses</a><br>
1. The codes to reproduce these experiments in the paper are written in R and should be executed in the corresponding software environment.<br> 
2. Generally, all these [codes](https://github.com/PANWANG2014/miRTIGO/tree/master/Source%20codes%20for%20analyses) are arranged into three parts as 'FUNCTIONS', 'INPUT DATA' and 'MAIN PROGRAM'. The users need to download and fill in the relevant input files before implementing corresponding analyses.<br>
3. Files required for the reproduction of the analyses can be broadly classified into three categories:<br>

* Files needed to execute the algorithms (basic data)<br>

	| Data file | Descriptions | 
	|:-------------:|:-------------| 
	| wMRE\_all.txt | _cumulative weighted context++ scores_ (weighted miRNA response elements, wMREs) between an miRNA and mRNA |
	| qMRE\_all.txt | number of target sites (quantitative miRNA response elements, qMREs) on one mRNA for one miRNA |
	| conserved\_qMRE.txt | number of conserved target sites (qMREs) on one mRNA for one miRNA|
	| mirna\_list.txt | miRNA identifiers used in the above files|
	| mrna\_list.txt | mRNA identifiers used in the above files|

	All these files are compiled from [TargetScan v7.0](http://www.targetscan.org/cgi-bin/targetscan/data_download.cgi?db=vert_70).

* Files needed to perform evaluation analyses (reference data)<br>
	* Experimentally confirmed MMIs<br>

	| Data file | Descriptions | validated MMI counts |
	|:-------------:|:-------------|:-----:|
	| V1 | [TarBase v7.0](http://carolina.imis.athena-innovation.gr/diana_tools/web/index.php?r=tarbasev8%2Findex) | 307,010 |
    | V2 | [miRTarbase v7.0](http://mirtarbase.mbc.nctu.edu.tw/php/index.php) | 380,639 |
	| V3 | high-confidence set compiled from [miRTarbase](http://carolina.imis.athena-innovation.gr/diana_tools/web/index.php?r=tarbasev8%2Findex), [miRecords](http://mirecords.umn.edu/miRecords) and [oncomiRDB](http://bioinfo.au.tsinghua.edu.cn/oncomirdb/) | 9642 |
    | V4 | [starBase v2.0](http://starbase.sysu.edu.cn/starbase2/index.php) | 10,028 |
	
    * Curated miRNA transfection experiments<br>

    | Data file | Descriptions |
	|:-------------:|:-------------|
	| T1 | It contains 106 unique miRNA transfections and was originally collected by [Li _et al_](https://academic.oup.com/bioinformatics/article/30/5/621/247621) |
	| T2 | It contains 23 unique miRNA transfections and was originally collected by [Gumienny R _et al_](https://academic.oup.com/nar/article/43/3/1380/2411948)|

	* Curated cancer-related miRNAs and genes<br>

	| Data file | Descriptions | Molecule counts |
	|:-------------:|:-------------|:-----:|
	| miRNA biomarkers | miRNAs obtained from [oncomiR](http://www.oncomir.org) that are simultaneously correlated with tumor development, tumor staging, tumor grade and patient survival | 288 |
	| oncomirs | high-confidence cancer associated miRNAs compiled from [MNDR v2.0](http://www.rna-society.org/mndr/) database | 399 |
    | Cancer genes | [COSMIC](https://cancer.sanger.ac.uk/census) database | 616 |<br>

    The reference data files are stored along with relevant source codes.

* Expression data files for relevant analysis<br>

	| Data file| Descriptions |
	|:-------------: |:-------------|
	| TCGA data | Paired miRNA-mRNA expression profiles from a total of 7,991 tumor samples belonging to 32 different cancer types |
	| NCI60 data | Paired miRNA-mRNA expression profiles from 59 NCI-60 cell lines |
    | Exosomal data | Paired miRNA-mRNA expression profiles of blood exosomes from prostate cancer patient and healthy person|
    | ESCC data | Paired miRNA-mRNA expression profiles from 119 ESCC (esophageal squamous cell carcinoma) patients and associated clinical information |

	The expression data files are provided in a compressed file [DATA.7z](https://www.dropbox.com/sh/aa0k59j39nftmo9/AAALFIiSpicrAEn8nRUJRjUWa?dl=0).<br> 

Detailed descriptions of the above files are provided in the **MATERIALS AND METHODS** of the paper.<br>
---
### <a name="6">4. Examples of predicted results</a><br>
miRACLe provides both population-level result and individual-level result, and it outputs the top 1% ranked MMIs in default. To show as an example, [Here](https://github.com/PANWANG2014/miRTIGO/blob/master/Predicted_results.7z) we provide the detailed population-level results of the top 5000 predicted MMIs for each of the 32 TCGA cancer types and the individual-level results of the top 5000 predictions for the NCI-60 cell panel.

---

