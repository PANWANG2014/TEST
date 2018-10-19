# miRTIGO: miRNA Target Inference by modeling Global molecular cOmpetition
<br>

## Table of contents
1. <a href="#1">Introduction</a><br>
2. <a href="#2">Purpose of the algorithm</a><br>
3. <a href="#3">Executing miRTIGO</a><br>
	3.1 <a href="#4">Files required</a><br>
	3.2 <a href="#5">Script Execution</a><br>
4. <a href="#6">Reproduction of the miRTIGO paper`s experiments</a><br>
5. <a href="#7">Examples of predicted results</a>



---
### <a name="1">1. Introduction</a>
miRTIGO is a novel miRNA target predictor, designed to identify sample-specific miRNA targets by  analyzing the same-sample miRNA-mRNA expression profiles. It decomposes the biological procedure of miRNA targeting into two independent stages: contacting and binding. By approximating the contacting stage by a random contact model, and the efficiency of the binding stage by sequence matching scores of RNAs, miRTIGO models endogenous RNA competition explicitly on a global scale and outputs an miRTIGO signature matrix to measure the relative activity of each individual miRNA-mRNA interaction (MMI).The figure below illustrates the rationale and workflow underlying miRTIGO.<br>

![test](https://github.com/PANWANG2014/miRTIGO/blob/master/Figure%201.jpg)<br>

---
### <a name="2">2. Purpose of this algorithm</a>
miRTIGO is developed to help researchers with insufficient programing skills to efficiently and accurately identify context-specific targets for miRNAs to explore the biological mechanisms of the miRNA-mediated post-transcriptional regulatory network when sample size is limited (e.g. rare tumor samples) or even in individual samples (e.g. single cells).

---
### <a name="3">3. Executing miRTIGO</a>
####  <a name="4">3.1 Files required</a>
In order to run the current version of miRTIGO, the users should provide two data files describing the expression levels of each miRNA and mRNA in the same sample, respectively. And one additional file describing the correspondence of samples between the miRNA and mRNA data files. All files are tab-delimited ASCII text files and must comply with the following specifications:

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

#### <a name="5">3.2 Script Execution</a><br>
miRTIGO is written in R. Thus the users first need to download and install the R software on the platform (refer to [R-project](https://www.r-project.org/) for details). [The code](https://github.com/PANWANG2014/miRTIGO/blob/master/miRTIGO_algorithm/miRTIGO.R) of miRTIGO consists of three parts, namely, 'FUNCTIONS', 'DATA INPUT' and 'MAIN PROGRAM'. The users only need to focus on the 'DATA INPUT' part and fill in the relevant files described as follows:<br>
	
>178 mrna = as.matrix(read.table("`mrna_list.txt`", head = TRUE, sep = "\t"))<br>
>179 mirna = as.matrix(read.table("`mirna_list.txt`", head = TRUE, sep = "\t"))<br>
>180 CWCS_matrix1 = as.matrix(read.table("`wMRE_all.txt`", head = TRUE, sep = "\t"))<br>

Those three files serve as the basis to define the sequence matching scores between miRNAs and mRNAs, which are compiled from TargetScan and provided [here](https://github.com/PANWANG2014/miRTIGO/blob/master/miRTIGO_algorithm/Sequence%20matching%20scores.7z). 
	
>188 name\_cancer = as.matrix(read.table("`Sample-to-sample.txt`", head = TRUE, sep = "\t"))<br>
>189 mirna\_cancer = as.matrix(read.table("`Input miRNA expression.txt`", head = FALSE, sep = "\t"))<br>
>190 mrna\_cancer = as.matrix(read.table("`Input mRNA expression.txt`", head = FALSE, sep = "\t"))<br>

Those three files should be provided by the users, which contains the paired expression profiles of miRNAs and mRNAs of the samples that they are interested in. Note that the input miRNA/mRNA file should be transformed into a non-negative matrix, in order for the main program to execute correctly.  

---
### <a name="6">4. Reproduction of the miRTIGO paper`s experiments</a><br>
1. The codes to reproduce these experiments in this paper are written in R and should be executed in the corresponding software environment.<br> 
2. Generally, all these [codes](https://github.com/PANWANG2014/miRTIGO/tree/master/Source_codes_for_experiments) are arranged into three parts as 'FUNCTIONS', 'INPUT DATA' and 'MAIN PROGRAM'. The users need to download and fill in the relevant input files before implementing corresponding analyses.<br>
3. Files required for the reproduction of the experiments can be broadly classified into three categories:<br>

* Files for executing the algorithms<br>

	| Data file | Descriptions | 
	|:-------------:|:-------------| 
	| wMRE\_all.txt | _cumulative weighted context++ scores_ (weighted miRNA response elements, wMREs) between an miRNA and mRNA |
	| qMRE\_all.txt | number of target sites (quantitative miRNA response elements, qMREs) on one mRNA for one miRNA |
	| conserved\_qMRE.txt | number of conserved target sites(qMREs) on one mRNA for one miRNA|
	| mirna\_list.txt | miRNA identifiers used in the above files|
	| mrna\_list.txt | mRNA identifiers used in the above files|

	All these files are compiled from [TargetScan v7.0](http://www.targetscan.org/cgi-bin/targetscan/data_download.cgi?db=vert_70).

* Files for evaluation analyses<br>
	* Experimentally confirmed MMIs<br>

	| Data file | Descriptions | MMI counts |
	|:-------------:|:-------------|:-----:|
	| _V1_ | Tarbase v7.0 | 307,010 |
    | _V2_ | miRTarbase v7.0 | 380,639 |
	| _V3_ | compiled from starBase v2.0 | 26,009 |
	| _V4_ | strong evidence-supported | 9,642 |
	| _V5_ | miRNA transfection-supported | 22,325 |
	| _V6_ | CLASH-supported | 17,293 |


	* Curated cancer-related miRNAs and genes<br>

	| Data file | Descriptions | Molecule counts |
	|:-------------:|:-------------|:-----:|
	| _Oncomirs_ | compiled from MNDR v2.0 database | 399 |
	| _miRNA biomarkers_ | miRNAs that are significantly correlated with tumor development, tumor staging, tumor grade and patient survival | 288 |
	| _Cancer genes_ | COSMIC database | 616 |

* Input data files<br>

	| Data file| Descriptions |
	|:-------------: |:-------------|
	| TCGA data | Paired miRNA-mRNA expression profiles from a total of 7,991 tumor samples belonging to 32 different cancer types |
	| NCI60 data | Paired miRNA-mRNA expression profiles from 130 samples belonging to 59 cell line types |

	All data files listed above are provided in a compressed file [DATA.7z](https://www.dropbox.com/sh/aa0k59j39nftmo9/AAALFIiSpicrAEn8nRUJRjUWa?dl=0). Detailed descriptions of these files are provided in the **Online Methods** of the paper.<br>

---
### <a name="7">5. Examples of predicted results</a><br>
miRTIGO provides both population-level result and sample-level result, and it outputs the top 1% ranked MMIs in default. To show as an example, [Here](https://github.com/PANWANG2014/miRTIGO/blob/master/Predicted_results.7z) we provide the detailed population-level results of the top 5000 predicted MMIs for each of the 32 TCGA cancer types and the sample-level results of the top 5000 predictions for the NCI-60 cell panel.
