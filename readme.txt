#########################################################################################
Software accompaniment to:

"Two-way Mixed-Effects Methods for Joint Association Analyses Using Both Host and Pathogen Genomes. 
M. Wang, F. Roux, C. Bartoli, C. H.-Chauveau, C. Meyer, H. Lee, D. Roby, M. S. McPeek, and J. Bergelson Proc. Natl. Acad. Sci. (direct submission). Vol. 115 (24), E5440-E5449, 2018."

The current version is still undergoing alpha testing, so it may not be able to handle input errors. Please read the demo.m file carefully before running your own analysis. 
We will release a beta version soon. Keep tuned!

#########################################################################################
Installation:
(1) Download the folder with the following Matlab functions:
Kinship_Calculate.m
Objectivefcn_ATOMM.m
ATOMM_NULL.m
ATOMM_Marginal_host.m
ATOMM_Marginal_pathogen.m
ATOMM_Interaction.m

(2)Open Matlab Console in your local computer. Choose Home -> Environment -- > Set Path -- > Add Folder ...
Find the directory of your downloaded folder (with Matlab functions), then choose Open -- > Save -- > Close

(3)Run demo.m code for an example. 
We also provide the example input files and output files in the folder.  

(4)If you want to run the analyses on your own data, please prepare the input files according to the following instruction. Then run the functions as we did in the demo.m.  
(You may have to change the file names and arguments if necessary). See demo.m for details. 

#########################################################################################

Input:

##########################################################################################
(a) host genotype: (the default filename is "sequence_host.txt")

1  1  0  1  1 ...
1  2  1  0  1 ...
1  3  1  0  0 ...
1  4  0  0  0 ...
..............
(1)(2)(3)(4)(5)

Each row of the genotype file specifies the genotype data at one host SNP. The file has in total N_h+2 columns, where N_h is the total number of host individuals. 

(1) the first column is the host chromosome ID.
(2) the second column is the host SNP ID. We require the SNPs to be indexed consecutively in integers from 1 to N, where N is the total number of host SNPs to be tested.
(3) the third column is the SNP profile for host individual 1. We currently allow only haploid individuals, so each SNP is encoded by either 0 or 1.
(4) the fourth column is the SNP profile for host individual 2. Same requirements as in the column 3.
(5) ... the remaining columns specify the SNP profiles for host individuals 3,...,N_h.

##########################################################################################
(b) pathogen genotype: (the default file name is "sequence_pathogen.txt")

1  1  0  0  1  ...
1  2  0  1  1  ...
1  3  0  0  1  ...
1  4  0  0  0  ...
..............
(1)(2)(3)(4)(5)

Each row of the genotype file specifies the genotype data at one pathogen SNP. The file has in total N_p+2 columns, where N_p is the total number of pathogen individuals.

(1) the first column is the pathogen chromosome ID.
(2) the second column is the pathogen SNP ID. We require the SNPs to be indexed consecutively in integers from 1 to N, where N is the total number of pathogen SNPs to be tested.
(3) the third column is the SNP profile for pathogen individual 1. We currently allow only haploid individuals, so each SNP is encoded by either 0 or 1.
(4) the fourth column is the SNP profile for pathogen individual 2. Same requirements as in the column 3.
(5) ... the remaining columns specify the SNP profiles for pathogen individuals 3,...,N_p.

##########################################################################################
(c) phenotype file from cross-factor design (the default file name is "phenotype.txt")

1  1  1 -1.8 ... 0.827
2  1  1 -1.6 ... -0.912
3  1  1  0.8 ... -0.315
...................
(1) (2) (3) (4) .. (last column)

Each row of the phenotype file specifies the observation for one host-pathogen pair. The file has in total 3+p columns, where p (p>=1) is the total number of covariates (including intercept). 
(1) the first column is the host ID, ranging from 1 to N_h.
(2) the second column is the pathogen ID, ranging from 1 to N_p.
(3) the third column is the intercept, where each entry is 1. Note that we always require the intercept to be included in the phenotype file.
(4) the columns 4 -- (2+p) are options. If the study includes additional covariates, please include them as these columns in the file. In the example input file we provided in "phenotype.txt", the study has two covariates (columns 3 and column 4).
(last column) the last column of the phenotype file is the phenotype measurement of interest. 

#########################################################################################
Output:

#########################################################################################
(a) heritability estimates (the default file name is "estimate.txt")

0.1324	0.1334	0.3566	0.3777
(1)	(2)	(3)	(4)

The heritability file contains 4 nonnegative numbers. Each of these four numbers should be between 0 and 1, and they sum up to 1. 

(1)The first number is the proportion of the phenotypical residual variance explained by host genome.
(2)The second number is the proportion of the phenotypical residual variance explained by pathogen genome.
(3)The third number is the proportion of the phenotypical residual variance explained by the interaction between host and pathogen genomes. 
(4)The fourth number is the proportion of the phenotypical residual variance explained by the i.i.d. noise. 

##########################################################################################

(b) marginal association tests in the host genome (the default file name is "marginal_host.txt")

1	1	0.018	0.894614049661
1	2	0.916	0.338634397151
1	3	0.001	0.978073667399
..............
(1) 	(2)	 (3) 	(4) 

Each row of the "marginal_host.txt" file corresponds to one host SNP. 

(1) the first column is the host chromosome ID.
(2) the second column is the host SNP ID. The first two columns of this output file is extracted from the input file "sequence_host.txt". One could specify either a subset of host SNPs, or all host SNPs to be tested. (See options in the ATOMM_marginal_host.m function and demo.m) 
(3) test statistic in the marginal test, which follows chisq with df 1 under the null.
(4) p-value for the marginal association test. 


##########################################################################################
(c) marginal association tests in the pathogen genome (the default file name is "marginal_pathogen.txt")

2	550	0.814	0.366801427723
2	551	0.150	0.698508936922
2	552	1.084	0.297754730683
..............
(1) 	(2)	 (3) 	(4) 

Each row of the "marginal_pathogen.txt" file corresponds to one pathogen SNP. 

(1) the first column is the pathogen chromosome ID.
(2) the second column is the pathogen SNP ID. The first two columns of this output file is extracted from the input file "sequence_pathogen.txt". One could specify either a subset of host SNPs, or all host SNPs to be tested. (See options in the ATOMM_marginal_pathogen.m function and demo.m) 
(3) test statistic in the marginal test, which follows chisq with df 1 under the null.
(4) p-value for the marginal association test. 

######################################################################################################
(d) interaction association tests between host and pathogen SNPs (the default file name is "interaction.txt")

1	20	1	50	0.124	0.724900402022
1	20	1	51	1.448	0.228810208035
1	20	1	52	0.003	0.957317302221
..............
1	27	1	50	0.230	0.631349297954
1	27	1	51	0.020	0.886994564991
1	27	1	52	0.940	0.332290276995
..............
(1) 	(2)	 (3) 	(4)     (5)     (6)

Each row of the "interaction.txt" file is the interaction test between one host SNP and one pathogen SNP. One could specify which host SNPs and pathogen SNPs to be tested.

<WARNING!> In general, performing interaction tests between all host SNPs and all pathogen SNPs would be time consuming, so we recommend users to perform interaction tests only for a subset of host SNPs and pathogen SNPs. One possibility is to conduct marginal association tests first, and then select significant host/pathogen SNPs for interaction tests. (see demo.m for details). 

(1) the first column is the host chromosome ID.
(2) the second column is the host SNP ID. The SNP IDs should be indexed consecutively by integers, starting from 1. 
(1) the third column is the pathogen chromosome ID.
(2) the fourth column is the pathogen SNP ID. The SNP IDs should be indexed consecutively by integers, starting from 1. 
(5) test statistic in the interaction test, which follows chisq with df 1 under the null.
(6) p-value for the interaction association test. 







