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
(2) the second column is the host SNP ID. We require the SNPs to be indexed consecutively in integer from 1 to N, where N is the total number of host SNPs to be tested.
(3) the third column is the SNP profile for host individual 1. We currently allow only haploid individuals, so each SNP is encoded by either 0 or 1.
(4) the fourth column is the SNP profile for host individual 2. Same requirements as the column 3.
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
(2) the second column is the pathogen SNP ID. We require the SNPs to be indexed consecutively in integer from 1 to N, where N is the total number of pathogen SNPs to be tested.
(3) the third column is the SNP profile for pathogen individual 1. We currently allow only haploid individuals, so each SNP is encoded by either 0 or 1.
(4) the fourth column is the SNP profile for pathogen individual 2. Same requirements as the column 3.
(5) ... the remaining columns specify the SNP profiles for pathogen individuals 3,...,N_p.

##########################################################################################
(c) phenotype file from cross-factor design (the default file name is "phenotype.txt")

1  1  1  0.367679175686317 .. 1.10935540806447
2  1  1  -1.55276219195733 .. -0.776654813319241
3  1  1  0.598996909239343 .. 0.513856783747087
..............
(1) (2) (3) (4) .. (last column)

Each row of the phenotype file specifies the observation for one host-pathogen pair. The file has in total 3+p columns, where p (p>=1) is the total number of covariates (including intercept). 
(1) the first column is the host ID, ranging from 1 to N_h.
(2) the second column is the pathogen ID, ranging from 1 to N_p.
(3) the third column is the intercept, where each entry is 1. Note that we always require the intercept to be included in the phenotype file.
(4) the columns 4 -- (2+p) are options. If the study concludes dditional covariates, please include them as columns in the file. In the example phenotype file "phenotype.txt", we include two covariates. 
(last column) the last column of the phenotype file is the phenotype measurement of interest. 


Output:

#########################################################################################
(a) marginal association tests in the host genome (the default file name is "marginal_host.txt")

1	1	0.660	0.416
1	2	0.452	0.502
1	3	3.527	0.060
..............
(1) 	(2)	 (3) 	(4) 

Each row of the "marginal_host.txt" file corresponds to one host SNP. 
(1) the first column is the host chromosome ID.
(2) the second column is the host SNP ID. The first two columns of this output file is extracted from the input file "sequence_host.txt". One could specify either a subset of host SNPs, or all host SNPs to be tested. (See options in the ATOMM_marginal_host.m function and demo.m) 
(3) test statistic in the marginal test, which follows chisq with df 1 under the null.
(4) p-value for the marginal association test. 

######################################################
(b) marginal association tests in the pathogen genome (the default file name is "marginal_pathogen.txt")

1	1	0.504	0.478
1	2	0.005	0.946
1	3	2.569	0.109
..............
(1) 	(2)	 (3) 	(4) 

Each row of the "marginal_pathogen.txt" file corresponds to one pathogen SNP. 

(1) the first column is the pathogen chromosome ID.
(2) the second column is the pathogen SNP ID. The first two columns of this output file is extracted from the input file "sequence_pathogen.txt". One could specify either a subset of host SNPs, or all host SNPs to be tested. (See options in the ATOMM_marginal_pathogen.m function and demo.m) 
(3) test statistic in the marginal test, which follows chisq with df 1 under the null.
(4) p-value for the marginal association test. 

######################################################
(c) interaction association tests between host and pathogen SNPs (the default file name is "interaction.txt")

1	1	1	1	0.149	0.699
1	1	1	2	0.149	0.700
1	1	1	3	0.150	0.699
..............
1	3	1	1	0.001	0.979
1	3	1	2	0.001	0.981
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







