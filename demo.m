Gh=dlmread('input/sequence_host.txt'); %% read in host genotype file.
Gp=dlmread('input/sequence_pathogen.txt');%% read in pathogen genotype file.
[Kinship_h,MAF_h]=Kinship_Calculate(Gh); %% calculate host genetic relatedness matrix of size N_h * N_h. The (i,j) entry in the matrix Kinship_h is the pair-wise empirical genetic relatedness between hosts i and j.
[Kinship_p,MAF_p]=Kinship_Calculate(Gp);%% calculate pathogen genetic relatedness matrix of size N_p * N_p. The (i,j) entry in the matrix Kinship_p is the pair-wise empirical genetic relatedness between pathogens i and j.

Phenotype=dlmread('input/phenotype.txt'); %% read in phentype file.

[herit,Sigma]=ATOMM_Null(Kinship_h,Kinship_p,Phenotype);%% estimate heritability under the null model; Results are saved in the "output/estimate.txt" file. In the matlab console, an window named Optimization PlotFcs will pop out. One should expect to observe the decreasement of the objective function. 

index=1:100; %% specify which host SNP to be included in the marginal association analysis. In this example, we only test the first 100 host SNPs from the "input/sequence_host.txt" file. If you want to include all host SNPs in the analysis, please use "index=1:N", where N is the total number of SNPs in the host genotype file.   

stat_host=ATOMM_Marginal_host(Phenotype,Sigma,Gh,Kinship_h,index); %% run marginal association tests for the specified host SNPs. Results are saved in the "stat_host" and "output/marginal_host.txt" file.

index=550:600; %%specify which pathogen SNPs to be included in the marginal association analysis. In this example, we only test the 550th -- 600th SNPs (rows) from the "input/sequence_pathogen.txt" file. If you want to include all pathogen SNPs in the analysis, please use "index=1:N", where N is the total number of SNPs in the pathogen genotype file.  

stat_pathogen=ATOMM_Marginal_pathogen(Phenotype,Sigma,Gp,Kinship_p,index); %% %% run marginal association tests for the specified pathogen SNPs. Results are saved in the "stat_pathogen" and "output/marginal_pathogen.txt" file.

index_h=20:30; %%specify which host SNPs to be included in the interaction tests. In this example, we only test the 20th - 30th SNPs from the "sequence_host.txt" file. In general, we would recommend users to run marginal tests first, and select a small subset (say, tens) of host SNPs with significant p-values. Then identify the IDs for these SNPS (i.e. the row number in the "input/sequence_host.txt" file) and specify "index_p" using these IDs. 

index_p=50:60; %%specify which host SNPs to be included in the interaction tests. In this example, we only test the 50th - 60th SNPs from the "sequence_pathogen.txt" file. In general, we would recommend users to run marginal tests first, and then select a small subset (say, tens) of pathogen SNPs with significant p-values. Then identify the IDs for these SNPS (i.e. the row number in the "input/sequence_pathogen.txt" file) and specify "index_h" using these IDs.

stat_interaction=ATOMM_Interaction(Phenotype,Sigma,Gh,Gp,Kinship_h,Kinship_p,index_h,index_p); %% run the interaction tests between the prespecified host SNPs and pathogen SNPs. The results are saved in "stat_interaction" and file "output/interaction.txt".
