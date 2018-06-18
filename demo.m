Gh=dlmread('sequence_host.txt');
Gp=dlmread('sequence_pathogen.txt');
[Kinship_h,MAF_h]=Kinship_Calculate(Gh);
[Kinship_p,MAF_p]=Kinship_Calculate(Gp);

Phenotype=dlmread('phenotype.txt');

[herit,Sigma]=ATOMM_Null(Kinship_h,Kinship_p,Phenotype);

index=1:10000; % the row index to be tested 
stat_host=ATOMM_Marginal_host(Phenotype,Sigma,Gh,Kinship_h,index);

index=1:1000; % the row index to be tested
stat_pathogen=ATOMM_Marginal_pathogen(Phenotype,Sigma,Gp,Kinship_p,index);

index_p=1:10;
index_h=1:1000;
stat_interaction=ATOMM_Interaction(Phenotype,Sigma,Gh,Gp,Kinship_h,Kinship_p,index_h,index_p);
