%%%%% Input: 
%%%%% Phenotype: Phenotype file
%%%%% Sigma: Phenotypic covariance matrix with three components 
%%%%% Gh: Pathogen genotype file
%%%%% Kinship_h: Pathogen empirical kinship matrix
%%%%% index: which pathogen SNPs to be tested

%%%%% Description: this function performs retrospective tests for marginal association of pathogen SNPs.


function[result]=ATOMM_Marginal_pathogen(Phenotype,Sigma,Gh,Kinship_h,index)

fileID = fopen('output/marginal_pathogen.txt','w');

nY=size(Phenotype);
n=nY(1);
nh=size(Gh);
nh=nh(2)-2;
vector_h=ones(1,nh);

X=Phenotype(:,3:(nY(2)-1));
Y=Phenotype(:,nY(2));
p=nY(2)-3;
beta=eye(p)/(X'*(Sigma\X))*(X'*(Sigma\Y));
mu=X*beta;


temG=zeros(n,nh);
for i =1:nh
temG(find(Phenotype(:,2)==i),i)=1;% second column of phenotype file is the pathogen ID
end

deno=(Y-mu)'*(Sigma\temG)*Kinship_h*(Sigma\temG)'*(Y-mu);
temY=(Y-mu)'*(Sigma\temG);
matrix_inv=pinv(Kinship_h);
result=zeros(length(index),3);

for i = 1:length(index)
g=Gh(index(i),3:(nh+2))';
%%sigma_g=g'*matrix_inv*g-(g'*matrix_inv*vector_h')^2/(vector_h*matrix_inv*vector_h');
sigma_g=g'*matrix_inv*g;
sigma_g=sigma_g/(nh-1);   
result(i,1:2)=Gh(index(i),1:2);
result(i,3)=(temY*g)^2/(deno*sigma_g);
result(i,4)=1-chi2cdf(result(i,3),1);
fprintf(fileID,'%d\t%d\t%.3f\t%.12f\n',result(i,1),result(i,2),result(i,3),result(i,4));

end

end
