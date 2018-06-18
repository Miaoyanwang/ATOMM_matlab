%%%%% Input: 
%%%%% Phenotype: Phenotype file
%%%%% Sigma: Phenotypic covariance matrix with three components 
%%%%% Gh: Host genotype file
%%%%% Gp: Pathogen genotype file
%%%%% Kinship_h: Host empirical kinship matrix
%%%%% Kinship_p: Pathogen empirical kinship matrix
%%%%% index_h: which host SNPs to be tested
%%%%% index_p: which pathogen SNPs to be tested

%%%%% Description: this function performs retrospective tests for interaction between specified host SNPs and pathogen SNPs. 

function[result]=ATOMM_Interaction(Phenotype,Sigma,Gh,Gp,Kinship_h,Kinship_p,index_h,index_p);

fileID=fopen('output/interaction.txt','w');

nY=size(Phenotype);
n=nY(1);
nh=size(Gh);
nh=nh(2)-2;
np=size(Gp);
np=np(2)-2;
p=nY(2)-3;

X=Phenotype(:,3:(nY(2)-1));
Y=Phenotype(:,nY(2));
vector=ones(1,n);

Kinship_htotal=Kinship_h(Phenotype(:,1),Phenotype(:,1));
Kinship_ptotal=Kinship_p(Phenotype(:,2),Phenotype(:,2));
Kinship_cross=Kinship_htotal.*Kinship_ptotal;
matrix_inv=pinv(Kinship_cross);
Sigma_inv=inv(Sigma);
result=zeros(length(index_h)*length(index_p),6);
ncount=0;


for i = 1:length(index_h)
for j = 1:length(index_p)
ncount=ncount+1;
gh=Gh(index_h(i),3:(nh+2));
gp=Gp(index_p(j),3:(np+2));
X_new=[X';gh(Phenotype(:,1));gp(Phenotype(:,2))]';
g=X_new(:,(p+1)).*X_new(:,(p+2));
beta=(X_new'*(Sigma\X_new))\(X_new'*(Sigma\Y));
mu=X_new*beta;
       
deno=(Y-mu)'*Kinship_cross*Sigma_inv*Kinship_cross*(Y-mu);
temp_vector=g'*matrix_inv*vector';
%sigma_g=g'*matrix_inv*g-temp_vector'/(vector*matrix_inv*vector')*temp_vector;
sigma_g=g'*matrix_inv*g;
sigma_g=sigma_g/(n-1);                                      
result(ncount,1:2)=Gh(index_h(i),1:2);
result(ncount,3:4)=Gp(index_p(j),1:2);                                      
result(ncount,5)=((Y-mu)'*Sigma_inv*g)^2/(deno*sigma_g);
result(ncount,6)=1-chi2cdf(result(ncount,5),1);
fprintf(fileID,'%d\t%d\t%d\t%d\t%.3f\t%.12f\n',result(ncount,1),result(ncount,2),result(ncount,3),result(ncount,4),result(ncount,5),result(ncount,6));
end
end
end
