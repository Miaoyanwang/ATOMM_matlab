function[Kinship,MAF]=Kinship_Calculate(G)
thresh=0.05;
dim_G=size(G);
nh=dim_G(2)-2;
nG=dim_G(1);

matrix=zeros(nh,nG);
MAF=ones(nG,1);


for i = 1:nG
g=G(i,3:(nh+2))';
MAF(i)=mean(g);
matrix(:,i)=(g-MAF(i))/sqrt(MAF(i)*(1-MAF(i)));
end


Kinship=matrix(:,find((MAF>thresh)&(MAF<(1-thresh))))*matrix(:,find((MAF>thresh)&(MAF<(1-thresh))))'/ length(find((MAF>thresh)&(MAF<(1-thresh))));
end
