%% Input: 
%% para: three coefficients in the variance-components model. 
%% Part1: host kinship 
%% Part2: pathogen kinship 
%% Part3: host-by-pathogen kinship 
%% X: covariates including intercept
%% Y: phenotype of interest
%% n: total number of observations

%% Description: This function specifies the objective function, i.e., negative log-likelihood (up to a constant)

function result=objectivefcn_ATOMM(para,Part1,Part2,Part3,X,Y,n)
thresh=10^(-6);
if(min(para)<thresh)
    result=Inf;
else
Sigma=para(1)*Part1+para(2)*Part2+para(3)*Part3+eye(n);
scale=sum(para)+1;
Sigma=Sigma/scale;
%%Sigma_inv=inv(Sigma); 
p=size(X);
p=p(2);
    
beta=eye(p)/(X'*(Sigma\X))*(X'*(Sigma\Y));
mu=X*beta;
sigma_t=(Y-mu)'*(Sigma\(Y-mu))/n;
 
Temp=chol(Sigma);
%%result=n/2*log(2*pi)+n/2*log(sigma_t)+sum(log(diag(Temp)))+n/2;%% correct minus log-likelihood function

result=2*sum(log(diag(Temp)))+n*log(sigma_t);

end


