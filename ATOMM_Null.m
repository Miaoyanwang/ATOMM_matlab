%%%%% Input: 
%%%%% Phenotype: Phenotype file
%%%%% Kinship_h: Host empirical kinship matrix
%%%%% Kinship_p: Pathogen empirical kinship matrix

%%%%% Description: this function performs estimation for two-may mixed effects model. There are four variance components in the model: random host effects, random pathogen effects, random host-pathogen interaction effects, and i.i.d. noise. 


function[herit,Sigma]=ATOMM_Null(Kinship_h,Kinship_p,Phenotype)
estID=fopen('output/estimate.txt','w');

nY=size(Phenotype);
n=nY(1);

Kinship_htotal=Kinship_h(Phenotype(:,1),Phenotype(:,1));
Kinship_ptotal=Kinship_p(Phenotype(:,2),Phenotype(:,2));
Kinship_cross=Kinship_htotal.*Kinship_ptotal;


X=Phenotype(:,3:(nY(2)-1));
Y=Phenotype(:,nY(2));

x0=[0.25,0.25,0.25];
fun=@(x)objectivefcn_ATOMM(x,Kinship_htotal,Kinship_ptotal,Kinship_cross,X,Y,n);
options = optimset('Display','iter','PlotFcns',@optimplotfval);
optnew = optimset(options,'TolX',0.1);
optnew1 = optimset(optnew,'TolFun',0.1);
[x,fval,exitflag,output] = fminsearch(fun,x0,optnew1);
x0=x;
[x,fval,exitflag,output] = fminsearch(fun,x0,optnew1);

herit=x/(sum(x)+1);
Sigma=herit(1)*Kinship_htotal+herit(2)*Kinship_ptotal+herit(3)*Kinship_cross+(1-sum(herit))*eye(n);


fprintf(estID,'%.4f\t%.4f\t%.4f\t%.4f\n',herit(1),herit(2),herit(3),1-sum(herit));
fclose(estID);

end
