g=read.table("sequence_host.txt")
g=as.matrix(g[,-c(1,2)])
mean=apply(g,1,mean)
center=(g-mean)/sqrt(mean*(1-mean))
K1=t(center)%*%center/dim(g)[1]

g=read.table("sequence_pathogen.txt")
g=as.matrix(g[,-c(1,2)])
mean=apply(g,1,mean)
center=(g-mean)/sqrt(mean*(1-mean))
K2=t(center)%*%center/dim(g)[1]

n1=dim(K1)[1]
n2=dim(K2)[1]
cross=cbind(rep(1:n1,n2),rep(1:n2,rep(n1,n2)))

K1=K1[cross[,1],cross[,1]]
K2=K2[cross[,2],cross[,2]]
K3=K1*K2
n=dim(K3)[1]

X=mvrnorm(1,mu=rep(0,n),Sigma=diag(1,n))
X=cbind(1,X)
Y=mvrnorm(1,mu=X%*%c(0.1,0.2),Sigma=0.1*K1+0.2*K2+0.3*K3+0.4*diag(1,n))
write.table(cbind(cross,X,Y),"phenotype.txt",row.names=F,col.names=F)
