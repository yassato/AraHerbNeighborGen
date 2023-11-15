###########################################
# numerical examples for phenotype values #
###########################################

library(tidyverse)
library(patchwork)

plt_f = function(b0,b1,b2,b12,pch) {
  f_star = 0.5-(b1/(2*b2))
  p2 = ggplot(NULL, aes(x=f_star,y=(b12+b2)*(2*f_star-1)+b0+b1)) + geom_point(pch=pch,size=3) +
    theme_classic() + ylab("Phenotypic value") + xlab("Frequency of AA genotypes") + xlim(0,1) + ylim(0,1) +
    geom_function(aes(x=1,y=1),fun=function(x,b0,b1,b2,b12) { (b12+b2)*(2*x-1)+b0+b1 }, args=list(b0,b1,b2,b12)) + 
    geom_function(aes(x=1,y=1),fun=function(x,b0,b1,b2,b12) { (b12-b2)*(2*x-1)+b0-b1 }, args=list(b0,b1,b2,b12), colour=grey(0.0,0.33)) +
    geom_function(aes(x=1,y=1),fun=function(x,b0,b1,b2,b12) { 2*b2*x*(2*x-1)+(b12-b2)*(2*x-1)+b0-b1+2*b1*x }, args=list(b0,b1,b2,b12), lwd=1.5,col=rgb(0,0,1,1))
  return(p2)
}

p1 = plt_f(b0 = 0.5, b1 = 0.0, b2 = 0.2, b12 = 0.0, pch=NA)
p1 = p1 + geom_text(aes(x=0.025,y=1.0),label="--- (black): AA genotypes",size=3,hjust=0) + 
  geom_text(aes(x=0.025,y=0.9),label="--- (gray): aa genotypes",colour=grey(0.5,1),size=3,hjust=0) +
  geom_text(aes(x=0.025,y=0.8),label="---  (blue): weighted mean",size=3,hjust=0,colour=rgb(0,0,1,1)) + 
  labs(subtitle=bquote(beta[0]*" = 0.5; "*beta[1]*" = 0.0; "*beta[2]*" = "*.(format(0.2,nsmall=1))))

p2 = plt_f(b0 = 0.5, b1 = 0.0, b2 = -0.2, b12 = 0.0, pch=NA)
p2 = p2 + labs(subtitle=bquote(beta[0]*" = 0.5; "*beta[1]*" = 0.0; "*beta[2]*" = "*.(format(-0.2,nsmall=1))))

p3 = plt_f(b0 = 0.5, b1 = 0.05, b2 = 0.2, b12 = 0.0, pch=NA)
p3 = p3 + labs(subtitle=bquote(beta[0]*" = 0.5; "*beta[1]*" = 0.05; "*beta[2]*" = "*.(format(0.2,nsmall=1))))

p4 = plt_f(b0 = 0.5, b1 = 0.05, b2 = -0.2, b12 = 0.0, pch=NA)
p4 = p4 + labs(subtitle=bquote(beta[0]*" = 0.5; "*beta[1]*" = 0.05; "*beta[2]*" = "*.(format(-0.2,nsmall=1))))

p = (p1 | p2) / (p3 | p4)

ggsave(p1,filename="../figs/FDSinbred.pdf",width=2.5,height=2.5)
# ggsave(p,filename="../figs/FDSinbredSupp.pdf",width=6,height=6)



#define functions
XiXj1 = function(id,dmat,range) {
  xi = Xi[id]
  if(xi==0) { xi = 1 } else { xi = xi }
  j = as.numeric(which((dmat[,id]<=range)&(dmat[,id]>0)))
  xj = Xi[j]
  xj[xj==0] = 1
  return(xi*sum(xj))
}

#set coefficients
J = 0.2 #= beta_2
h = 0.0001 #= beta_1

#set a field
rect = 100
N = rect*rect
Xi = sample(c(-1,1),N,replace=T)

smap = cbind(rep(1:rect,each=rect),
             rep(1:rect,rect)
)
dmat = as.matrix(dist(smap))

group = rep(1,rect*rect)


xixj = c()
for(i in 1:N) xixj = c(xixj, XiXj1(i,dmat=dmat,range=sqrt(2)))
Ei_t0 = (J*xixj+h*Xi)

#MCMC by Gibbs sampling
for(j in 1:100) {
  Xi_t0 = Xi
  
  perm = sample(1:N,N)
  for(i in perm) {
    Xi[i] = sample(c(-1,1),1)
    xi = Xi[i]
    if(xi==0) { xi = 1 }
    Ei = -(J*XiXj1(i,dmat=dmat,range=sqrt(8))+h*xi) #set the energy negative to minimize itself
    
    #Metropolis algorithm
    if(Ei_t0[i]<=Ei) {
      Xi[i] = Xi[i]; Ei_t0[i] = Ei
    } else if(runif(1,0,1)<exp(Ei-Ei_t0[i])) {
      Xi[i] = Xi[i]; Ei_t0[i] = Ei
    } else {
      Xi[i] = Xi_t0[i]
    }
  }
  E = mean(Ei_t0)
  image(matrix(Xi,rect,rect),main=paste(j,E),col=c("black","grey","white"))
}

d1 = expand.grid(x=c(1:rect),y=c(1:rect))
d1$z = Xi
ip1 = ggplot(d1,aes(x=x,y=y,fill=z)) + geom_tile(color="grey50") + 
  scale_fill_gradient(low="white", high="black") + 
  labs(subtitle=bquote(beta[0]*" = 0.0; "*beta[1]*" = 0.00; "*beta[2]*" = "*.(format(0.2,nsmall=1)))) +
  theme_classic() + ylab("") + xlab("") + theme(legend.position="none")

d2 = expand.grid(x=c(1:rect),y=c(1:rect))
d2$z = Xi
ip2 = ggplot(d2,aes(x=x,y=y,fill=z)) + geom_tile(color="grey50") + 
  scale_fill_gradient(low="white", high="black") + 
  labs(subtitle=bquote(beta[0]*" = 0.0; "*beta[1]*" = 0.00; "*beta[2]*" = "*.(format(-0.2,nsmall=1)))) +
  theme_classic() + ylab("") + xlab("") + theme(legend.position="none")

ip = (ip1 | ip2) + plot_annotation(tag_levels = "a")
ipp = (ip1 | ip2) / (p1 | p2) / (p3 | p4) + plot_annotation(tag_levels = "a")

# Figure S1
ggsave(ipp,filename="../figs/IsingFDS.pdf",width=6,height=8)


