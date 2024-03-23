#simulation distribution from infection to report 
#combine two gamma distribution (infection to onset (incubation) & onset to report)
comb.dist=function(mu1, sd1, mu2, sd2, n.g){
  p1=mu1^2/sd1^2; p2=mu1^1/sd1^2
  p3=mu2^2/sd2^2; p4=mu2^1/sd2^2
  n.sim=1e5
  temp=rgamma(n.sim, p1, p2)+rgamma(n.sim, p3, p4)
  x1=c(0,0:(n.g-2)+1.5); x2=c(1.5, 0:(n.g-2)+2.5)
  prob=as.vector(do.call('cbind', lapply(1:n.g, function(i){
    mean(temp>x1[i]&temp<=x2[i])
  })))
  prob=prob/sum(prob)
}

#incubation (SG, omicron): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10045676/#R5
i.mean=3; i.sd=1.55 
# onset to report (source: https://www.mdpi.com/2077-0383/11/12/3269)
n.g=21; g.mean=3.18; g.sd=sqrt(10.33)
g=comb.dist(i.mean, i.sd, g.mean, g.sd, n.g)
#onset to hospitalization (source: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0273134)
n.h=28; h.mean=5;h.sd=4.2
h=comb.dist(i.mean, i.sd, h.mean,h.sd, n.h)
#onset to death (source: https://www.mdpi.com/2077-0383/11/12/3269)
f.mean=20.57; f.sd=sqrt(180.62); n.f=60
f=comb.dist(i.mean, i.sd, f.mean,f.sd, n.f)

#viral shedding distribution (source: https://ehp.niehs.nih.gov/doi/full/10.1289/EHP10050)
v.mean=6; v.sd=7; n.v=56
#v.mean=4.7; v.sd=1.7; n.v=21
v=comb.dist(i.mean, i.sd, v.mean,v.sd, n.v)
  
#serial interval: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10045676/#R5
n.w=14; w <- EpiEstim::discr_si(1:n.w, 2.8, 1.5)

#multiple age group death
#corresponding onset to death distribution (source: https://www.mdpi.com/2077-0383/11/12/3269)
if(n.age>1){
  n.f=60
  if(n.age==2){f.mean=c(20.5, 23); f.sd=c(sqrt(265),sqrt(200))}
  if(n.age==3){f.mean=c(21.3, 22.5, 20); f.sd=sqrt(c(265, 260, 180))}
  f=matrix(0,n.f,n.age)
  for(i in 1:n.age)
  {
    f[,i]=comb.dist(i.mean, i.sd, f.mean[i],f.sd[i], n.f)
  }
}


