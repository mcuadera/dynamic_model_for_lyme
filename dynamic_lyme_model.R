#Modeling interactions of host abundance and nymphal density
#Based on the paper: 'Deer, predators, and the emergence of Lyme disease' by Levi et al.
#Dynamic model interpretation in R by Mervin Keith Cuadera
#03/12/2020

rm(list=ls(all=TRUE))
library(deSolve)

#Explanation of parameters:

#N_m+F_alt is the total population of hosts
#r is maximum intrinsic growth rate
#K is carrying capacity (individuals)
#aP is maximum predation rate (kills per km^2)  
#c is half saturation parameter (this is just half of aP because this is the same as max_rate/2)
#T_mt is probability of being biten by an infected nymph
#F_alt fraction of ticks biting alternative hosts, but also increase total host abundance (individuals)
#beta is the tick bite rate; (originally written as beta(Nm+F))
#b0 is half-saturation parameter representing the density of small mammals where half of ticks
#would be expected to feed. (half of ticksexpected to feed when total population is km^-2) 
#v is the birth rate (they vary this as either 1,1.5, 0.5 million larva born per km^2 annually)
#mu_L is the per-capita death rate of larval ticks
#mu_N is the per-capita death rate of nymphal ticks
#T_tm is the is successful infection of Borrelia from an infected hosts.


#Explanation of classes:

#N_m is functional group of small mammal host density
#S_m is susceptible small-mammal host population
#I_m is infected small-mammal host population
#I_t is the infected nymphs
#S_t is the larval ticks, which are all susceptible
#J_t is the uninfected nymphs

host_vector_model <- function(t,y,parms) {
  #define state variables
  N_m <- y[1]
  S_m <- y[2]
  I_m <- y[3]
  S_t <- y[4]
  I_t <- y[5]
  J_t <- y[6]
  #defining parameters used for the differential equation
  r     <- parms[1]
  K     <- parms[2]
  aP    <- parms[3] #can be 1,000-9,000
  c     <- parms[4]
  T_mt  <- parms[5]
  F_alt <- parms[6]
  beta  <- parms[7]
  b0   <- parms[8]
  v     <- parms[9] #can be 500,000 1,000,000 1,500,000
  mu_L  <- parms[10]
  mu_N  <- parms[11]
  T_tm  <- parms[12]

  #differential equations of state variables
    #order:
      #N_m
      #S_m
      #I_m
      #S_t
      #I_t
      #J_t
  dN = numeric(6)
  dN[1] = r*N_m*(1-N_m/K)-(aP*N_m^2)/(c^2+N_m^2)  #N_m 
  
  dN[2] = r*N_m*(1-N_m/K)-(T_mt*I_t*S_m)/(b0+N_m+F_alt)-S_m*((aP*N_m)/(c^2+N_m^2))  #S_m
  
  dN[3] = (T_mt*I_t*S_m)/(b0+N_m+F_alt)-I_m*(aP*N_m)/(c^2+N_m^2) #I_m 
  
  dN[4] = v-(N_m+F_alt)/(b0+N_m+F_alt)*S_t-mu_L*S_t #S_t
  
  dN[5] = (T_tm*I_m*S_t)/(b0+N_m+F_alt)-(N_m+F_alt)/(b0+N_m+F_alt)*I_t-mu_N*I_t #I_t
  
  dN[6] = ((S_m+F_alt+I_m*(1-T_tm))/(b0+N_m+F_alt))*S_t-((N_m+F_alt)/(b0+N_m+F_alt)*J_t)-mu_N*J_t #J_t
  
  jbad = which(  (dN<0)&(y<=0) );
  dN[jbad]=0; #this removes biologically impossible situations like negative population.
  
  return(list(dN))
}

r <- 2;K <- 10000;aP_0 <- 1000;c <- 2500;T_mt <- 0.9;F_alt <- 4120;beta <- 0.10;
b0 <- 80000;v <- 1500000; mu_L <- 0.2;mu_N <- 0.2;T_tm <- 0.9
#aP was varied from 1000 to 9000 on the paper. Thus I used aP_0 to define initial aP
#I will loop this about 80,000 times, where I define aP inside the loop function instead.

infected_nymphs_1_5 <- rep(0,1000) #using v = 1.5 million
infected_nymphs_1 <- rep(0,1000) #using v = 1 million
infected_nymphs_500 <- rep(0,1000) #using v = 500 million

#generating aP values (max predation rate)
aP_val <- rep(0,1000); aP_val[1] = aP_0;
for (i in 2:1000){
  aP_val[i] <- 10 + aP_val[i-1]
}

#this is with v = 1.5 million (tick birth rate)
for (i in 2:1000) {
  v <- 1500000
  aP <- aP_val[i-1];
  parms <- c(r,K,aP,c,T_mt,F_alt,beta,b0,v,mu_L,mu_N,T_tm);
  out<-ode(y=c(1,1,1,1,1,1), times= seq(0,200,by=.1), func=host_vector_model, parms)
  infected_nymphs_1_5[i-1] <- tail(out[1000,6])
}
#infected_nymphs_1_5[1] results in 1137883

#this is with v = 1 million
for (i in 2:1000) {
  v <- 1000000
  aP <- aP_val[i-1];
  parms <- c(r,K,aP,c,T_mt,F_alt,beta,b0,v,mu_L,mu_N,T_tm);
  out<-ode(y=c(1,1,1,1,1,1), times= seq(0,200,by=.1), func=host_vector_model, parms)
  infected_nymphs_1[i-1] <- tail(out[1000,6])
}

#this is with v = 0.5 million
for (i in 2:1000) {
  v <- 500000
  aP <- aP_val[i-1];
  parms <- c(r,K,aP,c,T_mt,F_alt,beta,b0,v,mu_L,mu_N,T_tm);
  out<-ode(y=c(1,1,1,1,1,1), times= seq(0,200,by=.1), func=host_vector_model, parms)
  infected_nymphs_500[i-1] <- tail(out[1000,6])
}

#plotting the results
plot(aP_val,infected_nymphs_1_5, type='l',lty=3,xlim = c(1000,9000), ylim = c(0,12*10^5),
     xlab = 'Asymptotic Predation Rate', ylab = 'Infected Nymphs (km^-2)')
points(aP_val,infected_nymphs_1, type = 'l',lty=1)
points(aP_val,infected_nymphs_500,type = 'l', lty=2)
legend("topright", legend = c("v = 1.5 million", 
                             "v = 1 million",
                             "v = 0.5 million"), lty = c(3,1,2), cex = 0.7)

#This is the code used inside the loops. Use for checking results.
  #aP <- 1000
  #parms <- c(r,K,aP,c,T_mt,F_alt,beta,b0,v,mu_L,mu_N,T_tm)  
  #out<-ode(y=c(1,1,1,1,1,1), times= seq(0,200,by=.1), func=host_vector_model, parms)
  #head(out)
  #tail(out[1000,6])
  #infected_nymphs_1_5[1] results in 1137883 confirming our loop works properly.

#Looking at how varying F affects nymphs at a constant birth rate and different predation rates. 
infected_nymphs_l <- rep(0,1000) #using aP = 4000
infected_nymphs_m <- rep(0,1000) #using aP = 6000
infected_nymphs_h <- rep(0,1000) #using v = 8000

#generating F values
F_0 <- 1000
F_val <- rep(0,1000); F_val[1] = F_0;
for (i in 2:1000){
  F_val[i] <- 10 + F_val[i-1]
}

#at low predation
for (i in 2:1000) {
  aP <- 4000
  F_alt <- F_val[i-1];
  parms <- c(r,K,aP,c,T_mt,F_alt,beta,b0,v,mu_L,mu_N,T_tm);
  out<-ode(y=c(1,1,1,1,1,1), times= seq(0,200,by=.1), func=host_vector_model, parms)
  infected_nymphs_l[i-1] <- tail(out[1000,6])
}

#at medium predation
for (i in 2:1000) {
  aP <- 6000
  F_alt <- F_val[i-1];
  parms <- c(r,K,aP,c,T_mt,F_alt,beta,b0,v,mu_L,mu_N,T_tm);
  out<-ode(y=c(1,1,1,1,1,1), times= seq(0,200,by=.1), func=host_vector_model, parms)
  infected_nymphs_m[i-1] <- tail(out[1000,6])
}

#at high predation
for (i in 2:1000) {
  aP <- 8000
  F_alt <- F_val[i-1];
  parms <- c(r,K,aP,c,T_mt,F_alt,beta,b0,v,mu_L,mu_N,T_tm);
  out<-ode(y=c(1,1,1,1,1,1), times= seq(0,200,by=.1), func=host_vector_model, parms)
  infected_nymphs_h[i-1] <- tail(out[1000,6])
}

#plotting the results
plot(F_val,infected_nymphs_l, type='l',lty=3,xlim = c(1000,9000),
     xlab = 'Density of Dilution Hosts', ylab = 'Infected Nymphs (km^-2)')
points(F_val,infected_nymphs_m, type = 'l',lty=1)
points(F_val,infected_nymphs_h,type = 'l', lty=2)
legend("topright", legend = c("aP = 4000", 
                              "aP = 6000",
                              "aP = 8000"), lty = c(3,1,2), cex = 0.7)
