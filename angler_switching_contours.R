#################################################################
# angler behavior in spatial fishery with two contrasting species
#################################################################

library(scales)
library(tidyr)
library(dplyr)

#
# Scenario 1: CW specialists
#

# one warm water species, one cold water species
# anglers can move anywhere - movement (and therefore fishing pressure) is dependent on preference and abundance

###################
# model parameters
###################

#---------------------
# site characteristics

years = 500        # years
lat = 10           # grid height


init.temps = seq(0, 10, length.out = lat) # temperatures across time and latitude - colder in the north
temps = array(NA, dim = c(lat, years))
temps[,1] <- init.temps

cc = years/5 # total temp increase over time in degrees

for(i in 1:lat) { temps[i,] <- temps[i,1] + seq(0,cc, length.out = years) }


#------------------------
# species characteristics

r = 0.1       # maximum pop growth rate

rho = seq(0, .5, length.out = 100)  # strength of competition

theta = 0.1   # temperature sensitivity
t_opt1 = 10   # optima for species 1
t_opt2 = 20   # optima for species 2

# differences from temperature optima
delta.t_opt1 = abs(temps - t_opt1)
delta.t_opt2 = abs(temps - t_opt2)

#-----------------------
# angler characteristics

ET = 100            # total fishing hours
b1 = seq(0, 1, length.out = 100)  # proportion of specialists

q1 = 0.001            # catchability of species 1
q2 = 0.05            # catchability of species 2


# all scenarios
ss = expand.grid(b1, rho)

####################
# process functions

#-----------
# ecological

# logistic growth function
Nt <- function(r, N) N + ((r * N) - (r*N^2))
Nt.c <- function(r, N, rho, N2) N + ((r * N) - (r*N^2 + rho*N*N2))

# temperature-dependent growth rate
rtemp <- function(temp, r, theta, opt_temp) r * (1 + theta * -abs(temp - opt_temp))

#rcomp - cold water r is supressed when lots of warm water
#rcomp <- function(N, r, rho) r * (1 + rho * -N)

#-------
# social

# site preference
#Nstar <- function(N1, N2, beta) N1*beta - N2*(1-beta)

# catch
catch <- function(q, E, N) q * E * N



#############
# simulations
#############

#---------------
# storage arrays
N.warm = array(NA, dim = c(lat, length(ss[,1]), years))
N.cold = array(NA, dim = c(lat, length(ss[,1]), years))

#N.star = array(NA, dim = c(lat, lon, years))

E.warm = array(NA, dim = c(lat, length(ss[,1]), years))
E.cold = array(NA, dim = c(lat, length(ss[,1]), years))

C.warm = array(NA, dim = c(lat, length(ss[,1]), years))
C.cold = array(NA, dim = c(lat, length(ss[,1]), years))


#------------------
# starting pop size
N.warm[,,1] = 0.5
N.cold[,,1] = 0.5


# simulate pop dynamics through time
for (t in 1:(years - 1)) {
  
  for(j in 1:length(ss[,1])){
    
    prop.cold = sum(N.cold[ , j, t]) / (sum(N.cold[ , j, t]) + sum(N.warm[ , j, t]))
    prop.warm = 1 - prop.cold
    
    E.cold[, j, t] <- prop.cold * (1-ss[j,1]) * ET * (N.cold[, j, t]/sum(N.cold[ , j, t])) + ET * ss[j,1] * (N.cold[, j, t]/sum(N.cold[ , j, t]))
    E.warm[, j, t] <- prop.warm * (1-ss[j,1]) * ET * (N.warm[, j, t]/sum(N.warm[ , j, t]))
    
    C.cold[, j, t] <- catch(q1, E.cold[, j, t], N.cold[, j, t])
    C.warm[, j, t] <- catch(q2, E.warm[, j, t], N.warm[, j, t])          
    
  }
  
  
  for (i in 1:lat) {
    for (j in 1:length(ss[,1])) {
      
      N.cold[i, j, t + 1] <- Nt.c(rtemp(temps[i, t], r, theta, t_opt1), N.cold[i, j, t], ss[j,2], N.warm[i, j, t]) - C.cold[i, j, t]
      N.warm[i, j, t + 1] <- Nt(rtemp(temps[i, t], r, theta, t_opt2), N.warm[i, j, t]) - C.warm[i, j, t]
      
      
    }
  }
}


# extinction years
N.colds = apply(N.cold, 3, function(x) colSums(x, na.rm = T))
ex.year = apply(N.colds, 1, function(x) which(x < 0.0001)[1])

ex.pred = cbind(ss, ex.year)
ex.pred = ex.pred %>% tidyr::spread(Var1, ex.year) %>%
  select(-Var2)

ex.pred = as.matrix(ex.pred)


png("ex_contour_CWspec.png", width = 6, height = 6, units = "in", res = 600)

filled.contour(rho, b1, ex.pred, zlim = range(100:220), color = function(n) hcl.colors(n, "YlOrRd", rev = FALSE), plot.axes = {
  axis(1)
  axis(2)
  contour(rho, b1, ex.pred, add = TRUE, lwd = 2, nlevels = 20)
},
plot.title={
  title(main="Time to Extinction",cex=5)
  title(xlab="Competition",cex.lab=1.5)
  mtext("Proportion of CW Specialists",2,cex=1.5,line=2.5,las=0)
}
)

dev.off()




#################################################################
# angler behavior in spatial fishery with two contrasting species
#################################################################

library(scales)
library(tidyr)
library(dplyr)

#
# Scenario 2: WW specialists
#

# one warm water species, one cold water species
# anglers can move anywhere - movement (and therefore fishing pressure) is dependent on preference and abundance

###################
# model parameters
###################

#---------------------
# site characteristics

years = 500        # years
lat = 10           # grid height


init.temps = seq(0, 10, length.out = lat) # temperatures across time and latitude - colder in the north
temps = array(NA, dim = c(lat, years))
temps[,1] <- init.temps

cc = years/5 # total temp increase over time in degrees

for(i in 1:lat) { temps[i,] <- temps[i,1] + seq(0,cc, length.out = years) }


#------------------------
# species characteristics

r = 0.1       # maximum pop growth rate

rho = seq(0, .5, length.out = 100)  # strength of competition

theta = 0.1   # temperature sensitivity
t_opt1 = 10   # optima for species 1
t_opt2 = 20   # optima for species 2

# differences from temperature optima
delta.t_opt1 = abs(temps - t_opt1)
delta.t_opt2 = abs(temps - t_opt2)

#-----------------------
# angler characteristics

ET = 100            # total fishing hours
b1 = seq(0, 1, length.out = 100)  # proportion of specialists

q1 = 0.001            # catchability of species 1
q2 = 0.05            # catchability of species 2


# all scenarios
ss = expand.grid(b1, rho)

####################
# process functions

#-----------
# ecological

# logistic growth function
Nt <- function(r, N) N + ((r * N) - (r*N^2))
Nt.c <- function(r, N, rho, N2) N + ((r * N) - (r*N^2 + rho*N*N2))

# temperature-dependent growth rate
rtemp <- function(temp, r, theta, opt_temp) r * (1 + theta * -abs(temp - opt_temp))

#rcomp - cold water r is supressed when lots of warm water
#rcomp <- function(N, r, rho) r * (1 + rho * -N)

#-------
# social

# site preference
#Nstar <- function(N1, N2, beta) N1*beta - N2*(1-beta)

# catch
catch <- function(q, E, N) q * E * N



#############
# simulations
#############

#---------------
# storage arrays
N.warm = array(NA, dim = c(lat, length(ss[,1]), years))
N.cold = array(NA, dim = c(lat, length(ss[,1]), years))

#N.star = array(NA, dim = c(lat, lon, years))

E.warm = array(NA, dim = c(lat, length(ss[,1]), years))
E.cold = array(NA, dim = c(lat, length(ss[,1]), years))

C.warm = array(NA, dim = c(lat, length(ss[,1]), years))
C.cold = array(NA, dim = c(lat, length(ss[,1]), years))


#------------------
# starting pop size
N.warm[,,1] = 0.5
N.cold[,,1] = 0.5


# simulate pop dynamics through time
for (t in 1:(years - 1)) {
  
  for(j in 1:length(ss[,1])){
    
    prop.cold = sum(N.cold[ , j, t]) / (sum(N.cold[ , j, t]) + sum(N.warm[ , j, t]))
    prop.warm = 1 - prop.cold
    
    E.cold[, j, t] <- prop.cold * (1-ss[j,1]) * ET * (N.cold[, j, t]/sum(N.cold[ , j, t])) 
    E.warm[, j, t] <- prop.warm * (1-ss[j,1]) * ET * (N.warm[, j, t]/sum(N.warm[ , j, t])) + ET * ss[j,1] * (N.warm[, j, t]/sum(N.warm[ , j, t]))
    
    C.cold[, j, t] <- catch(q1, E.cold[, j, t], N.cold[, j, t])
    C.warm[, j, t] <- catch(q2, E.warm[, j, t], N.warm[, j, t])          
    
  }
  
  
  for (i in 1:lat) {
    for (j in 1:length(ss[,1])) {
      
      N.cold[i, j, t + 1] <- Nt.c(rtemp(temps[i, t], r, theta, t_opt1), N.cold[i, j, t], ss[j,2], N.warm[i, j, t]) - C.cold[i, j, t]
      N.warm[i, j, t + 1] <- Nt(rtemp(temps[i, t], r, theta, t_opt2), N.warm[i, j, t]) - C.warm[i, j, t]
      
      
    }
  }
}


# extinction years
N.colds = apply(N.cold, 3, function(x) colSums(x, na.rm = T))
ex.year = apply(N.colds, 1, function(x) which(x < 0.0001)[1])

ex.pred = cbind(ss, ex.year)
ex.pred = ex.pred %>% tidyr::spread(Var1, ex.year) %>%
  select(-Var2)

ex.pred = as.matrix(ex.pred)
data_frame[,ncol(data_frame):1 ]

png("ex_contour_WWspec.png", width = 6, height = 6, units = "in", res = 600)

filled.contour(rho, b1, ex.pred[,ncol(ex.pred):1], zlim = range(100:220), color = function(n) hcl.colors(n, "YlOrRd", rev = FALSE), plot.axes = {
  axis(1)
  axis(2, at = seq(0,1,by = 0.2), labels = seq(1,0,by = -0.2))
  contour(rho, b1, ex.pred[,ncol(ex.pred):1], add = TRUE, lwd = 2, nlevels = 10)
  grid(col=rgb(1,1,1,0))
},
plot.title={
  title(main="Time to Extinction",cex=5)
  title(xlab="Competition",cex.lab=1.5)
  mtext("Proportion of WW Specialists",2,cex=1.5,line=2.5,las=0)
}
)

dev.off()





