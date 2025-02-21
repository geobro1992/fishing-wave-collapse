#################################################################
# angler behavior in spatial fishery with two contrasting species
#################################################################

library(scales)
library(tidyr)

#
# Scenario 1
#

# one warm water species, one cold water species
# anglers can move anywhere - movement (and therefore fishing pressure) is dependent on preference and abundance

###################
# model parameters
###################

#---------------------
# site characteristics

years = 100        # years
lat = 10           # grid height


init.temps = seq(0, 10, length.out = lat) # temperatures across time and latitude - colder in the north
temps = array(NA, dim = c(lat, years))
temps[,1] <- init.temps

cc = 20 # total temp increase over time in degrees

for(i in 1:lat) { temps[i,] <- temps[i,1] + seq(0,cc, length.out = years) }


#------------------------
# species characteristics

r = 0.1       # maximum pop growth rate

rho = c(0, 0.05)  # strength of competition

theta = 0.1   # temperature sensitivity
t_opt1 = 10   # optima for species 1
t_opt2 = 20   # optima for species 2

# critical temperature (where r < 0)
t_crit1 = t_opt1 + (1/theta)
t_crit2 = t_opt2 + (1/theta)

# differences from temperature optima
delta.t_opt1 = abs(temps - t_opt1)
delta.t_opt2 = abs(temps - t_opt2)
t_opt1.years = apply(delta.t_opt1, 1, which.min)
t_opt2.years = apply(delta.t_opt2, 1, which.min)

# differences from critical temperature 
delta.t_crit1 = abs(temps - t_crit1)
delta.t_crit2 = abs(temps - t_crit2)
t_crit1.years = apply(delta.t_crit1, 1, which.min)
t_crit2.years = apply(delta.t_crit2, 1, which.min)

#-----------------------
# angler characteristics

ET = c(0, 100)            # total fishing hours
b1 = c(0)               # proportion of warmwater specialists


q1 = 0.001            # catchability of species 1
q2 = 0.003            # catchability of species 2


# all scenarios
ss = expand.grid(b1, rho, ET)
ss[4,1] = 1.0

####################
# process functions

#-----------
# ecological

# logistic growth function
Nt <- function(r, N) N + ((r * N) - (r*N^2))
Nt.c <- function(r, N, rho, N2) N + ((r * N) - (r*N^2 + rho*N*N2))

# temperature-dependent growth rate
rtemp <- function(temp, r, theta, opt_temp) r * (1 + theta * -abs(temp - opt_temp))


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
    
    E.cold[, j, t] <- prop.cold * (1-ss[j,1]) * ss[j,3] * (N.cold[, j, t]/sum(N.cold[ , j, t])) 
    E.warm[, j, t] <- prop.warm * (1-ss[j,1]) * ss[j,3] * (N.warm[, j, t]/sum(N.warm[ , j, t])) + ss[j,3] * ss[j,1] * (N.warm[, j, t]/sum(N.warm[ , j, t]))
    
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



# plot 2X2 sims to show the effect of fishing and competition

# color code by latitude
colo = colorRampPalette(c("blue", "red"))


png("2x2_ww_specialists.png", width = 6, height = 6, units = "in", res = 600)

par(mfrow = c(2, 2), mar = c(2,2,2,2), oma = c(3, 3, 3, 3))

for (j in 1:length(ss[,1])) {
  
  plot(1:years, N.cold[1,j,], type = "l", ylim = c(0,1), col = "blue", lwd = 2, xaxt="n", yaxt="n", ylab = "", xlab = "")
  points(t_crit1.years[1], N.cold[1,j,t_crit1.years[1]], col = colo(10)[1], pch = 19, cex = 1.5)
  
  for (i in 2:lat) {
    
    lines(1:years, N.cold[i,j,], col = colo(10)[i], lwd = 2)
    points(t_crit1.years[i], N.cold[i,j,t_crit1.years[i]], col = colo(10)[i], pch = 19, cex = 1.5)
  }
}

mtext(text = "No Competition",side = 3,line = 0,outer = TRUE, cex = 1.5, adj = 0.15)
mtext(text = "Competition",side = 3,line = 0,outer = TRUE, cex = 1.5, adj = 0.8)
mtext(text = "Fishing",side = 2,line = -0.5,outer = TRUE, cex = 1.5, adj = 0.2)
mtext(text = "No fishing",side = 2,line = 1,outer = TRUE, cex = 1.5, adj = 0.85)

dev.off()



###########################
# both and early collapse


ec = vector()

for (i in 1:lat) {
  
  ec[i] = t_crit1.years[i] - which.max(N.cold[i,4,])
  
}




png("chasing_both_ww_specialists.png", width = 7, height = 6, units = "in", res = 600)


mat <- matrix(c(0,0,1,1,0,0,
                2,2,3,3,4,4), nrow = 2, byrow = TRUE)

op <- par(cex.main = 1.5, mar = c(5, 6, 2, 2) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5 , font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
layout(mat)

plot(1:years, N.cold[1,4,], type = "l", xlim = c(0,100), ylim = c(0,1), col = "blue", lwd = 2, ylab = "", xlab = "", axes = FALSE)
points(t_crit1.years[1]-1, N.cold[1,4,99], col = colo(10)[1], pch = 19, cex = 1.5)

for (i in 2:lat) {
  
  lines(1:years, N.cold[i,4,], col = colo(10)[i], lwd = 2)
  points(t_crit1.years[i], N.cold[i,4,t_crit1.years[i]], col = colo(10)[i], pch = 19, cex = 1.5)
}


axis(1)
axis(2) 

par(las = 0)
mtext("Time", side = 1, line = 3, cex = 1.5)
mtext("Abundance", side = 2, line = 3.7, cex = 1.5)




plot(1:years, N.warm[1,4,], type = "l", xlim = c(0,100), ylim = c(0,1), col = "blue", lwd = 2, ylab = "", xlab = "", axes = FALSE)

for (i in 2:lat) {
  
  lines(1:years, N.warm[i,4,], col = colo(10)[i], lwd = 2)
  
}


axis(1)
axis(2) 

par(las = 0)
mtext("Time", side = 1, line = 3, cex = 1.5)
mtext("Competitor", side = 2, line = 3.7, cex = 1.5)



plot(1:years, E.warm[1,4,], type = "l", xlim = c(0,100), ylim = c(0,20), col = "blue", lwd = 2, ylab = "", xlab = "", axes = FALSE)

for (i in 2:lat) {
  
  lines(1:years, E.warm[i,4,], col = colo(10)[i], lwd = 2)
  
}


axis(1)
axis(2) 

par(las = 0)
mtext("Time", side = 1, line = 3, cex = 1.5)
mtext("Fishing Effort", side = 2, line = 3.7, cex = 1.5)



plot(1:lat, rev(ec), xlim = c(0,10), ylim = c(0, 40), xlab = "", ylab = "", axes = FALSE,
     pch = 21, bg = rev(colo(10)), cex = 2, col = "black")

axis(1)
axis(2)

par(las = 0)
mtext("Latitude", 1, line = 3, cex = 1.5)
mtext(expression("Early Decline "*(T[crit] - T[max])), 2, line = 3, cex = 1.5, las = 0)

abline(h = mean(ec), lty = 2, lwd = 2)


dev.off()






png("all_scenarios_detail.png", width = 11, height = 11, units = "in", res = 600)


mat <- matrix(c(1,2,3,4,
                5,6,7,8,
                9,10,11,12,
                13,14,15,16), nrow = 4, byrow = F)

op <- par(cex.main = 1.5, mar = c(5, 6, 2, 2) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5 , font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
layout(mat)


###########################
# no fishing or competition

ec = vector()

for (i in 1:lat) {
  
  ec[i] = 0
  
}



plot(1:years, N.cold[1,1,], type = "l", xlim = c(0,100), ylim = c(0,1), col = "blue", lwd = 2, ylab = "", xlab = "", axes = FALSE)
points(t_crit1.years[1]-1, N.cold[1,1,99], col = colo(10)[1], pch = 19, cex = 1.5)

for (i in 2:lat) {
  
  lines(1:years, N.cold[i,1,], col = colo(10)[i], lwd = 2)
  points(t_crit1.years[i], N.cold[i,1,t_crit1.years[i]], col = colo(10)[i], pch = 19, cex = 1.5)
}


axis(1)
axis(2) 

par(las = 0)
mtext("Time", side = 1, line = 2.5, cex = 1.5)
mtext("Abundance", side = 2, line = 3.7, cex = 1.5)



plot.new()
plot.window( xlim=c(-5,5), ylim=c(-5,5) )
text(0,0,"No Fishing", cex=2)

plot.new()
plot.window( xlim=c(-5,5), ylim=c(-5,5) )
text(0,0,"No Competitor", cex=2)


plot(1:lat, rev(ec), xlim = c(0,10), ylim = c(0, 40), xlab = "", ylab = "", axes = FALSE,
     pch = 21, bg = rev(colo(10)), cex = 2, col = "black")

axis(1)
axis(2)

par(las = 0)
mtext("Latitude", 1, line = 2.5, cex = 1.5)
mtext("Premature Decline", 2, line = 3, cex = 1.5, las = 0)



###########################
# effort and early collapse

ec = vector()

for (i in 1:lat) {
  
  ec[i] = t_crit1.years[i] - which.max(N.cold[i,3,])
  
}



plot(1:years, N.cold[1,3,], type = "l", xlim = c(0,100), ylim = c(0,1), col = "blue", lwd = 2, ylab = "", xlab = "", axes = FALSE)
points(t_crit1.years[1]-1, N.cold[1,3,99], col = colo(10)[1], pch = 19, cex = 1.5)

for (i in 2:lat) {
  
  lines(1:years, N.cold[i,3,], col = colo(10)[i], lwd = 2)
  points(t_crit1.years[i], N.cold[i,3,t_crit1.years[i]], col = colo(10)[i], pch = 19, cex = 1.5)
}


axis(1)
axis(2) 

par(las = 0)
mtext("Time", side = 1, line = 2.5, cex = 1.5)
mtext("Abundance", side = 2, line = 3.7, cex = 1.5)




plot(1:years, E.cold[1,3,], type = "l", xlim = c(0,100), ylim = c(0,20), col = "blue", lwd = 2, ylab = "", xlab = "", axes = FALSE)

for (i in 2:lat) {
  
  lines(1:years, E.cold[i,3,], col = colo(10)[i], lwd = 2)
  
}


axis(1)
axis(2) 

par(las = 0)
mtext("Time", side = 1, line = 2.5, cex = 1.5)
mtext("Fishing Effort", side = 2, line = 3.7, cex = 1.5)


plot.new()
plot.window( xlim=c(-5,5), ylim=c(-5,5) )
text(0,0,"No Competitor", cex=2)


plot(1:lat, rev(ec), xlim = c(0,10), ylim = c(0, 40), xlab = "", ylab = "", axes = FALSE,
     pch = 21, bg = rev(colo(10)), cex = 2, col = "black")

axis(1)
axis(2)

par(las = 0)
mtext("Latitude", 1, line = 2.5, cex = 1.5)
mtext("Premature Decline", 2, line = 3, cex = 1.5, las = 0)

abline(h = mean(ec), lty = 2, lwd = 2)



###########################
# competititon and early collapse


ec = vector()

for (i in 1:lat) {
  
  ec[i] = t_crit1.years[i] - which.max(N.cold[i,2,])
  
}


plot(1:years, N.cold[1,2,], type = "l", xlim = c(0,100), ylim = c(0,1), col = "blue", lwd = 2, ylab = "", xlab = "", axes = FALSE)
points(t_crit1.years[1]-1, N.cold[1,2,99], col = colo(10)[1], pch = 19, cex = 1.5)

for (i in 2:lat) {
  
  lines(1:years, N.cold[i,2,], col = colo(10)[i], lwd = 2)
  points(t_crit1.years[i], N.cold[i,2,t_crit1.years[i]], col = colo(10)[i], pch = 19, cex = 1.5)
}


axis(1)
axis(2) 

par(las = 0)
mtext("Time", side = 1, line = 3, cex = 1.5)
mtext("Abundance", side = 2, line = 3.7, cex = 1.5)



plot.new()
plot.window( xlim=c(-5,5), ylim=c(-5,5) )
text(0,0,"No Fishing", cex=2)


plot(1:years, N.warm[1,2,], type = "l", xlim = c(0,100), ylim = c(0,1), col = "blue", lwd = 2, ylab = "", xlab = "", axes = FALSE)

for (i in 2:lat) {
  
  lines(1:years, N.warm[i,2,], col = colo(10)[i], lwd = 2)
  
}


axis(1)
axis(2) 

par(las = 0)
mtext("Time", side = 1, line = 3, cex = 1.5)
mtext("Competitor", side = 2, line = 3.7, cex = 1.5)



plot(1:lat, rev(ec), xlim = c(0,10), ylim = c(0, 40), xlab = "", ylab = "", axes = FALSE,
     pch = 21, bg = rev(colo(10)), cex = 2, col = "black")

axis(1)
axis(2)

par(las = 0)
mtext("Latitude", 1, line = 3, cex = 1.5)
mtext("Premature Decline", 2, line = 3, cex = 1.5, las = 0)

abline(h = mean(ec), lty = 2, lwd = 2)


###########################
# both and early collapse


ec = vector()

for (i in 1:lat) {
  
  ec[i] = t_crit1.years[i] - which.max(N.cold[i,4,])
  
}



plot(1:years, N.cold[1,4,], type = "l", xlim = c(0,100), ylim = c(0,1), col = "blue", lwd = 2, ylab = "", xlab = "", axes = FALSE)
points(t_crit1.years[1]-1, N.cold[1,4,99], col = colo(10)[1], pch = 19, cex = 1.5)

for (i in 2:lat) {
  
  lines(1:years, N.cold[i,4,], col = colo(10)[i], lwd = 2)
  points(t_crit1.years[i], N.cold[i,4,t_crit1.years[i]], col = colo(10)[i], pch = 19, cex = 1.5)
}


axis(1)
axis(2) 

par(las = 0)
mtext("Time", side = 1, line = 3, cex = 1.5)
mtext("Abundance", side = 2, line = 3.7, cex = 1.5)





plot(1:years, E.cold[1,4,], type = "l", xlim = c(0,100), ylim = c(0,20), col = "blue", lwd = 2, ylab = "", xlab = "", axes = FALSE)

for (i in 2:lat) {
  
  lines(1:years, E.cold[i,4,], col = colo(10)[i], lwd = 2)
  
}


axis(1)
axis(2) 

par(las = 0)
mtext("Time", side = 1, line = 3, cex = 1.5)
mtext("Fishing Effort", side = 2, line = 3.7, cex = 1.5)


plot(1:years, N.warm[1,4,], type = "l", xlim = c(0,100), ylim = c(0,1), col = "blue", lwd = 2, ylab = "", xlab = "", axes = FALSE)

for (i in 2:lat) {
  
  lines(1:years, N.warm[i,4,], col = colo(10)[i], lwd = 2)
  
}


axis(1)
axis(2) 

par(las = 0)
mtext("Time", side = 1, line = 3, cex = 1.5)
mtext("Competitor", side = 2, line = 3.7, cex = 1.5)



plot(1:lat, rev(ec), xlim = c(0,10), ylim = c(0, 40), xlab = "", ylab = "", axes = FALSE,
     pch = 21, bg = rev(colo(10)), cex = 2, col = "black")

axis(1)
axis(2)

par(las = 0)
mtext("Latitude", 1, line = 3, cex = 1.5)
mtext("Premature Decline", 2, line = 3, cex = 1.5, las = 0)

abline(h = mean(ec), lty = 2, lwd = 2)


dev.off()


