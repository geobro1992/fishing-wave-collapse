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

rho = c(0, 0.01)  # strength of competition

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
b1 = c(0)               # proportion of coldwater specialists


q1 = 0.001            # catchability of species 1
q2 = 0.003            # catchability of species 2


# all scenarios
ss = expand.grid(b1, rho, ET)
ss[3,1] = 1.0

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

    E.cold[, j, t] <- prop.cold * (1-ss[j,1]) * ss[j,3] * (N.cold[, j, t]/sum(N.cold[ , j, t])) + ss[j,3] * ss[j,1] * (N.cold[, j, t]/sum(N.cold[ , j, t]))
    E.warm[, j, t] <- prop.warm * (1-ss[j,1]) * ss[j,3] * (N.warm[, j, t]/sum(N.warm[ , j, t]))

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




# color code by latitude
colo = colorRampPalette(c("blue", "red"))


# comparing scenarios figure

fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
  
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
  
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100
  
  x1 <- switch(pos,
               topleft     =x[1] + sw, 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}


###################3
# conceptual figures
####################

png("ecological_methods_fig.png", width = 6, height = 6, units = "in", res = 600)

par(mar = c(3,3,3,3), oma = c(3,3,3,3))

layout(matrix(c(1,1,2,3), 2, 2, byrow = T)) 

# lat gradient
colo = colorRampPalette(c("blue", "red"))

par(mar = c(2,6,2,6))

plot(1:years, temps[1,], type = "l", ylim = c(0,30), col = "blue", lwd = 2, xaxt="n", yaxt="n", ylab = "", xlab = "")

for (i in 2:lat) {
  
  lines(1:years, temps[i,], col = colo(10)[i], lwd = 2)
  
}

mtext(text = "Time",side = 1,line = 1,cex = 1.5)
mtext(text = "Temperature",side = 2,line = 1.5, cex = 1.5)
text(70, 5, "Polar", col = colo(10)[1], cex = 1.3)
text(30, 25, "Equatorial", col = colo(10)[10], cex = 1.3)

corners = par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
par(xpd = TRUE) #Draw outside plot area
text(x = corners[2]+13, y = mean(corners[3:4]), "Lattitude", cex = 1.7, srt = 270)

arrows(-7, 5, -7, 25, length = 0.15, angle = 30, code = 2, lwd = 2,
       col = par("fg"), lty = NULL, oma = T)

arrows(109, 25, 109, 5, length = 0.15, angle = 30, code = 2, lwd = 2,
       col = par("fg"), lty = NULL, oma = T)

mtext(LETTERS[1], cex=2, outer = F, side = 3, line = 1) 


# temp sensitivity
ts = seq(0, 30, length.out = 1000)

par(mar = c(2,2,2,2))

plot(ts, rtemp(ts, r, theta, t_opt1), type = "l", ylim = c(-0.15,0.15), col = "dark blue", lwd = 2, xaxt="n", yaxt="n", ylab = "", xlab = "")
lines(ts, rtemp(ts, r, theta, t_opt2), col = "dark red", lwd = 2)
lines(ts, rep(0, length(ts)), col = "dark grey", lwd = 1.5, lty = "dashed")

mtext(text = "Temperature",side = 1,line = 1,cex = 1.2)
mtext(text = "Pop. growth rate (r)",side = 2,line = 1, cex = 1.2)
text(8, 0.13, "cold-adapted\nspecies", col = "dark blue", cex = 0.8)
text(22, 0.13, "warm-adapted\nspecies", col = "dark red", cex = 0.8)

mtext(LETTERS[2], cex=2, outer = F, side = 3, line = 1) 

# competition function
N2s = seq(0, 1, length.out = 100)
rhos = c(0, 0.1, 0.2, 0.3, 0.4, 0.5)

colo = colorRampPalette(c("black", "grey"))

plot(N2s, 1- (rhos[1] * N2s), type = "l", xlim = c(0, 1.3), ylim = c(0,1), col = "black", lwd = 2, xaxt="n", yaxt="n", ylab = "", xlab = "")

for (i in 2:length(rhos)) {
  lines(N2s, 1- (rhos[i] * N2s), col = colo(6)[i], lwd = 2)
}

mtext(text = "Abundance of competitor",side = 1,line = 1,cex = 1.2)
mtext(text = "Carrying capacity (K)",side = 2,line = 1, cex = 1.2)
text(1.15, 1, expression(rho == 0.0), col = colo(6)[1], cex = 0.8)
text(1.15, 0.9, expression(rho == 0.1), col = colo(6)[2], cex = 0.8)
text(1.15, 0.8, expression(rho == 0.2), col = colo(6)[3], cex = 0.8)
text(1.15, 0.7, expression(rho == 0.3), col = colo(6)[4], cex = 0.8)
text(1.15, 0.6, expression(rho == 0.4), col = colo(6)[5], cex = 0.8)
text(1.15, 0.5, expression(rho == 0.5), col = colo(6)[6], cex = 0.8)

mtext(LETTERS[3], cex=2, outer = F, side = 3, line = 1) 

dev.off()










png("all_scenarios_detail.png", width = 11, height = 11, units = "in", res = 600)


mat <- matrix(c(1,2,3,4,
                5,6,7,8,
                9,10,11,12,
                13,14,15,16), nrow = 4, byrow = F)

op <- par(cex.main = 1.5, mar = c(5, 6, 2, 2) + 0.1, mgp = c(4.5, 1, 0), cex.lab = 1.5 , font.lab = 2, cex.axis = 1.3, bty = "n", las = 1, oma = c(1,1,5,1))
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
mtext("Time", side = 1, line = 2.5, cex = 1.2)
mtext("Cold-adapted\n Abundance", side = 2, line = 3.7, cex = 1.2)
fig_label(paste("  ", LETTERS[1]), cex=2, region="plot") 
mtext("no fishing,\n no competitor", side = 3, line = 2.5, cex = 1.5)


plot.new()
plot.window( xlim=c(-5,5), ylim=c(-5,5) )
text(0,0,"No Fishing", cex=2)
fig_label(paste("  ", LETTERS[5]), cex=2, region="plot") 

plot.new()
plot.window( xlim=c(-5,5), ylim=c(-5,5) )
text(0,0,"No Competitor", cex=2)
fig_label(paste("  ", LETTERS[9]), cex=2, region="plot") 


plot(1:lat, rev(ec), xlim = c(0,10), ylim = c(0, 40), xlab = "", ylab = "", axes = FALSE,
     pch = 21, bg = rev(colo(10)), cex = 2, col = "black")

axis(1)
axis(2)

par(las = 0)
mtext("Latitude", 1, line = 2.5, cex = 1.2)
mtext("Premature Decline", 2, line = 3, cex = 1.2, las = 0)
fig_label(paste("  ", LETTERS[13]), cex=2, region="plot") 



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
mtext("Time", side = 1, line = 2.5, cex = 1.2)
mtext("Cold-adapted\n Abundance", side = 2, line = 3.7, cex = 1.2)
fig_label(paste("  ", LETTERS[2]), cex=2, region="plot") 
mtext("fishing,\n no competitor", side = 3, line = 2.5, cex = 1.5)




plot(1:years, E.cold[1,3,], type = "l", xlim = c(0,100), ylim = c(0,20), col = "blue", lwd = 2, ylab = "", xlab = "", axes = FALSE)

for (i in 2:lat) {
  
  lines(1:years, E.cold[i,3,], col = colo(10)[i], lwd = 2)
  
}


axis(1)
axis(2) 

par(las = 0)
mtext("Time", side = 1, line = 2.5, cex = 1.2)
mtext("Fishing Effort", side = 2, line = 3.7, cex = 1.2)
fig_label(paste("  ", LETTERS[6]), cex=2, region="plot") 


plot.new()
plot.window( xlim=c(-5,5), ylim=c(-5,5) )
text(0,0,"No Competitor", cex=2)
fig_label(paste("  ", LETTERS[10]), cex=2, region="plot") 


plot(1:lat, rev(ec), xlim = c(0,10), ylim = c(0, 40), xlab = "", ylab = "", axes = FALSE,
     pch = 21, bg = rev(colo(10)), cex = 2, col = "black")

axis(1)
axis(2)

par(las = 0)
mtext("Latitude", 1, line = 2.5, cex = 1.2)
mtext("Premature Decline", 2, line = 3, cex = 1.2, las = 0)

abline(h = mean(ec), lty = 2, lwd = 2)
fig_label(paste("  ", LETTERS[14]), cex=2, region="plot") 



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
mtext("Time", side = 1, line = 3, cex = 1.2)
mtext("Cold-adapted\n Abundance", side = 2, line = 3.7, cex = 1.2)
fig_label(paste("  ", LETTERS[3]), cex=2, region="plot") 
mtext("competitor,\n no fishing", side = 3, line = 2.5, cex = 1.5)



plot.new()
plot.window( xlim=c(-5,5), ylim=c(-5,5) )
text(0,0,"No Fishing", cex=2)
fig_label(paste("  ", LETTERS[7]), cex=2, region="plot") 


plot(1:years, N.warm[1,2,], type = "l", xlim = c(0,100), ylim = c(0,1), col = "blue", lwd = 2, ylab = "", xlab = "", axes = FALSE)

for (i in 2:lat) {
  
  lines(1:years, N.warm[i,2,], col = colo(10)[i], lwd = 2)
  
}


axis(1)
axis(2) 

par(las = 0)
mtext("Time", side = 1, line = 3, cex = 1.2)
mtext("Competitor\n Abundance", side = 2, line = 3.7, cex = 1.2)
fig_label(paste("  ", LETTERS[11]), cex=2, region="plot") 



plot(1:lat, rev(ec), xlim = c(0,10), ylim = c(0, 40), xlab = "", ylab = "", axes = FALSE,
     pch = 21, bg = rev(colo(10)), cex = 2, col = "black")

axis(1)
axis(2)

par(las = 0)
mtext("Latitude", 1, line = 3, cex = 1.2)
mtext("Premature Decline", 2, line = 3, cex = 1.2, las = 0)

abline(h = mean(ec), lty = 2, lwd = 2)
fig_label(paste("  ", LETTERS[15]), cex=2, region="plot") 


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
mtext("Time", side = 1, line = 3, cex = 1.2)
mtext("Cold-adapted\n Abundance", side = 2, line = 3.7, cex = 1.2)
fig_label(paste("  ", LETTERS[4]), cex=2, region="plot") 
mtext("both fishing\n and competitor", side = 3, line = 2.5, cex = 1.5)





plot(1:years, E.cold[1,4,], type = "l", xlim = c(0,100), ylim = c(0,20), col = "blue", lwd = 2, ylab = "", xlab = "", axes = FALSE)

for (i in 2:lat) {
  
  lines(1:years, E.cold[i,4,], col = colo(10)[i], lwd = 2)
  
}


axis(1)
axis(2) 

par(las = 0)
mtext("Time", side = 1, line = 3, cex = 1.2)
mtext("Fishing Effort", side = 2, line = 3.7, cex = 1.2)
fig_label(paste("  ", LETTERS[8]), cex=2, region="plot") 


plot(1:years, N.warm[1,4,], type = "l", xlim = c(0,100), ylim = c(0,1), col = "blue", lwd = 2, ylab = "", xlab = "", axes = FALSE)

for (i in 2:lat) {
  
  lines(1:years, N.warm[i,4,], col = colo(10)[i], lwd = 2)
  
}


axis(1)
axis(2) 

par(las = 0)
mtext("Time", side = 1, line = 3, cex = 1.2)
mtext("Competitor\n Abundance", side = 2, line = 3.7, cex = 1.2)
fig_label(paste("  ", LETTERS[12]), cex=2, region="plot") 



plot(1:lat, rev(ec), xlim = c(0,10), ylim = c(0, 40), xlab = "", ylab = "", axes = FALSE,
     pch = 21, bg = rev(colo(10)), cex = 2, col = "black")

axis(1)
axis(2)

par(las = 0)
mtext("Latitude", 1, line = 3, cex = 1.2)
mtext("Premature Decline", 2, line = 3, cex = 1.2, las = 0)

abline(h = mean(ec), lty = 2, lwd = 2)
fig_label(paste("  ", LETTERS[16]), cex=2, region="plot") 


dev.off()




####################################
# cummulative extinctions over time

c.ex = data.frame(n = rep(10:1, 2), scenario = rev(sort(rep(c("fishing", "no fishing"), 10))), y = NA)

# no fishing or competition

for (i in 1:lat) {
  
  c.ex[i, "y"] = t_crit1.years[i] 
  c.ex[i+10, "y"] = which.max(N.cold[i,3,])

}

c.ex$scenario = as.factor(c.ex$scenario)

library(ggplot2)
ggplot(c.ex, aes(y, n, color = scenario))+
  geom_line() +
  xlab("Year") +
  ylab("Cumulative Extinction") +
  theme_classic()


ggplot(c.ex, aes(x = y, y = n, group = scenario, 
                     color = scenario)) +
  geom_step(linewidth = 1.5, direction = "hv") +
  geom_point(size = 3) + 
  xlim(0, 100) +
  theme_classic()
