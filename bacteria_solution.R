
rm(list=ls())	#clean memory
graphics.off() 	#close plot windows



##########################
### Creating the model ###
##########################

### time ###
times1 = seq(0,50,by=0.5)

### initial values ###
state1 = c(Bact = 0.1, Sub = 100) 

### diferential equations ###
equations = function(t, state, pars) 			# returns rate of change 
		{ 
		with(as.list(c(state, pars)),{ 
				V1 = gmax*Sub/(Sub+ks)*Bact
				V2 = Bac_death_rate*Bact

				dBact = eff * V1 - V2
				dSub = - V1 + canibal_rate * V2 
				return(list(c(dBact, dSub), TOC = Bact + Sub)) 
				}) 
		} 

### parameters ###
pars1 = list(gmax = 0.5, eff = 0.5, ks = 0.5, Bac_death_rate = 0.02, canibal_rate=0.5)

### Solve the equations ###
library(deSolve)
out = ode(y=state1, times=times1, func = equations, parms=pars1)

### Plot the time course ###
plot(out)
matplot(out[,1], out[,-1], type = "l", lty = 1:3, lwd = c(2, 2, 1), 
 	col = "black", xlab = "time, hour", ylab = "mol C/m3") 

legend("topright", c("Bacteria", "Glucose", "TOC"), lty = 1:3, lwd = c(2, 2, 1))
par(mfrow=c(1,1))





##############
### Events ###
##############

### Creating the events ###
event_data = data.frame(
                    var = c("Sub"),
                    time = c(40),
                    value = c(20),
                    method = c("add")
                    )
### Running the simulation with events ###
out <- ode(y=state1, times=times1, func = equations, parms=pars1,
          events = list(data = event_data))
plot(out)
par(mfrow=c(1,1))





#########################
### Local sensitivity ###
#########################
library(FME)

### time ###
times1 = seq(0,50,by=0.5)

### Creating a function that accepts pars and the output is variable values ###
solveBact = function(pars) 
					{
					return(ode(y = state1, times = times1, func = equations, parms = pars)) 	
					}	 

### storing the local sensitivity information for Bact variable ###
SnsBact<- sensFun(func = solveBact, parms = pars1, sensvar = "Bact", varscale = 1) 
head(SnsBact)
### summary for Bact variable local sensitivity ###
summary(SnsBact)
plot(SnsBact)

### summary for all variables local sensitivity ###
summary(sensFun(solveBact, pars1, varscale = 1), var = TRUE)





##########################
### global sensitivity ###
##########################
library(FME)

### Creating a function that accepts pars and the output is variable values ###
solveBact = function(pars) 
					{
					return(ode(y = state1, times = times1, func = equations, parms = pars)) 	
					}	 

### Set the variability intervals for three parameters "gmax", "eff" and "Bac_death_rate" 
parRanges <- data.frame(min = c(0.4, 0.4, 0.0), max = c(0.6, 0.6, 0.02))
rownames(parRanges) <- c("gmax", "eff", "Bac_death_rate")
parRanges

### storing the global sensitivity for "Bac_death_rate" parameter
print(system.time( 
			sR <- sensRange(func = solveBact, parms = pars1, dist = "grid", 
			sensvar = c("Bact", "Sub"), parRange = parRanges[3,], num = 50) 
			))
head(summary(sR))

### print the results ###
par(mfrow=c(2, 2))
plot(summary(sR), xlab = "time, hour", ylab = "molC/m3", 
	legpos = "topright", mfrow = NULL) 
plot(summary(sR), xlab = "time, hour", ylab = "molC/m3", mfrow = NULL, 
	quant = TRUE, col = c("lightblue", "darkblue"), legpos = "topright")
mtext(outer = TRUE, line = -1.5, side = 3, "Sensitivity to Bac_death_rate", cex = 1.25)
par(mfrow = c(1, 1))





### global sensitivity for several parameters ###
Sens2 <- summary(sensRange(func = solveBact, parms = pars1, 
		dist = "latin", sensvar = "Bact", parRange = parRanges, num = 100))
plot(Sens2, main = "Sensitivity gmax,eff,Bac_death_rate", xlab = "time, hour", ylab = "molC/m3")





#########################
### fit model to data ###
#########################

### create fake data ###
pars_fake = pars1
pars_fake$gmax=.40
pars_fake$ks=.75
out_fake = ode(y=state1, times=times1, func = equations, parms=pars_fake)
plot(out_fake)
fake_data = cbind(c(10,20,30,40,50),out_fake[out_fake[,1]==10|out_fake[,1]==20|out_fake[,1]==30|out_fake[,1]==40|out_fake[,1]==50,2])
fake_data
### visualising fake data ###
out = ode(y=state1, times=times1, func = equations, parms=pars1)
plot(out[,1],out[,2])
points(fake_data,col='blue',pch=2)


### error function ###
Error_Function = function (pars_error)
	{
	out_error = ode(y=state1, times=times1, func = equations, parms=pars_error)
	p1 = abs( fake_data[1,2] - out_error[out_error[,1]==10,2])
	p2 = abs( fake_data[2,2] - out_error[out_error[,1]==20,2])
	p3 = abs( fake_data[3,2] - out_error[out_error[,1]==30,2])
	p4 = abs( fake_data[4,2] - out_error[out_error[,1]==40,2])
	p5 = abs( fake_data[5,2] - out_error[out_error[,1]==50,2])
	sum(p1,p2,p3,p4,p5)
	}
# Error_Function(pars1)



op = optim(
	par=pars1,	 
	method = "Nelder-Mead",
	fn = Error_Function,
	control = list(
			trace=1
			)
	)

# "Nelder-Mead" "BFGS"

op
cbind(pars1, pars_fake, pars_op=op$par)

pars_op = op$par
out_op = ode(y=state1, times=times1, func = equations, parms=pars_op)
par(mfrow=c(1,1))
plot(out_op[,1],out_op[,2])
points(fake_data,col='blue',pch=2)
Error_Function(pars_op)

