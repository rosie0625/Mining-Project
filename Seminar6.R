library(nloptr)     # package for nonlinear optimization
dataMatrix <- as.matrix(read.table('MiningTraining2016.csv', sep = ","))
events <- dataMatrix[,1]
head(dataMatrix)
# change the time scale to minutes
events <- events/60000000
plot(events, c(1:length(events)), col = 'blue', 
     type="s",ylab="n_events",xlab="time (minutes)",
     lwd=2,xlim=c(0,max(events)),ylim=c(0,length(events)))

#estimating parameters
dt <- 0.05
w <- 10 # experiment with window width to estimate intensity
n <- 50 # experiment with number of observations of intensity in the model
t0 <- ((n+w)*dt) # earliest time when regression can be fitted

#i0 <- (findInterval(t0, events) + 1) # earliest event number 
#currTime <- events[i0]
#tGrid <- seq(currTime - t0 + dt, currTime, by=dt) # grid at t0
#head(tGrid)
#eventsGrid <- findInterval(tGrid, events) 
#eventsGrid  # eventsGrid: vector of cumulative counts of events in w+n time intervals in tGrid
#cbind(tGrid=tGrid,eventsGrid=eventsGrid,events=head(events,60))[1:16,]
#N <- length(tGrid)
#intensity <- eventsGrid[(w+1):N] - eventsGrid[1:(N-w)]
#intensity <- intensity / (dt*w) # events per minute
#timeGrid <- tGrid[(N-n+1):N] # timeGrid contains t_i, i=1,...n
# Use pmax(x, 0.1) to avoid log(0)
#logIntensity <- log(pmax(intensity, 0.1))
#exampleData<-sample(1:20,12,replace=T)
#dim(exampleData)<-c(4,3)
#exampleData
#pmax(exampleData[,1],exampleData[,2],exampleData[,3])
regression <- function(tc, logRate, regressionTimes, returnError){
  # tc - time of collapse
  # logRate - logarithm of intensity calculated at times t_i, i= 1,...,n
  # regressionTimes - sequence (tc-t_i), i= 1,...,n
  # returnError - logical flag, TRUE if function is used in a minimization procedure; then it returns mean-squared residual
  #               if FALSE return coefficient p, the slope, of the model    
  
  logregressionTimes <- log(tc-regressionTimes)
  linModel <- lm(logRate ~ logregressionTimes)
  if(returnError){ #if we use this function in a minimization procedure
    err <- sqrt(mean(linModel$residuals^2))
    res <- err    
  }
  else{ #we will need the p coefficient after minimization
    p <- linModel$coefficients[2]
    names(p) <- NULL
    res <- p
  }
  return(res)
}

#res <- nloptr(x0=currTime+1, eval_f=regression, lb=currTime+0.1, ub=currTime+10, 
 #             opts=list("algorithm"="NLOPT_LN_COBYLA", "xtol_rel" = 1e-04),
#              logRate=logIntensity, regressionTimes=timeGrid, returnError=TRUE)

#tc <- res$solution
#p <- regression(res$solution, logIntensity, timeGrid, FALSE)

pp <- numeric(104)
intensities <- numeric(104)
ct <- numeric(104)
timeToShock <- numeric(104)
i0 <- (findInterval(t0, events) + 1)
for (i in 1:104){
  currTime <- events[i0]
  tGrid <- seq(currTime - t0 + dt, currTime, by=dt)
  eventsGrid <- findInterval(tGrid, events) 
  N <- length(tGrid)
  intensity <- eventsGrid[(w+1):N] - eventsGrid[1:(N-w)]
  intensity <- intensity / (dt*w) # events per minute
  timeGrid <- tGrid[(N-n+1):N]
  logIntensity <- log(pmax(intensity, 0.1))
  res <- nloptr(x0=currTime+1, eval_f=regression, lb=currTime+0.1, ub=currTime+10, 
                opts=list("algorithm"="NLOPT_LN_COBYLA", "xtol_rel" = 1e-04),
                logRate=logIntensity, regressionTimes=timeGrid, returnError=TRUE)
  
  tc <- res$solution
  p <- regression(res$solution, logIntensity, timeGrid, FALSE)
  pp[i] <- p
  intensities[i] <- intensity[length(intensity)]
  ct[i] <- currTime
  timeToShock[i] <- tc - currTime
  i0 <- i0 + 1
}
resTable <- data.frame(ct=ct, pp=pp, intensities=intensities, timeToShock=timeToShock)
head(resTable,15)

alarm <- as.numeric(resTable$timeToShock > 1 & resTable$timeToShock < 2)
resTable <- cbind(resTable, alarm)
lr <- glm(formula = resTable$alarm ~ resTable$pp + resTable$intensities + resTable$timeToShock, family=binomial(link="logit"))
 l <- exp(predict(lr))
 g <- l/(1+l)
 plot(g)
W <- w*dt # window width in minutes
t0 <- W
i0 <- findInterval(t0, events) + 1
intensities <- c()
for(i in i0:length(events))
{
  currTime <- events[i]
  tGrid <- c(currTime - W, currTime)
  pGrid <- findInterval(tGrid, events)
  intensity <- pGrid[2] - pGrid[1]
  intensity <- intensity / (dt*w) # events per minute
  intensities <- c(intensities, intensity)
}

plot(events[i0:length(events)], intensities,
     col = 'blue', type="s",
     xlim=c(0,max(events)), ylim=c(0, max(intensities)),
     xlab='time', ylab='intenlity', lwd=2)
# Draw the dedline for the alarm signal
abline(v=events[length(events)]-1, col='red', lwd=2)
grid()

