require(phaseR)
require(deSolve)
fullfitnessmodel <- function(t,y,parameters){
  I <- y[1]
  U <- y[2]
  tau <- parameters[1]
  d <- parameters[2]
  a <- parameters[3]
  B <- parameters[4]
  D <- parameters[5]
  q <- parameters[6]
  dy <- numeric(2)
  dy[1] <- (tau*(1/(exp(-(I+U)*B)+a)) - (d+D)*(I+U))*I
  dy[2] <- (1-tau)*(1/(exp(-(I+U)*B)+a))*I + ((1/(exp(-(I+U)*B)+a))*(1-(q*(I/(I+U))))-d*(I+U))*U
  list(dy)
}

params = c(1,1,.1,3,0,1)
ymatrix = matrix(c(5,5,3,7,9,1,7,3,1,9), nrow=5, ncol=2, byrow=TRUE);
TT = 25
RESOLUTION = 300
WIDTH = 5
HEIGHT = 5
#plot.new()
fullfitnessmodel.flowField <-
  flowField(fullfitnessmodel, parameters=params,x.lim=c(-0.05,11),y.lim=c(-0.05,11),
            points=15, add=FALSE,col='chartreuse4',xlab='Wb-Infected Mosquitoes',ylab='Wb-Uninfected Mosquitoes',main=NULL)
grid()
fullfitnessmodel.nullclines <-
  nullclines(fullfitnessmodel,x.lim=c(0,11),y.lim=c(0,11), lwd=2,lty=1, 
             colour=c("steelblue1", "springgreen4"), parameters = params, points=500)
y0 <- ymatrix
fullfitnessmodel.trajectory <-
  trajectory(fullfitnessmodel, y0=y0, t.end=TT,parameters=params,colour=c('turquoise','seagreen3','palegreen3','lightseagreen','deepskyblue','green'),lwd=3,pch=20,cex=2)
dev.off
time = seq(0,TT, by=0.1)
TMAT= matrix(0,nrow=length(time),ncol=dim(ymatrix)[1])
IMAT=TMAT
UMAT=IMAT
for (i in 1:dim(ymatrix)[1]) {
  states = ymatrix[i,]
  solutions = ode(states, time, fullfitnessmodel, params)
  TMAT[,i]=solutions[,1]
  IMAT[,i]=solutions[,2]
  UMAT[,i]=solutions[,3]
}
#plot.new()
#par(mar=c(4, 4, 2, 2))
par(mar=c(4, 4, 2, 2))
matplot(TMAT[,1],IMAT,type='l',lwd=2, lty=1.1,
        xlab='time',
        ylab='Wb-Infected Mosquitos',
        cex.lab=1.25,
        cex.axis=1.1,
        ylim=c(-0.05,11),
        xlim=c(-0.05, 10)
)

