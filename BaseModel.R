require(phaseR)
require(deSolve)
fullfitnessmodel <- function(t,y,parameters){
  I <- y[1]
  U <- y[2]
  tau <- parameters[1]
  d <- parameters[2]
  b <- parameters[3]
  D <- parameters[4]
  q <- parameters[5]
  dy <- numeric(2)
  dy[1] <- (tau*b - (d+D)*(I+U))*I
  dy[2] <- (1-tau)*b*I + (b*(1-(q*(I/(I+U))))-d*(I+U))*U
  list(dy)
}

params = c(1,1,1,0,1)
ymatrix = matrix(c(1,9,3,7,7,3,5,5,9,1,0,1), nrow=6, ncol=2, byrow=TRUE);
TT = 35
RESOLUTION = 300
WIDTH = 5
HEIGHT = 5
plot.new()
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
        xlim=c(-0.05, TT+0.5)
)

