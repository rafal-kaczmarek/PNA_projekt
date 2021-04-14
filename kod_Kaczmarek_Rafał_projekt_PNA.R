setwd("D:\\studia\\IV rok\\Programowanie narzêdzi analitycznych\\projekt")
library("resampledata")
library("maxLik")
library("rgl")

#Weibull
dane = read.csv(file="hard_drive.csv", header=T, sep=";")
View(dane)
dane1=dane[which(dane$lifetime>0),]
View(dane1)

x=dane1$lifetime
summary(x)
N=length(x)

hist(x, main='¯ywotnoœæ dysków twardych (w dniach)',prob=TRUE)


## Rozklad wykladniczy

N=length(x)

lnL=function(lambda){
  l=N*log(lambda)-lambda*sum(x)
  return(l)
}
gradient = function(lambda){
  g= N/lambda-sum(x)
  return(g)
}

hesjan = function(lambda){
  h=-N/lambda^2
  return(h)
}

wynik=maxNR(fn=lnL,grad=gradient,hess=hesjan,start=10)
summary(wynik)

curve(lnL(x), from=0, to=2300)
abline(v=wynik$estimate,col="green")

## Testowanie hipotezy zakadajacej rok bezawaryjnej pracy.

lambda= 1/365
(std.err.p=sqrt(-solve(wynik$hessian)))
z.test=(wynik$estimate-lambda)/std.err.p
print(z.test)
p.value=2*(1-pnorm(abs(z.test),mean=0,sd=1))
print(p.value)
#p-value bliskie 0, co oznacza, ¿e odrzucam hipoteze zerowa

1/wynik$estimate
#Z badania wynika, ¿e lambda=0.0032159 co oznacza, ¿e przecietny czas dzialania maszyny to 310,95 dni


## Rozklad Weibulla
lnL=function(parametry){
  k=parametry[1]
  lambda=parametry[2]
  ll=N*log(k)-k*N*log(lambda)+(k-1)*sum(log(x))-sum((x/lambda)^k)
  return(ll)
}

gradient = function(parametry){
  k=parametry[1]
  lambda=parametry[2]
  gr=rep(0,times=2)
  gr[1]=N/k-N*log(lambda)+sum(log(x))-sum(((x/lambda)^k)*log(x/lambda))
  gr[2]=-N*k/lambda+(k/lambda)*sum((x/lambda)^k)
  return(gr)
}

hesjan = function(parametry){
  k=parametry[1]
  lambda=parametry[2]
  h=matrix(0,2,2)
  h[1,1]=-N/(k^2)-sum(((x/lambda)^k)*(log(x/lambda))^2)
  h[2,2]=k*N/(lambda^2)-(k+1)*k/(lambda^(k+2))*sum((x)^k)
  h[1,2]=-N/lambda+1/lambda*sum((x/lambda)^k)+k/lambda*sum(((x/lambda)^k)*log(x/lambda))
  h[2,1]=h[1,2]
  return(h)
}
wynik=maxNR(fn=lnL,grad=gradient, hess=hesjan, start=c(10,10))
summary(wynik)

wynik=maxNR(fn=lnL, start=c(10,10))
summary(wynik)

#H0 k=1 czy jest wykladniczy
std.err.k=sqrt(-solve(wynik$hessian))[1,1]
z.test=(wynik$estimate[1]-1)/std.err.k
p.value=2*(1-pnorm(abs(z.test),mean=0,sd=1))
print(p.value)

#H0 k=2 - sprawdzam (rozk³ad Rayleigha) prawdopodobieñstwo roœnie liniowo z czasem,
std.err.k=sqrt(-solve(wynik$hessian))[1,1]
z.test=(wynik$estimate[1]-2)/std.err.k
p.value=2*(1-pnorm(abs(z.test),mean=0,sd=1))
print(p.value)

#sprawdzam czy rozkld z tego zjawiska mozna porownywac z rozkladem z zjawiska dotyczacego obrabiarek
#H0: k=34.05 and lambda=945.25
#test LR
lnL_U=wynik$maximum
lnL_R=lnL(c(34.05,945.25))
LRtest=2*(lnL_U-lnL_R)
g=2
p.value=1-pchisq(q=LRtest, df=g)
print(p.value)
#p-value = 0 < alpha=5% wiec odrzucam H0


#Dla sprawdzenia drugi sposob

install.packages("numDeriv")
library(numDeriv)

hess <- hessian(func=lnL, x=x)
hessc <- hessian(func=lnL, x=x, "complex")
all.equal(hess, hessc, tolerance = .Machine$double.eps)


vec <- x
weibull_loglik <- function(parm){
  n <- length(vec)
  gamma <- parm[1]
  lambda <- parm[2]
  loglik <- sum(dweibull(vec, shape=gamma, scale=lambda, log=TRUE))
  return(-loglik)
}
weibull <- nlm(weibull_loglik, p = c(1,1), hessian=TRUE)
summary(weibull)
weibull$estimate
weibull$hessian
weibull$gradient



