########################################################################
#### Packages and working directory                                 ####
########################################################################
library(fda)
library(fds)
library(fda.usc)
library(janitor)
#setwd("D:/documents/2020-2021 m.m/Pavasaris 2021/FDA/Pratybos/Week 9")

########################################################################
#### One sample pointwise-test                                      ####
########################################################################

# Test statistics, confidence interval and plots
source("Ztest.R")

########################################################################
# Age-specific mortality rates for Australia and Australian states.
# Example: compare for male data in Australia with female mean
# H0: mu(male) = mu(female)
# H1: mu(male) != mu(female)
opar <- par(mfrow=c(2,1))
plot(ausmale)
plot(ausfemale)
par(opar)

class(ausmale)

dta.A <- fdata2fd(fdata(ausmale))
t.sq <- seq(0,100, length=501)
mean.female <- mean.fd(fdata2fd(fdata(ausfemale)))
stat <- Z(x=dta.A, t.seq = t.sq, mu=mean.female)

# Tasmania data
# Example: compare for Tasmanian male data in Australia with Australian total male mean
# H0: mu(Tasmania) = mu(Australia)
# H1: mu(Tasmania) != mu(Australia)
plot(tasmale)

# NaN values in coefficient matrix
dta.T <- fdata2fd(fdata(tasmale))
dta.T

dta.T.coef <- remove_empty(dta.T$coefs, which = "cols", quiet = FALSE)
dta.T$fdnames$reps <- colnames(dta.T.coef)
dta.T.fd <- fd(dta.T.coef, dta.T$basis, dta.T$fdnames)
plot(dta.T.fd)

t.sq <- seq(0,100, length=501)
mean.male <- mean.fd(fdata2fd(fdata(ausmale)))
stat <- Z(x=dta.T.fd, t.seq = t.sq, mu=mean.male)

# Victoria data
# Example: compare for Victoria female data in Australia with Australian total female mean
# H0: mu(Victoria) = mu(Australia)
# H1: mu(Victoria) != mu(Australia)
plot(vicfemale)

# NaN values in coefficient matrix
dta <- fdata2fd(fdata(vicfemale))
dta

dta.coef <- remove_empty(dta$coefs, which = "cols", quiet = FALSE)
dta$fdnames$reps <- colnames(dta.coef)
dta.fd <- fd(dta.coef, dta$basis, dta$fdnames)
plot(dta.fd)

t.sq <- seq(0,100, length=501)
mean.female <- mean.fd(fdata2fd(fdata(ausfemale)))
stat <- Z(x=dta.fd, t.seq = t.sq, mu=mean.female)
########################################################################


########################################################################

# Generate data
set.seed(257)
J <- 20
N <- 80
epsilon <- matrix(NA, ncol=J, nrow=N+1)
for(i in 1:J) {
  gamma <- rexp(1, 0.5)
  epsilon[,i] <- rnorm(N+1, 0, gamma)
}
sigma <- sqrt(1.5)
t <- seq(0, 1, length = N+1)
#xj dimension NxJ
xj <- sin(pi*t)
y <- xj+sigma*epsilon
matplot(t, y, type="l")
lines(t, xj, col=2, lwd=2)

# Smooth data to functional generated series

timeRng   <- c(0,1)
nbasis    <- 51
genrbasis <- create.fourier.basis(timeRng,nbasis)

harmaccelLfd <- vec2Lfd(c(0, pi/2, pi), timeRng)
genrlambda <- 0.0000000001
genrfdPar  <- fdPar(genrbasis, harmaccelLfd, genrlambda)

genrfdSmooth <- smooth.basis(t, y, genrfdPar)
genrfdpen   <- genrfdSmooth$fd

plot(genrfdSmooth)
plotfit.fd(y, t, genrfdpen)

# Smooth to functionl data sin(pi*t) = xj

sinefun <- smooth.basis(t, xj, genrfdPar)
sinefunfd <- sinefun$fd

plot(sinefunfd)
plotfit.fd(xj, t, sinefunfd)

# Test hypothesis, that mean of generated curves is a sinefunfd
# H0: mu(curves) = xj
# H1: mu(curves) != xj
stat <- Z(x=genrfdpen, t.seq = t, mu=sinefunfd)

########################################################################
#### One sample pointwise-bootstrap-test                            ####
########################################################################

source("Zboottest.R")

# Test hypothesis, that mean of generated curves is a sinefunfd
# H0: mu(curves) = xj
# H1: mu(curves) != xj
stat <- Z.boot(x=genrfdpen, t.seq = t, mu=sinefunfd, replication = 500)

# H0: mu(Victoria) = mu(Australia)
# H1: mu(Victoria) != mu(Australia)
stat <- Z.boot(x=dta.fd, t.seq = t.sq, mu=mean.female)

# H0: mu(Tasmania) = mu(Australia)
# H1: mu(Tasmania) != mu(Australia)
stat <- Z.boot(x=dta.T.fd, t.seq = t.sq, mu=mean.male)


########################################################################
#### One sample L2_norm_based_test                                  ####
########################################################################

source("trace.R")
source("L2stat.R")

# Test hypothesis, that mean of generated curves is a sinefunfd
# H0: mu(curves) = xj
# H1: mu(curves) != xj
stat <- L2.stat(x=genrfdpen, t.seq = t, mu0=sinefunfd, replication = 500, method=1)
stat
stat <- L2.stat(x=genrfdpen, t.seq = t, mu0=sinefunfd, replication = 500, method=2)
stat

# H0: mu(Tasmania) = mu(Australia)
# H1: mu(Tasmania) != mu(Australia)
stat <- L2.stat(x=dta.T.fd, t.seq = t.sq, mu0=mean.male, replication = 500, method=1)
stat
stat <- L2.stat(x=dta.T.fd, t.seq = t.sq, mu0=mean.male, replication = 500, method=2)
stat

# H0: mu(Victoria) = mu(Australia)
# H1: mu(Victoria) != mu(Australia)
stat <- L2.stat(x=dta.fd, t.seq = t.sq, mu0=mean.female, replication = 500, method=1)
stat
stat <- L2.stat(x=dta.fd, t.seq = t.sq, mu0=mean.female, replication = 500, method=2)
stat


########################################################################
#### One sample F-type-test                                         ####
########################################################################

source("trace.R")
source("Fstat.R")

# Test hypothesis, that mean of generated curves is a sinefunfd
# H0: mu(curves) = xj
# H1: mu(curves) != xj
stat <- F.stat(x=genrfdpen, t.seq = t, mu0=sinefunfd, replication = 500, method=1)
stat
stat <- F.stat(x=genrfdpen, t.seq = t, mu0=sinefunfd, replication = 500, method=2)
stat

# H0: mu(Tasmania) = mu(Australia)
# H1: mu(Tasmania) != mu(Australia)
stat <- F.stat(x=dta.T.fd, t.seq = t.sq, mu0=mean.male, replication = 500, method=1)
stat
stat <- F.stat(x=dta.T.fd, t.seq = t.sq, mu0=mean.male, replication = 500, method=2)
stat

# H0: mu(Victoria) = mu(Australia)
# H1: mu(Victoria) != mu(Australia)
stat <- F.stat(x=dta.fd, t.seq = t.sq, mu0=mean.female, replication = 500, method=1)
stat
stat <- F.stat(x=dta.fd, t.seq = t.sq, mu0=mean.female, replication = 500, method=2)
stat

########################################################################
#### Two sample pointwise-test                                      ####
########################################################################

source("Ztwosample.R")


# Age-specific mortality rates for Australia and Australian states.
# Example: compare for male data in Australia with female mean
# H0: mu(male) = mu(female)
# H1: mu(male) != mu(female)
opar <- par(mfrow=c(2,1))
plot(ausmale)
plot(ausfemale)
par(opar)

dta.m <- fdata2fd(fdata(ausmale))
dta.f <- fdata2fd(fdata(ausfemale))
t.sq <- seq(0,100, length=501)

stat <- Ztwosample(x=dta.m, y=dta.f, t.seq = t.sq)
stat

# H0: mu(Tasmania) = mu(Australia)
# H1: mu(Tasmania) != mu(Australia)
stat <- Ztwosample(x=dta.T.fd, y = dta.m, t.seq = t.sq)
stat

# H0: mu(Victoria) = mu(Australia)
# H1: mu(Victoria) != mu(Australia)
stat <- Ztwosample(x=dta.fd, y = dta.f, t.seq = t.sq)
stat

########################################################################
#### Two sample L2-norm-based-test                                  ####
########################################################################

source("L2stattwosample.R")

# Age-specific mortality rates for Australia and Australian states.
# Example: compare for male data in Australia with female mean
# H0: mu(male) = mu(female)
# H1: mu(male) != mu(female)
stat <- L2.stat.twosample(x=dta.m, y=dta.f, t.seq = t.sq, method=1)
stat
stat <- L2.stat.twosample(x=dta.m, y=dta.f, t.seq = t.sq, method=2, replications=500)
stat$pvalue

# H0: mu(Tasmania) = mu(Australia)
# H1: mu(Tasmania) != mu(Australia)
stat <- L2.stat.twosample(x=dta.T.fd, y = dta.m, t.seq = t.sq, method=1)
stat
stat <- L2.stat.twosample(x=dta.T.fd, y = dta.m, t.seq = t.sq, replications = 500, method=2)
stat$pvalue

# H0: mu(Victoria) = mu(Australia)
# H1: mu(Victoria) != mu(Australia)
stat <- L2.stat.twosample(x=dta.fd, y = dta.f, t.seq = t.sq, method=1)
stat
stat <- L2.stat.twosample(x=dta.fd, y = dta.f, t.seq = t.sq, replications = 500, method=2)
stat$pvalue

########################################################################
#### Two sample F-type-test                                         ####
########################################################################

source("Fstattwosample.R")

# Age-specific mortality rates for Australia and Australian states.
# Example: compare for male data in Australia with female mean
# H0: mu(male) = mu(female)
# H1: mu(male) != mu(female)
stat <- F.stat.twosample(x=dta.m, y=dta.f, t.seq = t.sq, method=1)
stat
stat <- F.stat.twosample(x=dta.m, y=dta.f, t.seq = t.sq, method=2, replications=500)
stat$pvalue

# H0: mu(Tasmania) = mu(Australia)
# H1: mu(Tasmania) != mu(Australia)
stat <- F.stat.twosample(x=dta.T.fd, y = dta.m, t.seq = t.sq, method=1)
stat
stat <- F.stat.twosample(x=dta.T.fd, y = dta.m, t.seq = t.sq, replications = 500, method=2)
stat$pvalue

# H0: mu(Victoria) = mu(Australia)
# H1: mu(Victoria) != mu(Australia)
stat <- F.stat.twosample(x=dta.fd, y = dta.f, t.seq = t.sq, method=1)
stat
stat <- F.stat.twosample(x=dta.fd, y = dta.f, t.seq = t.sq, replications = 500, method=2)
stat$pvalue

########################################################################
#### Two sample permutation test                                    ####
########################################################################

#fda library

?tperm.fd

stat <- tperm.fd(dta.m, dta.f)
stat

stat <- tperm.fd(dta.T.fd, dta.m)
stat

stat <- tperm.fd(dta.fd, dta.f)
stat

########################################################################
#### Berkely growth data                                            ####
########################################################################

# This tests the difference between boys and girls heights in the
# Berkeley growth data.

# First set up a basis system to hold the smooths

knots  <- growth$age
norder <- 6
nbasis <- length(knots) + norder - 2
hgtbasis <- create.bspline.basis(range(knots), nbasis, norder, knots)

# Now smooth with a fourth-derivative penalty and a very small smoothing
# parameter

Lfdobj <- 4
lambda <- 1e-2
growfdPar <- fdPar(hgtbasis, Lfdobj, lambda)

hgtmfd <- smooth.basis(growth$age, growth$hgtm, growfdPar)$fd
hgtffd <- smooth.basis(growth$age, growth$hgtf, growfdPar)$fd

# Call tperm.fd

tres <- tperm.fd(hgtmfd,hgtffd)
tres


tsq <- seq(1,18, length=601)
tres <- Ztwosample(x=hgtmfd, y = hgtffd, t.seq = tsq)
tres


tres <- L2.stat.twosample(x=hgtmfd, y = hgtffd, t.seq = tsq, method=1)
tres
tres <- L2.stat.twosample(x=hgtmfd, y = hgtffd, t.seq = tsq, replications = 500, method=2)
tres$pvalue

tres <- F.stat.twosample(x=hgtmfd, y = hgtffd, t.seq = tsq, method=1)
tres
tres <- F.stat.twosample(x=hgtmfd, y = hgtffd, t.seq = tsq, replications = 500, method=2)
tres$pvalue
