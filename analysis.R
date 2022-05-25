library(fda)
library(fda.usc)
data("MontrealTemp")
library(ggplot2)
library(data.table)
library(funFEM)
library(magrittr)
library(janitor)
library(RCurl)
library(demography)
library(far)


dt <- data.table(MontrealTemp)
dt$year <- rownames(MontrealTemp)

dt <- melt.data.table(dt, id.vars = 'year')

ggplot2::ggplot(data = dt, 
                aes(x = variable, 
                    y = value, 
                    color = year, 
                    group = year)) + 
  geom_line() +
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
  theme(legend.position="right") +
  labs(x = "Month and day", y = "Temperature in celcius", title = 'Temperature in Celsius in Montreal 1961-1994')

dt[, `:=`(month = substr(variable, 0, 3))]
dt$month <- factor(x = dt$month, levels = c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec"))
dt_lookup <- data.table(month = c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec"),
                        month_num = c(1:12))
dt$month_num <- dt_lookup$month_num[match(dt$month, dt_lookup$month)]

dt_mth <- dt[, .(avg = mean(value)), by = .(year, month_num)]

ggplot2::ggplot(data = dt_mth, 
                aes(x = month_num, 
                    y = avg, 
                    color = year, 
                    group = year)) + 
  geom_line() +
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
  theme(legend.position="right") +
  labs(x = "Month", y = "Average temperature in celcius", title = 'Monthly average temperature in Celsius in Montreal 1961-1994')

ggplot(dt) + geom_boxplot(aes(x = month_num, y = value, group = month)) +
  labs(x = "Month", y = "Temperature in celcius", title = 'Temperature variation in Celsius in Montreal 1961-1994')

dt_mth$year %<>% as.numeric()
dt_mth_long <- dcast.data.table(data = dt_mth,
                                formula = month_num ~ year,
                                value.var = 'avg',
                                fun.aggregate = sum)

dt_mth_long_last5 <- dcast.data.table(data = dt_mth[year %in% seq(max(dt_mth$year)-4, max(dt_mth$year))],
                                      formula = month_num ~ year,
                                      value.var = 'avg',
                                      fun.aggregate = sum)
# data removing last 5 years
dt_mth_without_last5 <- dcast.data.table(data = dt_mth[year %in% seq(min(dt_mth$year), max(dt_mth$year)-4)],
                                      formula = month_num ~ year,
                                      value.var = 'avg',
                                      fun.aggregate = sum)
# data removing last 10 years
dt_mth_without_last10 <- dcast.data.table(data = dt_mth[year %in% seq(min(dt_mth$year), max(dt_mth$year)-9)],
                                         formula = month_num ~ year,
                                         value.var = 'avg',
                                         fun.aggregate = sum)

dtdt_mth_long_last10 <- dcast.data.table(data = dt_mth[year %in% seq(max(dt_mth$year)-9, max(dt_mth$year))],
                                       formula = month_num ~ year,
                                       value.var = 'avg',
                                       fun.aggregate = sum)


mt_temp_mth <- as.matrix(dt_mth_long[,-1])
# adding mt_temp_mth without last 5 years
mt_temp_mth_minus_5 <- as.matrix(dt_mth_without_last5[,-1])
# adding mt_temp_mth without last 5 years
mt_temp_mth_minus_10 <- as.matrix(dt_mth_without_last10[,-1])
mt_temp_mth_last5 <- as.matrix(dt_mth_long_last5[,-1])
mt_temp_mth_last10 <- as.matrix(dt_mth_long_last10[,-1])

####
#### Smoothing
####

mthRng  = c(0, 12)
mthbasis = create.fourier.basis(rangeval = mthRng, 
                                nbasis = 10)

# Define the harmonic acceleration operator

Lcoef        = c(0,(2*pi/diff(mthRng))^2,0)
harmaccelLfd = vec2Lfd(Lcoef, mthRng)

# Smooth data with lambda that minimizes GCV
lambda     = 0.013
fdParobj   = fdPar(mthbasis, harmaccelLfd, lambda)

# Total period
temp_mth_fd = smooth.basis(y = mt_temp_mth, 
                           fdParobj = fdParobj)
#Period without last 5 years
temp_mth_fd_without_5 = smooth.basis(y=mt_temp_mth_minus_5,
                                     fdParobj = fdParobj)
#Period without last 10 years
temp_mth_fd_without_10 = smooth.basis(y=mt_temp_mth_minus_10,
                                     fdParobj = fdParobj)

plot(temp_mth_fd, xlab = 'Month', ylab = 'Average monthly temperature')

# Last 5 years
temp_mth_last5_fd = smooth.basis(y = mt_temp_mth_last5, 
                                 fdParobj = fdParobj)
plot(temp_mth_last5_fd, xlab = 'Month', ylab = 'Average monthly temperature')

# Last 10 years
temp_mth_last10_fd = smooth.basis(y = mt_temp_mth_last10, 
                                  fdParobj = fdParobj)
plot(temp_mth_last10_fd, xlab = 'Month', ylab = 'Average monthly temperature')

# Elementary pointwise mean and standard deviation

meantemp   = mean.fd(temp_mth_fd$fd)
stddevtemp = std.fd(temp_mth_fd$fd)

lines(meantemp, lwd=4, lty=2, col=1)
lines(stddevtemp, lwd=4, lty=2, col=4)

#### Create fda object

temp_mth_fdata <- fdata(temp_mth_fd$fd)
plot(temp_mth_fdata)

#### Optimal smoothing

l <- seq(0.001,1,by = 0.001)
nb <- seq(4, 12, by = 2)

temp_mth_param <- optim.basis(temp_mth_fdata, 
                              numbasis = nb, 
                              lambda = l)
plot(temp_mth_param$fdata.est)
temp_mth_param$gcv.opt
temp_mth_param$lambda.opt
temp_mth_param$numbasis.opt

#### Derivatives

# First - velocity

temp_mth_d1 <- deriv.fd(temp_mth_fd$fd, 1)
plot(temp_mth_d1, xlab = "Month")

# Second - acceleration
temp_mth_d2 <- deriv.fd(temp_mth_fd$fd, 2)
plot(temp_mth_d2, xlab = "Month")

# phaseplanePlot(1, temp_mth_fd$fd)

#### Covariance

tempvarbifd <- var.fd(temp_mth_fd$fd)
# evaluate the variance-covariance surface and plot
mth <- seq(1,12,1)
tempvarmat <- eval.bifd(mth,mth,tempvarbifd)
# make a contour plot of the variance function
contour(mth, mth, tempvarmat,
        xlab="Month",
        ylab="Month",
        axes=TRUE)

####
#### Principal component analysis
####

nharm = 4
temp_mth_pcalist = pca.fd(temp_mth_fd$fd, nharm, centerfns = TRUE)

# attributes(temp_mth_pcalist)
# temp_mth_pcalist$harmonics

# plot(temp_mth_pcalist$harmonics)

plot.pca.fd(temp_mth_pcalist, harm = 1)
plot.pca.fd(temp_mth_pcalist, harm = 2)
plot.pca.fd(temp_mth_pcalist, harm = 3)
plot.pca.fd(temp_mth_pcalist, harm = 4)

# Rotate

temp_mth_rot_pcalist = varmx.pca.fd(temp_mth_pcalist)

plot.pca.fd(temp_mth_rot_pcalist, harm = 1)
plot.pca.fd(temp_mth_rot_pcalist, harm = 2)
plot.pca.fd(temp_mth_rot_pcalist, harm = 3)
plot.pca.fd(temp_mth_rot_pcalist, harm = 4)

# # Pca scores
# 
# rotpcascores = temp_mth_rot_pcalist$scores
# 
# plot(rotpcascores[,1], rotpcascores[,2], type="p", pch="o",
#      xlab="Rotated Harmonic I", ylab="Rotated Harmonic II")
# 
# namebank <- colnames(mt_temp_mth)
# text(rotpcascores[,1], rotpcascores[,2], labels=namebank, cex= 0.7, pos=3)
# 

####
#### Clustering
####

# weather data example
# temperature
# CanadianWeather_Temp <- CanadianWeather$dailyAv[,,"Temperature.C"]
# basis <- create.bspline.basis(c(0, 365), nbasis=21, norder=4) 
# fdobj <- smooth.basis(day.5, CanadianWeather_Temp, 
                      # basis, fdnames=list("Day", "Station", "Deg C"))$fd
res <- funFEM(temp_mth_fd$fd, 
              K=3)

par(mfrow=c(1,2))
plot.fd(temp_mth_fd$fd, col=res$cls, lwd=2, lty=1)
fdmeans <- temp_mth_fd$fd
fdmeans$coefs <- t(res$prms$my)
plot(fdmeans, col=1:max(res$cls), lwd=2)

####
#### Hypotheses testing
#### 
source("/custom-functions/Ztest.R")
source("/custom-functions/Zboottest.R")
source("/custom-functions/trace.R")
source("/custom-functions/L2stat.R")
source("/custom-functions/trace.R")
source("/custom-functions/Fstat.R")

#### Last 5 years

# H0: mu(Last 5 year average) = mu(Total average)
# H1: mu(Last 5 year average) != mu(Total average)

# NaN values in coefficient matrix
dta <- temp_mth_last5_fd$fd
dta

dta.coef <- remove_empty(dta$coefs, 
                         which = "cols", 
                         quiet = FALSE)
dta$fdnames$reps <- colnames(dta.coef)
dta.fd <- fd(dta.coef, 
             dta$basis, 
             dta$fdnames)
plot.fd(dta.fd, ylab = "Temperature")

t.sq <- seq(0,12, length = 61)
mean.total <- mean.fd(temp_mth_fd$fd)
mean_minus_5.total <- mean.fd(temp_mth_fd_without_5$fd)
mean_minus_10.total <- mean.fd(temp_mth_fd_without_10$fd)


# One sample pointwise t-test

par(mar=c(1,1,1,1))
stat <- Z(x = dta.fd, 
          t.seq = t.sq, 
          mu = mean_minus_5.total)

# One sample pointwise bootstrap t-test
par(mar=c(1,1,1,1))
stat_boot <- Z.boot(x = dta.fd, 
                    t.seq = t.sq, 
                    mu = mean_minus_5.total)

# One sample L2_norm_based_test 

stat_l2_1 <- L2.stat(x = dta.fd, 
                     t.seq = t.sq, 
                     mu0 = mean_minus_5.total, 
                     replication = 500, 
                     method = 1)
stat_l2_1$statistics

stat_l2_1$pvalue
stat_l2_1$params$df

stat_l2_1

stat_l2_2 <- L2.stat(x = dta.fd, 
                     t.seq = t.sq, 
                     mu0 = mean_minus_5.total, 
                     replication = 500, 
                     method = 2)
stat_l2_2$statistics
stat_l2_2$pvalue

# One sample F-type test 

stat_f1 <- F.stat(x = dta.fd, 
                  t.seq = t.sq, 
                  mu0 = mean_minus_5.total, 
                  replication = 500, 
                  method = 1)
stat_f1$statistics
stat_f1$pvalue
stat_f1$params$df1
stat_f1$params$df2

stat_f2 <- F.stat(x = dta.fd, 
                  t.seq = t.sq, 
                  mu0 = mean_minus_5.total, 
                  replication = 500, 
                  method=2)
stat_f2$statistics
stat_f2$pvalue

#### Last 10 years

# H0: mu(Last 10 year average) = mu(Total average)
# H1: mu(Last 10 year average) != mu(Total average)

# NaN values in coefficient matrix
dta <- temp_mth_last10_fd$fd
dta

dta.coef <- remove_empty(dta$coefs, 
                         which = "cols", 
                         quiet = FALSE)
dta$fdnames$reps <- colnames(dta.coef)
dta.fd <- fd(dta.coef, 
             dta$basis, 
             dta$fdnames)
plot(dta.fd)

t.sq <- seq(0,12, length = 61)
# mean.total <- mean_minus_10.fd(temp_mth_fd$fd)

# One sample pointwise t-test

stat <- Z(x = dta.fd, 
          t.seq = t.sq, 
          mu = mean_minus_10.total)

# One sample pointwise bootstrap t-test

stat_boot <- Z.boot(x = dta.fd, 
                    t.seq = t.sq, 
                    mu = mean_minus_10.total)

# One sample L2_norm_based_test 

stat_l2_1 <- L2.stat(x = dta.fd, 
                     t.seq = t.sq, 
                     mu0 = mean_minus_10.total, 
                     replication = 500, 
                     method = 1)
stat_l2_1$statistics
stat_l2_1$pvalue
stat_l2_1$params$df

stat_l2_2 <- L2.stat(x = dta.fd, 
                     t.seq = t.sq, 
                     mu0 = mean_minus_10.total, 
                     replication = 500, 
                     method = 2)
stat_l2_2$statistics
stat_l2_2$pvalue
stat_l2_2$params$df

# One sample F-type test 

stat_f1 <- F.stat(x = dta.fd, 
                  t.seq = t.sq, 
                  mu0 = mean_minus_10.total, 
                  replication = 500, 
                  method = 1)
stat_f1$statistics
stat_f1$pvalue
stat_f1$params$df1
stat_f1$params$df2

stat_f2 <- F.stat(x = dta.fd, 
                  t.seq = t.sq, 
                  mu0 = mean_minus_10.total, 
                  replication = 500, 
                  method=2)
stat_f2$statistics
stat_f2$pvalue

####
#### Useful notes
####

# Fourier basis, typically used for periodic data

# Before you can convert raw discrete data into a functional data object with
# these functions, you must specify a basis. A basis is a system of primitive
# functions that are combined linearly to approximate actual functions. An
# example is successive powers of an argument t, linear combinations of which
# form polynomials. A more useful example is the unit function 1 and successive pairs of sine and cosine functions with frequencies that are integer
# multiples of a base period that make up a Fourier series

####
#### Functional time series: Electricity example
####

# Such curves cannot be treated as being independent with the same distribution; they do not form a simple random sample of curves
# Temperature on a given day will be affected by temperature on previous day, so the curves will be dependent

dt_mth %<>% data.table()
dt_mth <- dt_mth[order(year, month_num)]
ts.monthly <- ts(dt_mth$avg, start = 1961, frequency = 12)
ts.monthly

ftemp <- far::as.fdata(ts.monthly, 
                       col = 1, 
                       p = 12, 
                       dates = 1961:1994, 
                       name = "Average temperature")
ftemp

multplot(ftemp$"Average temperature", 
         type = "l")

model1 <- far(data=ftemp, 
              y="Average temperature", 
              center=TRUE, 
              na.rm=FALSE)
coef.far(model1)
plot(model1)
names(model1)

source('/custom-functions/predict_far2.R')
pred <- predict.far2(model1, newdata = ftemp, label = 1962:1995)
multplot(pred$"Average temperature", type = "l")

pred1 <- predict.far2(model1, newdata = select.fdata(ftemp,date=1:33), label = 1962:1994)
# Real values
real1 <- select.fdata(ftemp,date=2:34)

errors1 <- pred1[[1]]-real1[[1]]
matplot(errors1, type = "l")
summary(errors1)
apply(errors1, 2, function(x) sqrt(sum(x^2)))

model2.cv <- far.cv(data=ftemp, y="Average temperature", ncv=2,
                    cvcrit = "Average temperature",
                    center=TRUE, na.rm = FALSE)
print(model2.cv)
k2 <- model2.cv$minL2[1]

model2 <- far(data=ftemp, y="Average temperature", kn=k2,
              center=TRUE, na.rm=FALSE)
print(model2)
coef.far(model2)
plot(model2)

pred2 <- predict.far2(model2, newdata = ftemp, label = 1962:1995)
multplot(pred$"Average temperature", type = "l")

pred2 <- predict.far2(model2, newdata = select.fdata(ftemp,date=1:33), label = 1962:1994)

errors2 <- pred2[[1]]-real1[[1]]
matplot(errors2, type = "l")
summary(errors2)
apply(errors2, 2, function(x) sqrt(sum(x^2)))
apply(errors1, 2, function(x) sqrt(sum(x^2)))

###
basis <- create.polygonal.basis((0:11)/12)
mtemp <- matrix(ts.monthly, ncol=34)
colnames(mtemp) <- 1961:1994
ftemp2 <- smooth.basis((0:11)/12, mtemp, basis)
plotfit.fd(mtemp, (0:11)/12, ftemp2$fd)

ele.pc <- pca.fd(ftemp2$fd, 1)
plot(ele.pc)

# dev.new()
plot(model1)
# dev.new()
plot(ele.pc$harmonics)

vj <- ele.pc$harmonics$coefs/sqrt(sum(ele.pc$harmonics$coefs^2))
data.frame(model1$v[[1]], vj)
