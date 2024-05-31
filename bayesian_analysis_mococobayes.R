.libPaths(c("/home/mococobayes/R/x86_64-pc-linux-gnu-library/4.4",.libPaths()))

library(rstan) 
library(brms)
library(modelr)

# Pour optimiser Stan:
Sys.setenv(LC_ALL="en_US.UTF-8")
CORES = parallel::detectCores()
options(mc.cores = CORES)

home = getwd()

# loading data 
data2 <- read.csv("data2.csv")
dataOcculted <- read.csv("dataOcculted.csv")
dataJ <- read.csv("dataJ.csv")
dataE <- read.csv("dataE.csv")

## Second experiment
data2$segment     <- as.factor(data2$segment)
data2$repeated       <- as.factor(data2$repeated)

## Joystick data
dataJ$segment       <- as.factor(dataJ$segment)
dataJ$day           <- as.factor(dataJ$day)
dataJ$repeated       <- as.factor(dataJ$repeated)

## Eyetracking data
dataE$segment       <- as.factor(dataE$segment)
dataE$day           <- as.factor(dataE$day)
dataE$repeated       <- as.factor(dataE$repeated)

## Occulted data
dataOcculted$repeatedSegment     <- as.factor(dataOcculted$repeatedSegment)
dataOcculted$occulted            <- as.factor(dataOcculted$occulted)
dataOcculted$subjects            <- as.factor(dataOcculted$subjects)

# The participant 4 didn't do the task correctly
dataOcculted$MI[dataOcculted$subjects==4] <- NaN
dataOcculted$xcorr[dataOcculted$subjects==4] <- NaN
dataOcculted$RMSE[dataOcculted$subjects==4] <- NaN
dataOcculted$spectral_diff[dataOcculted$subjects==4] <- NaN

# priors
pr <- set_prior("normal(0, 30)", class = "b")
prst = set_prior("student_t(1, 0, 0.1)", class = "b")

## First experiment analysis
print("_________  First experiment analysis  ____________")
fJRMSE <- bf(RMSE ~ segment*day + segmentLength + (1 + segmentLength + segment*day  | subjects))
print('RMSE ~ segment*day + segmentLength + (1 + segmentLength + segment*day  | subjects)')
mJRMSE <- brm(fJRMSE , data = dataJ, prior = pr, sample_prior = TRUE, save_pars = save_pars(all=TRUE), iter=5000)
saveRDS(mJRMSE, "./savings/mJRMSE_good.RDS")

fERMSE <- bf(RMSE ~ segment*day + segmentLength + (1 + segmentLength + segment*day  | subjects))
print('RMSE ~ segment*day + segmentLength + (1 + segmentLength + segment*day  | subjects)')
mERMSE <- brm(fERMSE , data = dataE, prior = pr, sample_prior = TRUE, control=list(adapt_delta=0.99), save_pars = save_pars(all=TRUE), iter=5000)
saveRDS(mERMSE, "./savings/mERMSE_good.RDS")

fJpupil <- bf(pupil_0.057 ~ segment*day + segmentLength + (1 + segmentLength + segment*day  | subjects))
print('pupil_0.057 ~ segment*day + segmentLength + (1 + segmentLength + segment*day  | subjects)')
mJpupil <- brm(fJpupil , data = dataJ, prior = prst, sample_prior = TRUE, control=list(adapt_delta=0.95), save_pars = save_pars(all=TRUE), iter=5000)
saveRDS(mJpupil, "./savings/mJpupil_good.RDS")

fJFF <- bf(FF ~ segment*day + segmentLength + (1 + segmentLength + segment*day  | subjects))
print('FF ~ segment*day + segmentLength + (1 + segmentLength + segment*day  | subjects)')
mJFF <- brm(fJFF , data = dataJ, prior = pr, sample_prior = TRUE, save_pars = save_pars(all=TRUE), iter=5000)
saveRDS(mJFF, "./savings/mJFF_good.RDS")

fEFF <- bf(FF ~ segment*day + segmentLength + (1 + segmentLength + segment*day  | subjects))
print('FF ~ segment*day + segmentLength + (1 + segmentLength + segment*day  | subjects)')
mEFF <- brm(fEFF , data = dataE, prior = pr, sample_prior = TRUE, control=list(adapt_delta=0.95), save_pars = save_pars(all=TRUE), iter=5000)
saveRDS(mEFF, "./savings/mEFF_good.RDS")

## Mediation analysis
# On Joystick
fJFF_repeated <- bf(FF ~ repeated*day + segmentLength + (1 + segmentLength + repeated*day  | subjects))
fJRMSE_FF <- bf(RMSE ~ repeated*day + FF + segmentLength + (1 + FF + segmentLength + repeated*day  | subjects))
print('mediation : fJFF_repeated + fJRMSE_FF + set_rescor(FALSE))')
mJRMSEmediation <- brm(fJFF_repeated + fJRMSE_FF + set_rescor(FALSE), data = dataJ, refresh = 0)
saveRDS(mJRMSEmediation, "./savings/mJRMSEmediation_good.RDS")
print('mediation data saved')

# On eyetracking
fEFF_repeated <- bf(FF ~ repeated*day + segmentLength + (1 + segmentLength + repeated*day  | subjects))
fERMSE_FF <- bf(RMSE ~ repeated*day + FF + segmentLength + (1 + FF + segmentLength + repeated*day  | subjects))
print('mediation : fEFF_repeated + fERMSE_FF + set_rescor(FALSE))')
mERMSEmediation <- brm(fEFF_repeated + fERMSE_FF + set_rescor(FALSE), data = dataE, refresh = 0)
saveRDS(mERMSEmediation, "./savings/mERMSEmediation_good.RDS")
print('mediation data saved')

## Second experiment analysis
print("_________  Second experiment analysis  ____________")
f2RMSE <- bf(RMSE ~ segment + segment_length + (1 + segment_length + segment  | subjects))
print('RMSE ~ segment + segment_length + (1 + segment_length + segment  | subjects)')
m2RMSE <- brm(f2RMSE , data = data2, prior = pr, sample_prior = TRUE, save_pars = save_pars(all=TRUE), iter=5000)
saveRDS(m2RMSE, "./savings/m2RMSE_good.RDS")

f2FF <- bf(FF ~ segment + segment_length + (1 + segment_length + segment  | subjects))
print('FF ~ segment + segment_length + (1 + segment_length + segment  | subjects)')
m2FF <- brm(f2FF , data = data2, prior = pr, sample_prior = TRUE, save_pars = save_pars(all=TRUE), iter=5000)
saveRDS(m2FF, "./savings/m2FF_good.RDS")

f2pupil <- bf(pupil_uncorrected_1 ~ segment + segment_length + (1 + segment_length + segment  | subjects))
print('pupil_uncorrected_1 ~ segment + segment_length + (1 + segment_length + segment  | subjects)')
m2pupil <- brm(f2pupil , data = data2, prior = prst, sample_prior = TRUE, control=list(adapt_delta=0.95), save_pars = save_pars(all=TRUE), iter=5000)
saveRDS(m2pupil, "./savings/m2pupil_good.RDS")

# mediation analysis
f2FF_repeated <- bf(FF ~ repeated + segment_length + (1 + segment_length + repeated  | subjects))
f2RMSE_FF     <- bf(RMSE ~ repeated + FF + segment_length + (1 + segment_length + repeated  | subjects))
print('mediation alaysis : f2FF_repeated + f2RMSE_FF + set_rescor(FALSE)')
m2RMSEmediation <- brm(f2FF_repeated + f2RMSE_FF + set_rescor(FALSE), data = data2, refresh = 0)
saveRDS(m2RMSEmediation, "./savings/m2RMSEmediation_good.RDS")


## Occulatation analysis
print("_________  Occulatation analysis  ____________")
fOccXCORR <- bf(xcorr ~ occulted*repeatedSegment + segmentLength + (1 + segmentLength + occulted*repeatedSegment  | subjects))
print('xcorr ~ occulted*repeatedSegment + segmentLength + (1 + segmentLength + occulted*repeatedSegment  | subjects)')
mOccXCORR <- brm(fOccXCORR , data = dataOcculted, prior = pr, sample_prior = TRUE, save_pars = save_pars(all=TRUE), iter=5000)
saveRDS(mOccXCORR, "./savings/m2_occult_XCORR_good.RDS")

print('job finished')