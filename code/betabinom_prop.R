# library(glmmTMB)
library(aod)
library(dplyr)
options(scipen=999)
library("bbmle") 

# 1) prep the dataframe
path <- '/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/'
count = read.table(paste0(path, 'CLOCK/output/output/scRNA/beta_binomial/count_CSF_W5_10.txt', sep=''))
N = as.numeric(colSums(count))
# tissue=c(rep('CSF', 4),rep('PBMC', 4))
tissue=c(rep('CSF', 8))
# always CHECK THE SHIFT ON ONE
states = c('BSE', 'BSE', 'BSE', 'BSE', 'BSE', 'W5_10', 'BSE', 'W5_10')
# CSF: W5_10:
# PBMC: W5_10: ('BSE', 'BSE', 'BSE', 'W5_10', 'BSE', 'BSE', 'BSE', 'W5_10', 'BSE', 'BSE', 'W5_10')
# PBMC: YR1: ('BSE', 'BSE', 'YR1', 'BSE', 'YR1', 'BSE', 'BSE', 'BSE', 'BSE', 'BSE')
# PBMC: YR2: ('BSE', 'YR2', 'YR2', 'YR2', 'YR2', 'BSE', 'BSE', 'YR2', 'YR2', 'BSE', 'YR2', 'BSE', 'YR2', 'BSE', 'YR2', 'BSE', 'YR2', 'BSE')

# 2) run beta-binomial test 
BetaBinomialReg = function(count,N,C){
  res = lapply(c(1:length(count[,1])),function(i){
    x = as.numeric(count[i,])                          # add constant + 5
    data = data.frame(x=x,N=N,C=as.factor(C))
    # fm1 <- glmmTMB(cbind(x, N - x) ~ C, family=betabinomial, data=data)
    # res = summary(fm1)$coefficients
    # p = res$cond[2, 4]
    ## Beta-binomial model
    fm1 <- betabin(cbind(x, N - x) ~ C, ~ 1, data)
    res = summary(fm1)@Coef
    p = res[2, 4]
    fc = mean(x[C==T]/N[C==T])/mean(x[C==F]/N[C==F]) 
    return(c(p,fc))
  })
  res = do.call(cbind,res)
  return(res)
}

# 3) save results
# CSF: BSE vs TREATMENT 
res1 = BetaBinomialReg(count[,tissue=='CSF'], N[tissue=='CSF'], 
                       states[tissue=='CSF']=='W5_10')
write.csv(res1, paste0(path, 'CLOCK/output/output/scRNA/beta_binomial/binom_CSF_W5_10.csv'), row.names=FALSE)
# PBMC: BSE vs TREATMENT
res2 = BetaBinomialReg(count[,tissue=='PBMC'], N[tissue=='PBMC'], 
                       states[tissue=='PBMC']=='YR2')
write.csv(res2, paste0(path, 'CLOCK/output/output/scRNA/beta_binomial/binom_PBMC_YR2.csv'), row.names=FALSE)


# 4.1) replace NaNs manually if appeared:
x = as.numeric(count[4,tissue=='CSF'])
N[tissue=='CSF']
C <- states[tissue=='CSF']=='W5_W10'
data = data.frame(x=x,N=N,C=as.factor(C))
betbin <- betabin(cbind(x, N - x) ~ C, ~ 1, data)
res3 = summary(betbin)@Coef
p2 <- res3[2, 4]
p2
res1[1, 4] <- 0.99
res1[2, 4] <- 1

# 4.2) replace NaNs manually if appeared:
x = as.numeric(count[25,tissue=='PBMC'])
N[tissue=='PBMC']
C <- states[tissue=='PBMC']=='YR1'
data = data.frame(x=x,N=N,C=as.factor(C))
glmm_betbin <- glmmTMB(cbind(x, N - x) ~ C,
                       family=betabinomial, data=data)
res3 <- summary(glmm_betbin)$coefficients
p2 <- res3$cond[2, 4]
p2
res2[1, 25] <- p2

# 5) adjust p-val
stats <- read.csv(paste0(path, 'CLOCK/output/output/scRNA/beta_binomial/not_adjusted/binom_CSF_W5_10.csv'))
stats <- stats[,-4]  # for CSF W5_10 only
stats <- as.data.frame(t(stats))
stats <- stats %>% 
  mutate(V1 = p.adjust(p = V1, method = "BH")) 
stats <- as.data.frame(t(stats))
write.csv(stats, paste0(path, 'CLOCK/output/output/scRNA/beta_binomial/adjusted/binom_CSF_W5_10.csv'), row.names=FALSE)







