library(VGAM)
library(aod)
library(dplyr)
options(scipen=999)

path <- '/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/romans/'
count = read.table(paste0(path, 'CLOCK/data/count.txt', sep=''))
N = as.numeric(colSums(count))
tissue=c(rep('CSF', 4),rep('PBMC', 4))
states = c('BSE', 'YR1', 'BSE', 'YR1', 'BSE', 'YR1', 'BSE', 'YR1')

BetaBinomialReg = function(count,N,C){
  res = lapply(c(1:length(count[,1])),function(i){
    x = as.numeric(count[i,])
    data = data.frame(x=x,N=N,C=as.factor(C))
    fm1 <- betabin(cbind(x, N - x) ~ C, ~ 1, data)
    res = summary(fm1)@Coef
    p = res[2, 4]
    fc = mean(x[C==T]/N[C==T])/mean(x[C==F]/N[C==F]) 
    return(c(p,fc))
  })
  res = do.call(cbind,res)
  return(res)
}

# MS vs HC (CSF only)
res1 = BetaBinomialReg(count[,tissue=='CSF'], N[tissue=='CSF'], 
                       states[tissue=='CSF']=='YR1')
write.csv(res1, paste0(path, 'CLOCK/data/binom_CSF.csv'), row.names=FALSE)
# MS vs HC (PBMC only)
res2 = BetaBinomialReg(count[,tissue=='PBMC'], N[tissue=='PBMC'], 
                       states[tissue=='PBMC']=='YR1')
write.csv(res2, paste0(path, 'CLOCK/data/binom_PBMC.csv'), row.names=FALSE)

