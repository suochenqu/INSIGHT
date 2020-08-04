# data matrix, with three columns containing input copies,
# number of positive reactions and number of negative reactions
dat_mx <- matrix(c(10, 4, 8,
                   50, 11, 1 ,
                   100, 6, 0), ncol=3, byrow=TRUE)
colnames(dat_mx) <- c('copies', 'success', 'fail')
dat_mx

# loglikelihood function
loglikelihood <- function(p){
  prob_fail <- exp(-dat_mx[,1] * p)
  sum(log(prob_fail) * dat_mx[,'fail'] + log(1 - prob_fail) * dat_mx[, 'success'])
}

# maximise likelihood estimation of p, the probability that a single molecule can be amplified
ret <- optimise(loglikelihood, c(0,1), maximum=TRUE)
p <- ret$maximum

ps = seq(0,1,by=0.0001)
loglik <- sapply(ps, loglikelihood)
setNA <- function (v, a) {v[is.na(v)] <- a; v}
# confidence interval for p constructed using log likelihood ratio test (chisq with 1 dof under the null)
CI_p <- range(ps[setNA(loglik > ret$objective - qchisq(0.95, df=1) / 2, FALSE)]) 

plot(ps, loglik, type='l') # visualise log likelihood
abline(v=p, lty=1, col='red') # MLE
abline(v=CI_p, lty=2, col='red') # plot confidence interval for p

# use MLE and CI of p to compute MLE and CI of LoD-95
LoD <- round(-log(0.05)/p, 1)
CI_LoD <- round(rev(-log(0.05)/CI_p), 1)

cat('LoD95 = ', LoD, ', CI = (', CI_LoD[1], ', ', CI_LoD[2], ')\n', sep='')
