#### auxiliary functions ####
setNA <- function (v, a) {
  v[is.na(v)] <- a; v
}

# loglikelihood function
# dat_mx should have three columns containing input copies,
# number of positive reactions and number of negative reactions
logL <- function(dat_mx, p){
  prob_fail <- exp(-dat_mx[,1] * p)
  sum(log(prob_fail) * dat_mx[,3] +
        log(1 - prob_fail) * dat_mx[, 2])
}

#### LoD for a single stage ####

LoD <- function(dat_mx, level=0.95){

  # maximise likelihood estimation of p, the probability that
  # a single molecule can be positively tested
  func <- function(p) logL(dat_mx, p)
  ret <- optimise(func, c(0,1), maximum=TRUE)
  p <- ret$maximum

  # confidence interval for p constructed using log likelihood
  # ratio test (chisq with 1 dof under the null)
  ps <- seq(0, 1, by=0.0001)
  loglik <- sapply(ps, func)
  threshold <- ret$objective - qchisq(level, df=1) / 2
  CI_p <- range(ps[setNA(loglik > threshold, FALSE)])

  plot(ps, loglik, type='l') # visualise log likelihood
  abline(v=p, lty=1, col='red') # MLE
  abline(v=CI_p, lty=2, col='red') # plot confidence interval for p

  # use MLE and CI of p to compute MLE and CI of LoD-95
  LoD <- round(-log(1-level)/p, 1)
  CI_LoD <- round(rev(-log(1-level)/CI_p), 1)

  return(list(prob=p, CI_prob=CI_p, LoD=LoD, CI_LoD=CI_LoD))
}


#### Two stage LoD ####

LoD_ts <- function(mx1, mx2, level=0.95){
  alpha <- 1 - level
  ratio <- sum(mx1[,-1]) / sum(mx2[,-1])
  alpha2 <- alpha / (ratio + 1)
  alpha1 <- alpha2 * ratio
  ret1 <- LoD(mx1, level=1-alpha1)
  ret2 <- LoD(mx2, level=1-alpha2)
  p <- ret1$prob; CI_p <- ret1$CI_prob
  q <- ret2$prob; CI_q <- ret2$CI_prob

  # compute combined detection probability and LoD
  s <- 1-(1-p)*(1-q)
  CI_s <- 1-(1-CI_p)*(1-CI_q)

  # compute combined LoD
  LoD_s <- round(-log(0.05)/s, 1)
  CI_LoD_s <- round(rev(-log(0.05)/CI_s), 1)

  return(list(prob=s, CI_prob=CI_s, LoD=LoD_s, CI_LoD=CI_LoD_s))
}

#### main ####
# Stage 1 result data matrix
mx1 <- matrix(c(10, 4, 8,
               50, 171, 5,
               100, 168, 1), ncol=3, byrow=TRUE)
colnames(mx1) <- c('copies', 'success', 'fail')
mx1

# Stage 2 result data matrix
mx2 <- matrix(c(10, 0, 4,
                50, 1, 0), ncol=3, byrow=TRUE)
colnames(mx2) <- c('copies', 'success', 'fail')
mx2

# Stage 1 LoD95
LoD(mx1)

# Stage 1 + 2 LoD95
LoD_ts(mx1, mx2)
