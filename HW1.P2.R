# Problem 2
# change alpha with 0.05 and 0.01 here
alpha <- 0.01
p0 <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
n0 <- c(10, 25, 50, 100, 250)

rej.lrt <- matrix(0, nrow = length(n0), ncol = length(p0))
rej.sc <- matrix(0, nrow = length(n0), ncol = length(p0))
rej.wl <- matrix(0, nrow = length(n0), ncol = length(p0))

# run three tests and compute type 1 error
for (i in 1:length(n0)) {
  for (j in 1:length(p0)) {
    set.seed(123)
    y0 <- rbinom(1000, n0[i], p0[j]) 
    for (k in 1:length(y0)) {
      
      # loglikelihood and finding MLE
      loglik <- function(p, y, n) log(choose(n, y)) + y*log(p) +(n-y)*log(1-p)
      negloglik <- function (p, y, n) - loglik(p, y, n)
      out.nlm <- nlm(negloglik, 0.01, y = y0[k], n = n0[i], hessian = TRUE)
      phat <- out.nlm$estimate
      se <- sqrt(1/out.nlm$hessian)
      
      # LRT for H0: p=p0, and compute p-value
      lrt.stat <- 2*(loglik(phat, y0[k], n0[i]) - loglik(p0[j], y0[k], n0[i]))
      pval.lrt <- 1 - pchisq(lrt.stat, df=1)
      rej1 <- ifelse(pval.lrt < alpha, 1, 0)
      rej.lrt[i,j] <- rej1 + rej.lrt[i,j]
      
      # score test and compute p-value
      score <- function(p, y, n) y/p - (n-y)/(1-p)
      info <- function(p, n) n/p/(1-p)
      sc.stat <- score(p0[j], y0[k], n0[i])/sqrt(info(p0[j], n0[i]))
      pval.sc <- 2*(1-pnorm(abs(sc.stat)))
      # 1 - pchisq(sc.stat^2, df=1)  (should be equal to pval.sc)
      rej2 <- ifelse(pval.sc < alpha, 1, 0)
      rej.sc[i,j] <- rej2 + rej.sc[i,j]
      
      # wald test and compute p-value
      z <- (phat-p0[j])/se
      pval.wl <- 2*(1-pnorm(abs(z)))
      rej3 <- ifelse(pval.wl < alpha, 1, 0)
      rej.wl[i,j] <- rej3 + rej.wl[i,j]
    }
  }
}

############################## result analysis I ##############################
# result analysis I: sum up all errors of all p0 and n0
numrej.lrt <- sum(rej.lrt)
numrej.sc <- sum(rej.sc)
numrej.wl <- sum(rej.wl)

# when alpha = 0.05, run this code only
perc.lrt1 <- (numrej.lrt)/(6*5*1000)
perc.sc1 <- (numrej.sc)/(6*5*1000)
perc.wl1 <- (numrej.wl)/(6*5*1000)
alpha1.rej1 <- round(c(0.05, perc.lrt1, perc.sc1, perc.wl1), 5)
print(alpha1.rej1)

# when alpha = 0.01, skip part above and run code here
perc.lrt2 <- (numrej.lrt)/(6*5*1000)
perc.sc2 <- (numrej.sc)/(6*5*1000)
perc.wl2 <- (numrej.wl)/(6*5*1000)
alpha2.rej2 <- round(c(0.01, perc.lrt2, perc.sc2, perc.wl2), 5)
print(alpha2.rej2)

# summarize above results
data.frame(alpha1.rej1, alpha2.rej2, row.names = c("alpha", "lrt", "sc", "wl"))

############################## result analysis II ##############################
# result analysis II: sum up errors of all p0 for each n0
library(ggplot2)

# when alpha = 0.05
# summarize result of errors
err1 <- matrix(0, nrow = 5, ncol = 4)
for(i in 1:length(n0)) {
  err1[i,1] <- n0[i]
  err1[i,2] <- (rowSums(rej.lrt)/(6*1000))[i]
  err1[i,3] <- (rowSums(rej.sc)/(6*1000))[i]
  err1[i,4] <- (rowSums(rej.wl)/(6*1000))[i]
}
colnames(err1) <- c("n0", "lrt", "sc", "wl")
err1.df <- data.frame(err1)

# plot n0 vs error with alpha=0.05
ggplot(err1.df, aes(x=n0)) +
  geom_point(aes(y=lrt), color = "blue3") +
  geom_point(aes(y=sc), color = "violetred3") + 
  geom_point(aes(y=wl), color = "gold2") +
  geom_line(aes(y=lrt), linetype = 2, color = "blue3") +
  geom_line(aes(y=sc), linetype = 6, color = "violetred3") + 
  geom_line(aes(y=wl), color = "gold2") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", size = 1) +
  labs(x = "n", y = "type 1 error", title = "sample size vs type 1 error (alpha=0.05)") + 
  scale_colour_manual(name = "tests", 
                      labels = c("lrt", "sc", "wl"),
                      values = c("lrt"="blue3", 
                                 "sc"="violetred3",
                                 "wl"="gold2")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18))


# when alpha = 0.01
# summarize result of errors
err2 <- matrix(0, nrow = 5, ncol = 4)
for(i in 1:length(n0)) {
  err2[i,1] <- n0[i]
  err2[i,2] <- (rowSums(rej.lrt)/(6*1000))[i]
  err2[i,3] <- (rowSums(rej.sc)/(6*1000))[i]
  err2[i,4] <- (rowSums(rej.wl)/(6*1000))[i]
}
colnames(err2) <- c("n0", "lrt", "sc", "wl")
err2.df <- data.frame(err2)

# plot n0 vs error with alpha=0.01
ggplot(err2.df, aes(x=n0)) +
  geom_point(aes(y=lrt), color = "blue3") +
  geom_point(aes(y=sc), color = "violetred3") + 
  geom_point(aes(y=wl), color = "gold2") +
  geom_line(aes(y=lrt), linetype = 2, color = "blue3") +
  geom_line(aes(y=sc), linetype = 6, color = "violetred3") + 
  geom_line(aes(y=wl), color = "gold2") +
  geom_hline(yintercept = 0.01, linetype = "dashed", color = "red", size = 1) +
  labs(x = "n", y = "type 1 error", title = "sample size vs type 1 error (alpha=0.01)") + 
  scale_colour_manual(name = "tests", 
                      labels = c("lrt", "sc", "wl"),
                      values = c("blue3", 
                                 "violetred3",
                                 "gold2")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18)) 
  

