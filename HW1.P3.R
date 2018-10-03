# Problem 3
install.packages("binom")
install.packages("ggpubr")
library(binom)

# change alpha with 0.05 and 0.01 here
alpha <- 0.01
p0 <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
n0 <- c(10, 25, 50, 100, 250)

coverage.lrt <- matrix(0, nrow = length(p0), ncol = length(n0))
coverage.sc <- matrix(0, nrow = length(p0), ncol = length(n0))
coverage.wl <- matrix(0, nrow = length(p0), ncol = length(n0))

# create function to compute wald test CI
prop.wald <- function(y, n, alpha) {
  phat <- y/n
  se <- sqrt(phat*(1-phat)/n)
  z <- qnorm(c(alpha/2, 1-alpha/2))
  ci <- phat + z*se
  return(ci)
}

# compute three tests' CI
for (i in 1:length(p0)) {
  for (j in 1:length(n0)) {
    set.seed(123)
    y0 <- rbinom(1000, n0[j], p0[i]) 
    # create a matrix cover that has 1000 simulations rows and 3 tests columns
    cover <- matrix(0, 1000, 3)
    
    for (k in 1:length(y0)) {
      
      # LRT CI
      l <- binom.lrt(y0[k], n0[j], conf.level = 1-alpha)[5]
      u <- binom.lrt(y0[k], n0[j], conf.level = 1-alpha)[6]
      ci.lrt <- c(l, u)
      cover[k,1] <- (ci.lrt[1] < p0[i]) & (p0[i] < ci.lrt[2])
      
      # score test CI
      ci.sc <- prop.test(y0[k], n0[j], conf.level = 1-alpha, correct = TRUE)$conf
      cover[k,2] <- (ci.sc[1] < p0[i]) & (p0[i] < ci.sc[2]) 
      
      # wald test CI
      ci.wl <- prop.wald(y0[k], n0[j], alpha)
      cover[k,3] <- (ci.wl[1] < p0[i]) & (p0[i] < ci.wl[2])
      
    }
   coverage.lrt[i,j] <- apply(cover, 2, sum)[1]
   coverage.sc[i,j] <- apply(cover, 2, sum)[2]
   coverage.wl[i,j] <- apply(cover, 2, sum)[3]
  } 
}

############################## result analysis ##############################
# create matrix that has p0 and 3 tests under different n0
coverage.n1 <- matrix(0, nrow = length(p0), 4, 
                      dimnames = list(c("1","2","3","4","5","6"),
                                      c("p0", "lrt", "sc", "wl")))
coverage.n2 <- matrix(0, nrow = length(p0), 4, 
                      dimnames = list(c("1","2","3","4","5","6"),
                                      c("p0", "lrt", "sc", "wl")))
coverage.n3 <- matrix(0, nrow = length(p0), 4, 
                      dimnames = list(c("1","2","3","4","5","6"),
                                      c("p0", "lrt", "sc", "wl")))
coverage.n4 <- matrix(0, nrow = length(p0), 4, 
                      dimnames = list(c("1","2","3","4","5","6"),
                                      c("p0", "lrt", "sc", "wl")))
coverage.n5 <- matrix(0, nrow = length(p0), 4, 
                      dimnames = list(c("1","2","3","4","5","6"),
                                      c("p0", "lrt", "sc", "wl")))

# n=10
for(i in 1:length(p0)) {
  coverage.n1[i,1] <- p0[i]
  coverage.n1[i,2] <- coverage.lrt[i,1]/1000
  coverage.n1[i,3] <- coverage.sc[i,1]/1000
  coverage.n1[i,4] <- coverage.wl[i,1]/1000
}
# n=25
for(i in 1:length(p0)) {
  coverage.n2[i,1] <- p0[i]
  coverage.n2[i,2] <- coverage.lrt[i,2]/1000
  coverage.n2[i,3] <- coverage.sc[i,2]/1000
  coverage.n2[i,4] <- coverage.wl[i,2]/1000
}
# n=50
for(i in 1:length(p0)) {
  coverage.n3[i,1] <- p0[i]
  coverage.n3[i,2] <- coverage.lrt[i,3]/1000
  coverage.n3[i,3] <- coverage.sc[i,3]/1000
  coverage.n3[i,4] <- coverage.wl[i,3]/1000
}
# n=100
for(i in 1:length(p0)) {
  coverage.n4[i,1] <- p0[i]
  coverage.n4[i,2] <- coverage.lrt[i,4]/1000
  coverage.n4[i,3] <- coverage.sc[i,4]/1000
  coverage.n4[i,4] <- coverage.wl[i,4]/1000
}
# n=250
for(i in 1:length(p0)) {
  coverage.n5[i,1] <- p0[i]
  coverage.n5[i,2] <- coverage.lrt[i,5]/1000
  coverage.n5[i,3] <- coverage.sc[i,5]/1000
  coverage.n5[i,4] <- coverage.wl[i,5]/1000
}

cov.n1 <- data.frame(coverage.n1)
cov.n2 <- data.frame(coverage.n2)
cov.n3 <- data.frame(coverage.n3)
cov.n4 <- data.frame(coverage.n4)
cov.n5 <- data.frame(coverage.n5)

library(ggplot2)
# when alpha = 0.05
g1 <- ggplot(cov.n1, aes(x=p0)) + 
        geom_point(aes(y=lrt), color = "blue3") +
        geom_point(aes(y=sc), color = "violetred3") + 
        geom_point(aes(y=wl), color = "gold2") +
        geom_line(aes(y=lrt), linetype = 2, color = "blue3") +
        geom_line(aes(y=sc), linetype = 6, color = "violetred3") +
        geom_line(aes(y=wl), color = "gold2") +
        geom_hline(yintercept = 0.95, linetype = "dashed", color = "red", size = 1) +
        labs(x = "p0", y = "coverage", title = "p0 vs coverage under n=10") + 
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
g2 <- ggplot(cov.n2, aes(x=p0)) + 
        geom_point(aes(y=lrt), color = "blue3") +
        geom_point(aes(y=sc), color = "violetred3") + 
        geom_point(aes(y=wl), color = "gold2") +
        geom_line(aes(y=lrt), linetype = 2, color = "blue3") +
        geom_line(aes(y=sc), linetype = 6, color = "violetred3") +
        geom_line(aes(y=wl), color = "gold2") +
        geom_hline(yintercept = 0.95, linetype = "dashed", color = "red", size = 1) +
        labs(x = "p0", y = "coverage", title = "p0 vs coverage under n=25") + 
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
g3 <- ggplot(cov.n3, aes(x=p0)) + 
        geom_point(aes(y=lrt), color = "blue3") +
        geom_point(aes(y=sc), color = "violetred3") + 
        geom_point(aes(y=wl), color = "gold2") +
        geom_line(aes(y=lrt), linetype = 2, color = "blue3") +
        geom_line(aes(y=sc), linetype = 6, color = "violetred3") +
        geom_line(aes(y=wl), color = "gold2") +
        geom_hline(yintercept = 0.95, linetype = "dashed", color = "red", size = 1) +
        labs(x = "p0", y = "coverage", title = "p0 vs coverage under n=50") + 
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
g4 <- ggplot(cov.n4, aes(x=p0)) + 
        geom_point(aes(y=lrt), color = "blue3") +
        geom_point(aes(y=sc), color = "violetred3") + 
        geom_point(aes(y=wl), color = "gold2") +
        geom_line(aes(y=lrt), linetype = 2, color = "blue3") +
        geom_line(aes(y=sc), linetype = 6, color = "violetred3") +
        geom_line(aes(y=wl), color = "gold2") +
        geom_hline(yintercept = 0.95, linetype = "dashed", color = "red", size = 1) +
        labs(x = "p0", y = "coverage", title = "p0 vs coverage under n=100") + 
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
g5 <- ggplot(cov.n5, aes(x=p0)) + 
        geom_point(aes(y=lrt), color = "blue3") +
        geom_point(aes(y=sc), color = "violetred3") + 
        geom_point(aes(y=wl), color = "gold2") +
        geom_line(aes(y=lrt), linetype = 2, color = "blue3") +
        geom_line(aes(y=sc), linetype = 6, color = "violetred3") +
        geom_line(aes(y=wl), color = "gold2") +
        geom_hline(yintercept = 0.95, linetype = "dashed", color = "red", size = 1) +
        labs(x = "p0", y = "coverage", title = "p0 vs coverage under n=250") + 
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))

library(ggpubr)
g.alpha1 <- ggarrange(g1,g2,g3,g4,g5, nrow = 3, ncol = 2,
          labels = c("fig1", "fig2", "fig3", "fig4", "fig5"))
annotate_figure(g.alpha1, top = text_grob("Coversage Under Sample Size (alpha = 0.05)", face = "bold", size = 15))

# when alpha = 0.01
library(ggplot2)
f1 <- ggplot(cov.n1, aes(x=p0)) + 
  geom_point(aes(y=lrt), color = "blue3") +
  geom_point(aes(y=sc), color = "violetred3") + 
  geom_point(aes(y=wl), color = "gold2") +
  geom_line(aes(y=lrt), linetype = 2, color = "blue3") +
  geom_line(aes(y=sc), linetype = 6, color = "violetred3") +
  geom_line(aes(y=wl), color = "gold2") +
  geom_hline(yintercept = 0.99, linetype = "dashed", color = "red", size = 1) +
  labs(x = "p0", y = "coverage", title = "p0 vs coverage under n=10") + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
f2 <- ggplot(cov.n2, aes(x=p0)) + 
  geom_point(aes(y=lrt), color = "blue3") +
  geom_point(aes(y=sc), color = "violetred3") + 
  geom_point(aes(y=wl), color = "gold2") +
  geom_line(aes(y=lrt), linetype = 2, color = "blue3") +
  geom_line(aes(y=sc), linetype = 6, color = "violetred3") +
  geom_line(aes(y=wl), color = "gold2") +
  geom_hline(yintercept = 0.99, linetype = "dashed", color = "red", size = 1) +
  labs(x = "p0", y = "coverage", title = "p0 vs coverage under n=25") + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
f3 <- ggplot(cov.n3, aes(x=p0)) + 
  geom_point(aes(y=lrt), color = "blue3") +
  geom_point(aes(y=sc), color = "violetred3") + 
  geom_point(aes(y=wl), color = "gold2") +
  geom_line(aes(y=lrt), linetype = 2, color = "blue3") +
  geom_line(aes(y=sc), linetype = 6, color = "violetred3") +
  geom_line(aes(y=wl), color = "gold2") +
  geom_hline(yintercept = 0.99, linetype = "dashed", color = "red", size = 1) +
  labs(x = "p0", y = "coverage", title = "p0 vs coverage under n=50") + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
f4 <- ggplot(cov.n4, aes(x=p0)) + 
  geom_point(aes(y=lrt), color = "blue3") +
  geom_point(aes(y=sc), color = "violetred3") + 
  geom_point(aes(y=wl), color = "gold2") +
  geom_line(aes(y=lrt), linetype = 2, color = "blue3") +
  geom_line(aes(y=sc), linetype = 6, color = "violetred3") +
  geom_line(aes(y=wl), color = "gold2") +
  geom_hline(yintercept = 0.99, linetype = "dashed", color = "red", size = 1) +
  labs(x = "p0", y = "coverage", title = "p0 vs coverage under n=100") + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
f5 <- ggplot(cov.n5, aes(x=p0)) + 
  geom_point(aes(y=lrt), color = "blue3") +
  geom_point(aes(y=sc), color = "violetred3") + 
  geom_point(aes(y=wl), color = "gold2") +
  geom_line(aes(y=lrt), linetype = 2, color = "blue3") +
  geom_line(aes(y=sc), linetype = 6, color = "violetred3") +
  geom_line(aes(y=wl), color = "gold2") +
  geom_hline(yintercept = 0.99, linetype = "dashed", color = "red", size = 1) +
  labs(x = "p0", y = "coverage", title = "p0 vs coverage under n=250") + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

library(ggpubr)
g.alpha2 <- ggarrange(f1,f2,f3,f4,f5, nrow = 3, ncol = 2,
                      labels = c("fig1", "fig2", "fig3", "fig4", "fig5"))
annotate_figure(g.alpha2, top = text_grob("Coversage Under Sample Size (alpha = 0.01)", face = "bold", size = 15))


