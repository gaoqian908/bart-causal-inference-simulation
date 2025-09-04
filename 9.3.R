# R code for 4 methods in 3 response surfaces 
# You can change the iteration manually
niters <- 1000
set.seed(24015941)

# Generate data
n <- 500
p <- 5
rho <- 0.5
Sigma <- outer(1:p, 1:p, function(i, j) rho ^ abs(i - j))
X <- MASS :: mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
colnames(X) <- paste0("X", 1:p)

# Set propensity score
logit_ps <- sin(X[,1] ^ 2) + log(abs(X[,2] + 1)) - 0.5 *X[,3]*X[,4]
ps <- 1 / (1 + exp(-logit_ps))
Z <- rbinom(n, 1, ps)

# Good overlap
good_overlap_index <- which(ps >= 0.1 & ps <= 0.9)
X_good <- X[good_overlap_index,]
Z_good <- Z[good_overlap_index]
ps_good <- ps[good_overlap_index]
n_good <- length(good_overlap_index)

xt <- as.matrix(cbind(X_good,Z_good))
xp1 <- xt[Z_good,]
xp2 <- xp1
xp1[,ncol(xt)] <- 1
xp2[,ncol(xt)] <- 0
xp <- rbind(xp1,xp2)

nt <- sum(xt[,"Z_good"])
N <- nrow(X_good)
Xmat <- as.matrix(X_good)

# Calculate propensity score
df_ps <- data.frame(Z = Z_good, X_good)
ps_model <- glm(Z ~ ., data = df_ps, family = binomial(link = "logit"))
ps_est <- predict(ps_model, type = "response")

# Creat place to save index
results.a <- matrix(0,niters,13)
results.b <- matrix(0,niters,13)
results.c <- matrix(0,niters,13)
dimnames(results.a) <- list(NULL,c("b.te","b.cov","b.cil","r.te","r.cov","r.cil","aipw.te","aipw.cov","aipw.cil","ipw.te","ipw.cov","ipw.cil","tau.est"))
dimnames(results.b) <- list(NULL,c("b.te","b.cov","b.cil","r.te","r.cov","r.cil","aipw.te","aipw.cov","aipw.cil","ipw.te","ipw.cov","ipw.cil","tau.est"))
dimnames(results.c) <- list(NULL,c("b.te","b.cov","b.cil","r.te","r.cov","r.cil","aipw.te","aipw.cov","aipw.cil","ipw.te","ipw.cov","ipw.cil","tau.est"))

dimnames(results.a) <- list(NULL,c("b.te","b.cov","b.cil","r.te","r.cov","r.cil","aipw.te","aipw.cov","aipw.cil","ipw.te","ipw.cov","ipw.cil","tau.est"))
dimnames(results.b) <- list(NULL,c("b.te","b.cov","b.cil","r.te","r.cov","r.cil","aipw.te","aipw.cov","aipw.cil","ipw.te","ipw.cov","ipw.cil","tau.est"))
dimnames(results.c) <- list(NULL,c("b.te","b.cov","b.cil","r.te","r.cov","r.cil","aipw.te","aipw.cov","aipw.cil","ipw.te","ipw.cov","ipw.cil","tau.est"))

precision.a <- matrix(0,niters,4)
dimnames(precision.a) <- list(NULL,c("bart","reg","aipw","ipw"))
precision.b <- matrix(0,niters,4)
dimnames(precision.b) <- list(NULL,c("bart","reg","aipw","ipw"))
precision.c <- matrix(0,niters,4)
dimnames(precision.c) <- list(NULL,c("bart","reg","aipw","ipw"))
```

```{r}
for(i in 1:niters){
  # Print the status for each iteration
  cat("Processing iteration", i, "of", niters, "\n")
  if(i<=50){set.seed(37 + i*5)}
  if(i>50){set.seed(73 + i*5)}
  
  # RESPONSE SURFACES (3 versions)
  
  # YA-linear+homogeneous treatment effect
  sigy <- 3
  tau <- 4
  betaA <- c(0.5, 1.5, -2, 1, 0.5, rep(0, p - 4))
  yahat <- cbind(rep(1, N), Xmat) %*% betaA
  YA0 <- rnorm(N, yahat, sigy)
  YA1 <- rnorm(N, yahat+tau, sigy)   # homogeneous treatment effect
  YA <- YA1
  YA[Z_good==0] <- YA0[Z_good==0]
  tauAis <- YA1[Z_good==1] - YA0[Z_good==1]
  
  # YB-non-linear+homogeneous treatment effect
  sigy <- 3
  tau <- 4
  fB <- function(X) {
    tanh(X[,1])+2*X[,2]^2-3*X[,3]*X[,4]+cos(X[,5])}
  yb0hat <- fB(Xmat)
  yb1hat <- yb0hat + tau   # homogeneous treatment effect
  offsetB <- mean(yb1hat[Z_good == 1] - yb0hat[Z_good == 1]) - tau
  yb1hat <- yb1hat - offsetB   # Center the treatment effect around 4
  
  YB0 <- rnorm(n_good, mean = yb0hat, sd = sigy)
  YB1 <- rnorm(n_good, mean = yb1hat, sd = sigy)
  YB <- YB1
  YB[Z_good == 0] <- YB0[Z_good == 0]
  tauBis <- yb1hat[Z_good==1] - yb0hat[Z_good==1]
  
  
  # YC-non-linear+heterogeneous treatment effect
  sigy <- 3
  tauC_func <- function(X) {
    sin(X[,1]*X[,2]) +0.5*X[,3] - 0.5*X[,4]*X[,5]}
  tauc <- tauC_func(Xmat)
  fC <- fB
  yc0hat <- fC(Xmat)
  yc1hat <- yc0hat + tauc  # heterogeneous treatment effect
  offsetC <- mean(yc1hat[Z_good == 1] - yc0hat[Z_good == 1]) - 4
  yc1hat <- yc1hat - offsetC   # Center the treatment effect around 4
  
  YC0 <- rnorm(n_good, mean = yc0hat, sd = sigy)
  YC1 <- rnorm(n_good, mean = yc1hat, sd = sigy)
  YC <- YC1
  YC[Z_good == 0] <- YC0[Z_good == 0]
  tauCis <- yc1hat[Z_good == 1] - yc0hat[Z_good == 1]
  
  
  # generate sample treatment effects
  tauAs <- mean(YA1[Z_good==1] - YA0[Z_good==1])
  tauBs <- mean(YB1[Z_good==1] - YB0[Z_good==1])
  tauCs <- mean(YC1[Z_good==1] - YC0[Z_good==1])
  
  taus2 <- c(tauAs,tauBs,tauCs)
  
  results.a[i,13] <- taus2[1]
  results.b[i,13] <- taus2[2]
  results.c[i,13] <- taus2[3]
  
  #BART
  library(dbarts)
  bart2a <- dbarts::bart(x.train =xt,   y.train=YA,  x.test=xp,  verbose = FALSE)
  bart2b <- dbarts::bart(x.train=xt,   y.train=YB,  x.test=xp,  verbose = FALSE)
  bart2c <- dbarts::bart(x.train=xt,   y.train=YC,  x.test=xp,  verbose = FALSE)
  
  # a
  tmp <- apply(bart2a$yhat.test[,1:nt]-bart2a$yhat.test[,(nt+1):(2*nt)],1,mean)
  tmpa <- mean(tmp)
  results.a[i,1] <- tauAs-tmpa
  sd <- sqrt(var(tmp))
  ci <- c(tmpa-1.96*sd,tmpa+1.96*sd)
  results.a[i,2] <- (ci[1]<tauAs & ci[2]>tauAs)*1
  results.a[i,3] <- ci[2]-ci[1]
  
  tmp2 <- apply(bart2a$yhat.test[,1:nt]-bart2a$yhat.test[,(nt+1):(2*nt)],2,mean)
  precision.a[i,1] <- sqrt(mean((tmp2-tauAis)^2))
  
  # b
  tmp <- apply(bart2b$yhat.test[,1:nt]-bart2b$yhat.test[,(nt+1):(2*nt)],1,mean)
  tmpb <- mean(tmp)
  results.b[i,1] <- tauBs-tmpb
  sd <- sqrt(var(tmp))
  ci <- c(tmpb-1.96*sd,tmpb+1.96*sd)
  results.b[i,2] <- (ci[1]<tauBs & ci[2]>tauBs)*1
  results.b[i,3] <- ci[2]-ci[1]
  
  tmp2 <- apply(bart2b$yhat.test[,1:nt]-bart2b$yhat.test[,(nt+1):(2*nt)],2,mean)
  precision.b[i,1] <- sqrt(mean((tmp2-tauBis)^2))
  
  # c
  tmp <- apply(bart2c$yhat.test[,1:nt]-bart2c$yhat.test[,(nt+1):(2*nt)],1,mean)
  tmpc <- mean(tmp)
  results.c[i,1] <- tauCs-tmpc
  sd <- sqrt(var(tmp))
  ci <- c(tmpc-1.96*sd,tmpc+1.96*sd)
  results.c[i,2] <- (ci[1]<tauCs & ci[2]>tauCs)*1
  results.c[i,3] <- ci[2]-ci[1]
  
  tmp2 <- apply(bart2c$yhat.test[,1:nt]-bart2c$yhat.test[,(nt+1):(2*nt)],2,mean)
  precision.c[i,1] <- sqrt(mean((tmp2-tauCis)^2))
  
  # Outcome regression
  #a
  data2 <- cbind.data.frame(X_good,Trt=Z_good,YA=YA,YB=YB, YC = YC)
  tmpp <- summary(lm(data2[,c("YA","Trt","X1","X2","X3","X4","X5")]))
  tmp <- tmpp$coef[2,1:2]
  results.a[i,4] <- tauAs-tmp[1]
  ci <- c(tmp[1]-1.96*tmp[2],tmp[1]+1.96*tmp[2])
  results.a[i,5] <- (ci[1]<tauAs & ci[2]>tauAs)*1
  results.a[i,6] <- ci[2]-ci[1]
  precision.a[i,2] <- sqrt(mean((tmp[1]-tauAis)^2))
  
  #b
  tmpp <- summary(lm(data2[,c("YB","Trt","X1","X2","X3","X4","X5")]))
  tmp <- tmpp$coef[2,1:2]
  results.b[i,4] <- tauBs-tmp[1]
  ci <- c(tmp[1]-1.96*tmp[2],tmp[1]+1.96*tmp[2])
  results.b[i,5] <- (ci[1]<tauBs & ci[2]>tauBs)*1
  results.b[i,6] <- ci[2]-ci[1]
  precision.b[i,2]  <-  sqrt(mean((tmp[1]-tauBis)^2))
  
  #c
  tmpp <- summary(lm(data2[,c("YC","Trt","X2","X3","X4","X5")]))
  tmp <- tmpp$coef[2,1:2]
  results.c[i,4] <- tauCs-tmp[1]
  ci <- c(tmp[1]-1.96*tmp[2],tmp[1]+1.96*tmp[2])
  results.c[i,5] <- (ci[1]<tauCs & ci[2]>tauCs)*1
  results.c[i,6] <- ci[2]-ci[1]
  precision.c[i,2] <- sqrt(mean((tmp[1]-tauCis)^2))
  
  # IPW
  # a
  ipw_te <- mean(Z_good * YA / ps_est - (1 - Z_good) * YA / (1 - ps_est))
  results.a[i,10] <- tauAs - ipw_te
  ipw_var <- var(Z_good * YA / ps_est - (1 - Z_good) * YA / (1 - ps_est)) / n_good
  ci <- c(ipw_te - 1.96 * sqrt(ipw_var), ipw_te + 1.96 * sqrt(ipw_var))
  results.a[i,11] <- (ci[1] < tauAs & ci[2] > tauAs) * 1
  results.a[i,12] <- ci[2] - ci[1]
  precision.a[i,4] <- sqrt(mean((ipw_te - tauAis)^2))
  
  # b
  ipw_te <- mean(Z_good * YB / ps_est - (1 - Z_good) * YB / (1 - ps_est))
  results.b[i,10] <- tauBs - ipw_te
  ipw_var <- var(Z_good * YB / ps_est - (1 - Z_good) * YB / (1 - ps_est)) / n_good
  ci <- c(ipw_te - 1.96 * sqrt(ipw_var), ipw_te + 1.96 * sqrt(ipw_var))
  results.b[i,11] <- (ci[1] < tauBs & ci[2] > tauBs) * 1
  results.b[i,12] <- ci[2] - ci[1]
  precision.b[i,4] <- sqrt(mean((ipw_te - tauBis)^2))
  
  # c
  ipw_te <- mean(Z_good * YC / ps_est - (1 - Z_good) * YC / (1 - ps_est))
  results.c[i,10] <- tauCs - ipw_te
  ipw_var <- var(Z_good * YC / ps_est - (1 - Z_good) * YC / (1 - ps_est)) / n_good
  ci <- c(ipw_te - 1.96 * sqrt(ipw_var), ipw_te + 1.96 * sqrt(ipw_var))
  results.c[i,11] <- (ci[1] < tauCs & ci[2] > tauCs) * 1
  results.c[i,12] <- ci[2] - ci[1]
  precision.c[i,4] <- sqrt(mean((ipw_te - tauCis)^2))
  
  # AIPW
  # Create test sets: everyone as treated and everyone as control 
  xt1 <- as.data.frame(X_good)
  xt0 <- as.data.frame(X_good)
  
  # a
  outcome_treated_modelA <- lm(YA[Z_good == 1] ~ ., data = as.data.frame(X_good[Z_good == 1,]))
  outcome_control_modelA <- lm(YA[Z_good == 0] ~ ., data = as.data.frame(X_good[Z_good == 0,]))
  mu1_hat <- predict(outcome_treated_modelA, newdata = xt1)
  mu0_hat <- predict(outcome_control_modelA, newdata = xt0)
  
  # Compute AIPW estimate
  aipw_te <- mean(
    mu1_hat - mu0_hat +
      Z_good * (YA - mu1_hat) / ps_est -
      (1 - Z_good) * (YA - mu0_hat) / (1 - ps_est))
  
  results.a[i, 7] <- tauAs - aipw_te
  aipw_var <- var(mu1_hat - mu0_hat + Z_good * (YA - mu1_hat) / 
                    ps_est - (1 - Z_good) * (YA - mu0_hat) / (1 - ps_est)) / n_good
  ci <- c(aipw_te - 1.96 * sqrt(aipw_var), aipw_te + 1.96 * sqrt(aipw_var))
  results.a[i, 8] <- (ci[1] < tauAs & ci[2] > tauAs) * 1
  results.a[i, 9] <- ci[2] - ci[1]
  precision.a[i, 3] <- sqrt(mean((aipw_te - tauAis)^2))
  
  # b
  outcome_treated_modelB <- lm(YB[Z_good == 1] ~ ., data = as.data.frame(X_good[Z_good == 1,]))
  outcome_control_modelB <- lm(YB[Z_good == 0] ~ ., data = as.data.frame(X_good[Z_good == 0,]))
  mu1_hat <- predict(outcome_treated_modelB, newdata = xt1)
  mu0_hat <- predict(outcome_control_modelB, newdata = xt0)
  
  # AIPW estimate
  aipw_te <- mean(
    mu1_hat - mu0_hat +Z_good * (YB - mu1_hat) / ps_est -
      (1 - Z_good) * (YB - mu0_hat) / (1 - ps_est))
  
  results.b[i, 7] <- tauBs - aipw_te
  aipw_var <- var(mu1_hat - mu0_hat + Z_good * (YB - mu1_hat) / ps_est - 
                    (1 - Z_good) * (YB - mu0_hat) / (1 - ps_est)) / n_good
  ci <- c(aipw_te - 1.96 * sqrt(aipw_var), aipw_te + 1.96 * sqrt(aipw_var))
  results.b[i, 8] <- (ci[1] < tauBs & ci[2] > tauBs) * 1
  results.b[i, 9] <- ci[2] - ci[1]
  precision.b[i, 3] <- sqrt(mean((aipw_te - tauBis)^2))
  
  # c
  outcome_treated_modelC <- lm(YC[Z_good == 1] ~ ., data = as.data.frame(X_good[Z_good == 1,]))
  outcome_control_modelC <- lm(YC[Z_good == 0] ~ ., data = as.data.frame(X_good[Z_good == 0,]))
  mu1_hat <- predict(outcome_treated_modelC, newdata = xt1)
  mu0_hat <- predict(outcome_control_modelC, newdata = xt0)
  
  # AIPW estimate
  aipw_te <- mean(
    mu1_hat - mu0_hat +Z_good * (YC - mu1_hat) / ps_est -
      (1 - Z_good) * (YC - mu0_hat) / (1 - ps_est))
  
  results.c[i, 7] <- tauCs - aipw_te
  aipw_var <- var(mu1_hat - mu0_hat + Z_good * (YC - mu1_hat) / ps_est - 
                    (1 - Z_good) * (YC - mu0_hat) / (1 - ps_est)) / n_good
  ci <- c(aipw_te - 1.96 * sqrt(aipw_var), aipw_te + 1.96 * sqrt(aipw_var))
  results.c[i, 8] <- (ci[1] < tauCs & ci[2] > tauCs) * 1
  results.c[i, 9] <- ci[2] - ci[1]
  precision.c[i, 3] <- sqrt(mean((aipw_te - tauCis)^2))
}

library(ggplot2)
library(reshape2)
library(patchwork)

# Function to create combined plot for a given results and precision matrix
plot_simulation_results <- function(results, precision, label, n, p, niters) {
  # Create data frames
  bias_df <- data.frame(
    Iteration = 1:niters,
    BART = results[,1],
    Reg = results[,4],
    IPW = results[,10],
    AIPW = results[,7]
  )
  
  coverage_df <- data.frame(
    Iteration = 1:niters,
    BART = results[,2],
    Reg = results[,5],
    IPW = results[,11],
    AIPW = results[,8]
  )
  
  ci_df <- data.frame(
    Iteration = 1:niters,
    BART = results[,3],
    Reg = results[,6],
    IPW = results[,12],
    AIPW = results[,9]
  )
  
  precision_df <- data.frame(
    Iteration = 1:niters,
    BART = precision[,1],
    Reg = precision[,2],
    IPW = precision[,4],
    AIPW = precision[,3]
  )
  
  # Melt data frames for ggplot
  bias_long <- melt(bias_df, id.vars = "Iteration", variable.name = "Method", value.name = "Bias")
  coverage_long <- melt(coverage_df, id.vars = "Iteration", variable.name = "Method", value.name = "Coverage")
  ci_long <- melt(ci_df, id.vars = "Iteration", variable.name = "Method", value.name = "CI_Length")
  precision_long <- melt(precision_df, id.vars = "Iteration", variable.name = "Method", value.name = "RMSE")
  
  # Create plots
  p1 <- ggplot(bias_long, aes(x = Method, y = Bias, fill = Method)) +
    geom_boxplot() +
    theme_minimal() +
    ggtitle(paste0("Bias of Estimators (", label, ")")) +
    ylab("Bias") +
    xlab("") +
    theme(legend.position = "none")
  
  p2 <- ggplot(coverage_long, aes(x = Method, y = Coverage, fill = Method)) +
    geom_bar(stat = "summary", fun = mean, position = "dodge") +
    theme_minimal() +
    ggtitle(paste0("Coverage Probability of 95% CI (", label, ")")) +
    ylab("Coverage") +
    xlab("") +
    scale_y_continuous(limits = c(0,1)) +
    theme(legend.position = "none")
  
  p3 <- ggplot(ci_long, aes(x = Method, y = CI_Length, fill = Method)) +
    geom_boxplot() +
    theme_minimal() +
    ggtitle(paste0("95% CI Length (", label, ")")) +
    ylab("Length of Confidence Interval") +
    xlab("") +
    theme(legend.position = "none")
  
  p4 <- ggplot(precision_long, aes(x = Method, y = RMSE, fill = Method)) +
    geom_boxplot() +
    theme_minimal() +
    ggtitle(paste0("RMSE of Estimators (", label, ")")) +
    ylab("Root Mean Squared Error") +
    xlab("") +
    theme(legend.position = "none")
  
  # Combine plots and add overall title with n, p, niters
  combined <- (p1 | p2) / (p3 | p4) +
    plot_annotation(
      title = paste("Simulation Results:", label),
      subtitle = paste("Sample size (n):", n, " | Number of predictors (p):", p, " | Iterations:", niters),
      theme = theme(plot.title = element_text(size = 16, face = "bold"),
                    plot.subtitle = element_text(size = 12))
    )
  
  return(combined)
}
combined_plot_a <- plot_simulation_results(results.a, precision.a, "YA", n, p, niters)
combined_plot_b <- plot_simulation_results(results.b, precision.b, "YB", n, p, niters)
combined_plot_c <- plot_simulation_results(results.c, precision.c, "YC", n, p, niters)

print(combined_plot_a)
print(combined_plot_b)
print(combined_plot_c)

# Generate tables  
get_summary_table <- function(results, precision, label) {
  data.frame(
    Surface = label,
    Method = c("BART", "Reg", "IPW", "AIPW"),
    Bias = c(mean(results[, "b.te"]),
             mean(results[, "r.te"]),
             mean(results[, "ipw.te"]),
             mean(results[, "aipw.te"])),
    Coverage = c(mean(results[, "b.cov"]),
                 mean(results[, "r.cov"]),
                 mean(results[, "ipw.cov"]),
                 mean(results[, "aipw.cov"])),
    CI_Length = c(mean(results[, "b.cil"]),
                  mean(results[, "r.cil"]),
                  mean(results[, "ipw.cil"]),
                  mean(results[, "aipw.cil"])),
    RMSE = c(mean(precision[, "bart"]),
             mean(precision[, "reg"]),
             mean(precision[, "ipw"]),
             mean(precision[, "aipw"]))
  )
}
summary_YA <- get_summary_table(results.a, precision.a, "YA")
summary_YB <- get_summary_table(results.b, precision.b, "YB")
summary_YC <- get_summary_table(results.c, precision.c, "YC")


print(summary_YA)
print(summary_YB)
print(summary_YC)







