############################################################
# Block 2 – Part (a): Portfolio return characteristics
# Assets: spl, kru (WIG20)
############################################################

# Packages
suppressPackageStartupMessages({
  require(zoo)
  require(xts)
  require(knitr)
  require(ggplot2)
  library(rmgarch)
  require(tseries)   # jarque.bera.test
  require(moments)   # skewness, kurtosis
  require(forecast)  # Acf (nicer ACF plot)
  require(rugarch)
  require(copula)
})


Sys.setlocale("LC_ALL", "en_US.UTF-8")

# ----------------------------------------------------------
# 1) Load data
# ----------------------------------------------------------
filename <- "wig20.Rdata"
load(filename)   # expects object named 'data' (zoo) inside

# check
stopifnot(exists("data"))
stopifnot(inherits(data, "zoo") || inherits(data, "xts"))

# ----------------------------------------------------------
# 2) Select sample and two assets: spl, kru
# ----------------------------------------------------------
startDate <- as.Date("2005-01-01")
endDate   <- as.Date("2050-01-01")

y <- window(data, start = startDate, end = endDate)

# select stocks and remove NAs
y <- na.omit(y[, c("spl","kru")])

# ----------------------------------------------------------
# 3) Compute returns and portfolio
# dy: percent log returns
# ----------------------------------------------------------
dy <- 100 * diff(log(y))
dy <- na.omit(dy)

# portfolio weights
w <- c(0.5, 0.5)

# rP as 1-column zoo (IMPORTANT: keep 2D so colnames works)
rP <- zoo(
  matrix(w[1]*coredata(dy[,1]) + w[2]*coredata(dy[,2]), ncol = 1),
  order.by = index(dy)
)
colnames(rP) <- "rP"

R  <- as.numeric(coredata(rP))  # numeric vector for stats/ACF/QQ

# ----------------------------------------------------------
# (a) Outputs required by the assignment
# time series plot, moments (annualized mean/sd), QQ, ACF, ACF squares
# ----------------------------------------------------------

# --- (a1) Time series plot
par(mfrow=c(1,1), cex=0.9, bty="l")
plot(rP, main="Portfolio log returns rP (%, daily): 0.5*spl + 0.5*kru",
     xlab="", ylab="Return (%)")
abline(h=0, lty=2)

# --- (a2) Moments with annualized mean & sd
# Use 252 trading days convention
# --- (a2) Moments with annualized mean & sd (teacher-style vertical table)

dates_r <- index(rP)   # use portfolio return dates (rP is zoo)
Nyear   <- 365 / as.numeric(mean(diff(dates_r)))  # teacher-style annualization

R <- as.numeric(coredata(rP))  # portfolio daily returns (%)

JB <- tseries::jarque.bera.test(R)

mom <- data.frame(
  value = c(
    Nyear,
    mean(R, na.rm=TRUE) * Nyear,                 # annualized mean
    sd(R,   na.rm=TRUE) * sqrt(Nyear),           # annualized sd
    min(R,  na.rm=TRUE),                         # daily min
    max(R,  na.rm=TRUE),                         # daily max
    moments::skewness(R, na.rm=TRUE),
    moments::kurtosis(R, na.rm=TRUE),
    as.numeric(JB$statistic),
    as.numeric(JB$p.value)
  )
)

rownames(mom) <- c("Nyear", "mu (annual)", "sigma (annual)",
                   "min (daily)", "max (daily)",
                   "skewness", "kurtosis", "JB stat", "JB p-value")

knitr::kable(mom, digits = 3, caption="(a) Portfolio moments (daily %, annualized mean/sd)")

# --- (a3) QQ plot (standardized returns vs Normal)
Rstar <- (R - mean(R, na.rm=TRUE)) / sd(R, na.rm=TRUE)
qqnorm(Rstar, main="QQ plot: standardized portfolio returns vs Normal")
qqline(Rstar, col="red", lwd=2)

# --- (a4) ACF and ACF of squared returns
par(mfrow=c(1,1), cex=0.9, bty="l")
forecast::Acf(R,   lag.max=60, main="ACF: portfolio returns")
forecast::Acf(R^2, lag.max=60, main="ACF: squared portfolio returns")
par(mfrow=c(1,1))

# ----------------------------------------------------------
# Optional: quick textual interpretation (for your speech)
# ----------------------------------------------------------
cat("\n--- Quick interpretation hints ---\n")
cat("* If QQ tails deviate from the line and JB p-value is small -> non-normal / fat tails.\n")
cat("* If ACF of returns is weak but ACF of squared returns is significant -> volatility clustering -> motivates GARCH.\n")
cat("----------------------------------\n")


############################################################
# (b) Best univariate GARCH model for portfolio returns
# Requires: rP (zoo 1-col) or R (numeric) from part (a)
############################################################

# -----------------------------
# 0) Prepare data
# -----------------------------
# If you have rP (zoo), use it; otherwise use R.
if (exists("rP")) {
  x <- as.numeric(coredata(rP))
} else if (exists("R")) {
  x <- as.numeric(R)
} else {
  stop("Cannot find rP or R. Run part (a) first to create portfolio returns.")
}

x <- na.omit(x)

# Optional: keep a zoo version for plotting with dates
if (exists("rP")) {
  x_zoo <- zoo(x, order.by = index(rP)[!is.na(coredata(rP))])
} else {
  x_zoo <- zoo(x, order.by = 1:length(x))
}

# -----------------------------
# 1) Candidate model set
# -----------------------------
# Variance models (you can add/remove)
var_models <- c("sGARCH", "iGARCH", "eGARCH", "gjrGARCH")  # common set
# Distributions
dist_models <- c("norm", "std", "ged")           # Normal, Student-t, GED

# Mean model: constant mean is usually enough for returns
mean_spec <- list(armaOrder = c(0,0), include.mean = TRUE)

# Helper: fit one model safely
fit_one <- function(x, vmodel, dist) {
  spec <- ugarchspec(
    mean.model      = mean_spec,
    variance.model  = list(model = vmodel, garchOrder = c(1,1)),
    distribution.model = dist
  )
  fit <- tryCatch(
    ugarchfit(spec = spec, data = x, solver = "hybrid"),
    error = function(e) NULL
  )
  fit
}

# -----------------------------
# 2) Fit all candidates and collect selection criteria
# -----------------------------
fits <- list()
sel  <- data.frame(
  Model = character(),
  Dist  = character(),
  LL    = numeric(),
  AIC   = numeric(),
  BIC   = numeric(),
  stringsAsFactors = FALSE
)

k <- 0
for (vm in var_models) {
  for (dm in dist_models) {
    k <- k + 1
    cat(sprintf("Fitting %s(1,1) with %s ...\n", vm, dm))
    f <- fit_one(x, vm, dm)
    key <- paste0(vm, "_", dm)
    fits[[key]] <- f
    
    if (!is.null(f)) {
      ic <- infocriteria(f)  # returns AIC, BIC, etc.
      sel <- rbind(sel, data.frame(
        Model = vm,
        Dist  = dm,
        LL    = as.numeric(likelihood(f)),
        AIC   = as.numeric(ic["Akaike", 1]),
        BIC   = as.numeric(ic["Bayes", 1]),
        stringsAsFactors = FALSE
      ))
    }
  }
}

# If nothing fit, stop
if (nrow(sel) == 0) stop("All GARCH fits failed. Check your data or packages.")

# -----------------------------
# 3) Model selection: choose optimal (by BIC, also show AIC)
# -----------------------------
sel_BIC <- sel[order(sel$BIC), ]
sel_AIC <- sel[order(sel$AIC), ]

cat("\n=== Model selection table (sorted by BIC) ===\n")
print(knitr::kable(sel_BIC, digits = 4))

best_row <- sel_BIC[1, ]
best_key <- paste0(best_row$Model, "_", best_row$Dist)
best_fit <- fits[[best_key]]

cat("\nBest model by BIC:\n")
cat(sprintf("  %s(1,1) with %s distribution\n", best_row$Model, best_row$Dist))

# -----------------------------
# 4) Full parameter estimates table (coef, se, t, p)
# -----------------------------
# rugarch gives a nice matrix from fit@fit$matcoef
param_mat <- best_fit@fit$matcoef
param_tbl <- data.frame(
  Parameter = rownames(param_mat),
  Estimate  = param_mat[,1],
  StdError  = param_mat[,2],
  tValue    = param_mat[,3],
  pValue    = param_mat[,4],
  row.names = NULL
)

cat("\n=== Parameter estimates (full table) ===\n")
print(knitr::kable(param_tbl, digits = 6))

# -----------------------------
# 5) Conditional standard deviation plot (sigma_t)
# -----------------------------
sigma_t <- sigma(best_fit)

sigma_z <- zoo(
  matrix(as.numeric(sigma_t), ncol = 1),
  order.by = index(x_zoo)
)
colnames(sigma_z) <- "sigma_t"

par(mfrow=c(1,1), cex=0.9, bty="l")
plot(sigma_z, type="l",
     main = sprintf("Conditional SD (sigma_t): %s(1,1) - %s", best_row$Model, best_row$Dist),
     xlab="", ylab="sigma_t")

# -----------------------------
# 6) (Optional) Quick diagnostic: standardized residuals ACF
# -----------------------------
# Useful to show model is reasonable (not required, but can help)
z <- residuals(best_fit, standardize = TRUE)
par(mfrow=c(1,2), cex=0.85, bty="l")
acf(as.numeric(z),    lag.max=60, main="ACF: std. residuals")
acf(as.numeric(z)^2,  lag.max=60, main="ACF: (std. residuals)^2")
par(mfrow=c(1,1))

############################################################
# (c) Estimate the best copula for two assets (spl, kru)
# Output: (1) Copula comparison table (LL values)
#         (2) Scatter: realizations (u,v) vs simulation from best copula
############################################################

# ----------------------------------------------------------
# 0) Data: need dy (zoo/xts) with columns "spl" and "kru"
#     dy should be daily log returns in %, as you computed:
#     dy <- 100*diff(log(y))
# ----------------------------------------------------------
stopifnot(exists("dy"))
stopifnot(all(c("spl","kru") %in% colnames(dy)))

R1 <- as.numeric(coredata(dy[,"spl"]))
R2 <- as.numeric(coredata(dy[,"kru"]))

# align & remove NA (important!)
ok <- complete.cases(R1, R2)
R1 <- R1[ok]; R2 <- R2[ok]
dates_uv <- index(dy)[ok]

# ----------------------------------------------------------
# 1) Fit univariate GARCH models to each asset
#    (consistent with your b: heavy tails -> Student-t)
#    Here: eGARCH(1,1) with Student-t errors ("std")
# ----------------------------------------------------------
spec_uni <- ugarchspec(
  mean.model     = list(armaOrder=c(0,0), include.mean=TRUE),
  variance.model = list(model="eGARCH", garchOrder=c(1,1)),
  distribution.model = "std"
)

fit1 <- ugarchfit(spec = spec_uni, data = R1, solver = "hybrid")
fit2 <- ugarchfit(spec = spec_uni, data = R2, solver = "hybrid")

# Standardized residuals
z1 <- as.numeric(residuals(fit1, standardize = TRUE))
z2 <- as.numeric(residuals(fit2, standardize = TRUE))

# Student-t degrees of freedom (shape)
nu1 <- as.numeric(coef(fit1)["shape"])
nu2 <- as.numeric(coef(fit2)["shape"])

# ----------------------------------------------------------
# 2) PIT transform to uniforms u,v in (0,1)
# PIT to uniforms using Student-t CDF
u1 <- rugarch::pdist("std", z1, shape = nu1)
u2 <- rugarch::pdist("std", z2, shape = nu2)

# avoid exact 0/1 (ML can fail if boundaries appear)
eps <- 1e-6
u1 <- pmin(pmax(u1, eps), 1-eps)
u2 <- pmin(pmax(u2, eps), 1-eps)

U <- cbind(u1, u2)
colnames(U) <- c("U_spl", "U_kru")

# ----------------------------------------------------------
# 3) Fit candidate copulas and compare LL
#    Elliptical: Gaussian, t
#    Archimedean: Clayton, Gumbel, Frank
# ----------------------------------------------------------
fit_cop_safe <- function(cop, U) {
  out <- tryCatch(
    fitCopula(cop, data = U, method = "ml"),
    error = function(e) NULL
  )
  out
}

cand <- list(
  Gaussian = normalCopula(param=0.2, dim=2, dispstr="un"),
  tCopula  = tCopula(param=0.2, dim=2, dispstr="un", df=6, df.fixed=FALSE),
  Clayton  = claytonCopula(param=1, dim=2),
  Gumbel   = gumbelCopula(param=1.2, dim=2),
  Frank    = frankCopula(param=1, dim=2)
)

fits <- lapply(cand, fit_cop_safe, U = U)

# LL table
LLtab <- data.frame(
  Copula = names(fits),
  LL     = sapply(fits, function(f) if (is.null(f)) NA_real_ else as.numeric(logLik(f))),
  stringsAsFactors = FALSE
)
LLtab <- LLtab[order(-LLtab$LL), ]  # larger LL is better

cat("\n=== Copula comparison (Log-Likelihood, higher is better) ===\n")
print(knitr::kable(LLtab, digits = 3))

best_name <- LLtab$Copula[1]
best_fit  <- fits[[best_name]]
best_cop  <- best_fit@copula

cat("\nBest copula by LL: ", best_name, "\n", sep="")

# ----------------------------------------------------------
# 4) Simulation from best copula vs realizations (scatter plots)
# ----------------------------------------------------------
set.seed(123)
n <- nrow(U)
U_sim <- rCopula(n, best_cop)

par(mfrow=c(1,2), cex=0.9, bty="l", mar=c(4,4,3,1))

plot(U[,1], U[,2],
     pch=16, cex=0.4,
     main="Realizations (PIT uniforms)",
     xlab="U_spl", ylab="U_kru")

plot(U_sim[,1], U_sim[,2],
     pch=16, cex=0.4,
     main=paste0("Simulation from best copula: ", best_name),
     xlab="U_spl (sim)", ylab="U_kru (sim)")

par(mfrow=c(1,1))

# ----------------------------------------------------------
# (Optional) Quick sanity check: Kendall's tau (real vs simulated)
# ----------------------------------------------------------
tau_real <- cor(U[,1], U[,2], method="kendall")
tau_sim  <- cor(U_sim[,1], U_sim[,2], method="kendall")
cat(sprintf("\nKendall tau (real) = %.3f, (sim) = %.3f\n", tau_real, tau_sim))

############################################################
# (d) VaR / ES – Multi-step Simulation (Final Fix)
############################################################
# === DEFINITION OF PARAMETERS ===
H_set <- c(1, 10)
probs <- c(0.01, 0.05)
Nsim  <- 50000  # <--- Setting Nsim here

simulate_paths <- function(method, H, Nsim){
  
  # --- 1. Historical Simulation ---
  if(method=="Historical"){
    boot <- matrix(sample(R, Nsim*H, replace=TRUE), nrow=Nsim, ncol=H)
    return(boot)
  }
  
  # --- 2. Normal Simulation ---
  if(method=="Normal"){
    return(matrix(rnorm(Nsim*H, mean(R), sd(R)), nrow=Nsim, ncol=H))
  }
  
  # --- 3. t-Student Simulation ---
  if(method=="t"){
    require(MASS)
    fit <- tryCatch(
      suppressWarnings(fitdistr(R, "t")), 
      error=function(e) list(estimate=c(m=mean(R),s=sd(R),df=10))
    )
    m_t <- fit$estimate["m"]
    s_t <- fit$estimate["s"]
    v_t <- fit$estimate["df"]
    return(matrix(m_t + s_t * rt(Nsim*H, df=v_t), nrow=Nsim, ncol=H))
  }
  
  # --- 4. EWMA Simulation ---
  if(method=="EWMA"){
    spec_ewma <- ugarchspec(variance.model=list(model="iGARCH", garchOrder=c(1,1)), 
                            mean.model=list(armaOrder=c(0,0), include.mean=TRUE), 
                            fixed.pars=list(omega=0, alpha1=0.06, beta1=0.94))
    fit_ewma <- ugarchfit(spec_ewma, data=R, solver="hybrid")
    
    # 强制转为 numeric
    mu_last  <- as.numeric(fitted(fit_ewma)[length(R)])
    sig_last <- as.numeric(sigma(fit_ewma)[length(R)])
    
    return(matrix(rnorm(Nsim*H, mean=mu_last, sd=sig_last), nrow=Nsim, ncol=H))
  }
  
  # --- 5. GARCH (Univariate) Simulation ---
  if(method=="GARCH"){
    sim <- ugarchsim(best_garch, n.sim=H, m.sim=Nsim, startMethod="sample", rseed=1:Nsim)
    return(t(sim@simulation$seriesSim)) 
  }
  
  # --- 6. DCC-GARCH Simulation ---
  if(method=="DCC"){
    spec_uni <- ugarchspec(variance.model=list(model=best_row$Model, garchOrder=c(1,1)), 
                           mean.model=list(armaOrder=c(0,0), include.mean=TRUE), 
                           distribution.model=best_row$Dist)
    spec_dcc <- dccspec(multispec(replicate(2, spec_uni)), dccOrder=c(1,1), distribution="mvt")
    
    # 使用 solnp 求解器
    fit_dcc  <- dccfit(spec_dcc, data=dy, solver="solnp") 
    
    sim_dcc <- dccsim(fit_dcc, n.sim=H, m.sim=Nsim, startMethod="sample", rseed=1:Nsim)
    draws <- sim_dcc@msim$simX 
    
    Rp_mat <- matrix(NA, nrow=Nsim, ncol=H)
    for(m_idx in 1:Nsim){
      Rp_mat[m_idx, ] <- draws[[m_idx]] %*% w
    }
    return(Rp_mat)
  }
  
  # --- 7. Copula Simulation (关键修复点) ---
  if(method=="Copula"){
    UVsim <- rCopula(Nsim*H, best_cop)
    
    # [FIX]: 必须使用 as.numeric() 剥离时间索引，否则会与 rCopula 向量冲突
    mu1 <- as.numeric(fitted(fit1)[length(R)])
    sig1 <- as.numeric(sigma(fit1)[length(R)])
    mu2 <- as.numeric(fitted(fit2)[length(R)])
    sig2 <- as.numeric(sigma(fit2)[length(R)])
    
    # 现在这里是纯数字运算，不会报错了
    r1_sim <- mu1 + sig1 * qt(UVsim[,1], df=nu1)
    r2_sim <- mu2 + sig2 * qt(UVsim[,2], df=nu2)
    
    Rp_vec <- w[1]*r1_sim + w[2]*r2_sim
    return(matrix(Rp_vec, nrow=Nsim, ncol=H, byrow=TRUE))
  }
  
  stop(paste("Unknown method:", method))
}

# ----------------------------------------------------------
# 执行计算循环
# ----------------------------------------------------------
VaR_ES_Result <- data.frame()

cat("Starting Simulation (Nsim =", Nsim, ")... This may take a while.\n")

model_list <- c("Historical", "Normal", "t", "EWMA", "GARCH", "DCC", "Copula")

for(m in model_list){
  cat("Processing Method:", m, "...\n")
  for(H in H_set){
    
    # 1. 获取模拟路径 (Nsim x H)
    paths <- simulate_paths(m, H, Nsim)
    
    # 2. 计算 H 天的累积收益
    # Topic 5 Line 153: RdrawsC <- apply(Rdraws, 2, cumsum) -> 我们取最后一天累计值
    # 简单的加总 (log returns add up)
    Rp_H <- rowSums(paths)
    
    # 3. 计算 VaR 和 ES
    for(p in probs){
      # VaR: quantile
      VaR_val <- quantile(Rp_H, p, na.rm=TRUE)
      
      # ES: conditional mean
      # Topic 5 Line 157: ESHt[h] = mean(temp[1:M0])
      ES_val  <- mean(Rp_H[Rp_H <= VaR_val], na.rm=TRUE)
      
      VaR_ES_Result <- rbind(VaR_ES_Result, data.frame(
        Method = m,
        Horizon = paste0(H, "-day"),
        Prob = paste0(p*100, "%"),
        VaR = as.numeric(VaR_val),
        ES  = as.numeric(ES_val)
      ))
    }
  }
}

# ----------------------------------------------------------
# 输出结果表格
# ----------------------------------------------------------
cat("\n=== (d) VaR / ES Table (Based on Topic 5,6,7 Logic) ===\n")
print(knitr::kable(VaR_ES_Result, digits=4, align='c'))

library(ggplot2)

# 确保 VaR_ES_Result 里的数据是数值型
plot_data <- VaR_ES_Result
plot_data$VaR <- as.numeric(plot_data$VaR)
plot_data$ES  <- as.numeric(plot_data$ES)

# 为了画图好看，把 Horizon 的因子顺序固定一下
plot_data$Horizon <- factor(plot_data$Horizon, levels = c("1-day", "10-day"))

# 画图：VaR 对比图
p1 <- ggplot(plot_data, aes(x = Method, y = VaR, fill = Prob)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~Horizon, scales = "free_y") +  # 分面展示 1天 和 10天
  theme_minimal() +
  labs(title = "VaR Comparison across Models",
       subtitle = "DCC and Copula show significantly higher tail risk at 10-day horizon",
       y = "VaR (%)", x = "") +
  scale_fill_manual(values = c("1%" = "darkred", "5%" = "steelblue")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 画图：ES 对比图 (通常 ES 更能说明尾部风险)
p2 <- ggplot(plot_data, aes(x = Method, y = ES, fill = Prob)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~Horizon, scales = "free_y") +
  theme_minimal() +
  labs(title = "Expected Shortfall (ES) Comparison",
       subtitle = "Note the extreme risk predicted by DCC-GARCH for 10-day",
       y = "Expected Shortfall (%)", x = "") +
  scale_fill_manual(values = c("1%" = "darkred", "5%" = "steelblue")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 显示图表
print(p1)
print(p2)

############################################################
# Block 4 – Part (e): Final Perfect Backtesting
# Models: Historical, EWMA, DCC
# Tests: Kupiec, Christoffersen 1 (Ind), Christoffersen 2 (CC), Pelletier, Frey-McNeil
############################################################
library(rugarch)

# --- C. DCC-GARCH (Refit for Backtesting) ---
spec_uni <- ugarchspec(
  variance.model = list(model = best_row$Model, garchOrder = c(1,1)),
  mean.model     = list(armaOrder = c(0,0), include.mean = TRUE),
  distribution.model = best_row$Dist
)

spec_dcc <- dccspec(
  multispec(replicate(2, spec_uni)),
  dccOrder = c(1,1),
  distribution = "mvt"
)

fit_dcc <- dccfit(spec_dcc, data = dy, solver = "solnp")

# ----------------------------------------------------------
# 1) 准备数据 (Risk Time Series)
# ----------------------------------------------------------
W_roll <- 250
R_vec  <- as.numeric(R)
dates  <- index(rP)
alphas <- c(0.01, 0.05)

# ES 辅助函数
es_factor_norm <- function(p) dnorm(qnorm(p))/p
es_factor_std  <- function(p, nu) {
  q <- qt(p, df=nu)
  ((nu + q^2)/(nu - 1)) * (dt(q, df=nu)/p)
}

Risk_TS <- data.frame(Date = dates)

cat("Generating Risk Series (VaR & ES)...\n")

# --- A. Historical Simulation (With Lag) ---
for(a in alphas){
  hs_var_raw <- zoo::rollapply(R_vec, W_roll, function(x) quantile(x, a, na.rm=T), align="right", fill=NA)
  hs_es_raw  <- zoo::rollapply(R_vec, W_roll, function(x) mean(x[x<=quantile(x,a,na.rm=T)], na.rm=T), align="right", fill=NA)
  Risk_TS[[paste0("HS_VaR_", a)]] <- c(NA, head(hs_var_raw, -1))
  Risk_TS[[paste0("HS_ES_", a)]]  <- c(NA, head(hs_es_raw, -1))
}

# --- B. EWMA (In-sample Sigma) ---
spec_ewma <- ugarchspec(variance.model=list(model="iGARCH", garchOrder=c(1,1)), 
                        mean.model=list(armaOrder=c(0,0), include.mean=T), 
                        fixed.pars=list(omega=0, alpha1=0.06, beta1=0.94))
fit_ewma_bt <- ugarchfit(spec_ewma, data=R_vec)
sig_ewma_ts <- sigma(fit_ewma_bt); mu_ewma_ts <- fitted(fit_ewma_bt)

for(a in alphas){
  Risk_TS[[paste0("EWMA_VaR_", a)]] <- mu_ewma_ts + sig_ewma_ts * qnorm(a)
  Risk_TS[[paste0("EWMA_ES_", a)]]  <- mu_ewma_ts - sig_ewma_ts * es_factor_norm(a)
}

# --- C. DCC-GARCH ---
if(!exists("fit_dcc")) stop("请先运行 Part (d) 拟合 DCC 模型")
H_list <- rcov(fit_dcc); T_len <- dim(H_list)[3]
sig_dcc_ts <- numeric(T_len)
for(i in 1:T_len) sig_dcc_ts[i] <- sqrt(as.numeric(t(w) %*% H_list[,,i] %*% w))
df_dcc <- coef(fit_dcc)["[Joint]mshape"]; if(is.na(df_dcc)) df_dcc <- Inf

for(a in alphas){
  if(is.infinite(df_dcc)){
    Risk_TS[[paste0("DCC_VaR_", a)]] <- sig_dcc_ts * qnorm(a)
    Risk_TS[[paste0("DCC_ES_", a)]]  <- -sig_dcc_ts * es_factor_norm(a)
  } else {
    Risk_TS[[paste0("DCC_VaR_", a)]] <- sig_dcc_ts * qt(a, df=df_dcc)
    Risk_TS[[paste0("DCC_ES_", a)]]  <- -sig_dcc_ts * es_factor_std(a, df_dcc)
  }
}

# ----------------------------------------------------------
# 2) 执行所有测试 (5 Tests)
# ----------------------------------------------------------
start_idx <- W_roll + 2
valid_idx <- start_idx:length(R_vec)
actual_bt <- R_vec[valid_idx]

models_list <- c("HS", "EWMA", "DCC")
bt_results  <- data.frame()

cat("Running Backtests...\n")

for(mod in models_list){
  for(a in alphas){
    var_vec <- Risk_TS[valid_idx, paste0(mod, "_VaR_", a)]
    es_vec  <- Risk_TS[valid_idx, paste0(mod, "_ES_", a)]
    
    # 1. Kupiec (Unconditional Coverage)
    vt <- VaRTest(alpha=a, actual=actual_bt, VaR=var_vec, conf.level=0.95)
    p_kupiec <- vt$uc.LRp
    
    # 2. Christoffersen 1 (Independence)
    # LR_ind = LR_cc - LR_uc (df=1)
    stat_ind <- max(0, vt$cc.LRstat - vt$uc.LRstat)
    p_ind    <- pchisq(stat_ind, df=1, lower.tail=FALSE)
    
    # 3. Christoffersen 2 (Conditional Coverage)
    # LR_cc tests both coverage and independence (df=2)
    p_cc     <- vt$cc.LRp
    
    # 4. Pelletier (Duration)
    vd <- VaRDurTest(alpha=a, actual=actual_bt, VaR=var_vec, conf.level=0.95)
    p_pell <- vd$LRp
    
    # 5. Frey-McNeil (ES Backtest)
    # Test if mean of (Exceedance - ES) is zero
    exceed_idx <- which(actual_bt < var_vec)
    n_exceed   <- length(exceed_idx)
    
    if(n_exceed >= 3) { # 需要至少3个点才能做t-test
      # 残差: 实际值 - ES预测值 (理论上均值应为0)
      resid_es <- actual_bt[exceed_idx] - es_vec[exceed_idx]
      # 双侧 t-test, H0: mean = 0
      fm_test <- t.test(as.numeric(resid_es), mu=0)
      p_fm    <- fm_test$p.value
    } else {
      p_fm <- NA
    }
    
    # 汇总
    # 判定函数 (Pass if p > 0.05)
    pass <- function(p) {
      if(is.na(p)) return("N/A")
      if(p > 0.05) return("Pass") else return("Fail")
    }
    
    bt_results <- rbind(bt_results, data.frame(
      Model = mod,
      Alpha = paste0(a*100, "%"),
      Exceed = vt$actual.exceed,
      
      # Test 1: Kupiec
      Kupiec_p = p_kupiec,
      Res_Kup = pass(p_kupiec),
      
      # Test 2: Christ 1 (Ind)
      Christ1_Ind_p = p_ind,
      Res_Ind = pass(p_ind),
      
      # Test 3: Christ 2 (CC)
      Christ2_CC_p = p_cc,
      Res_CC = pass(p_cc),
      
      # Test 4: Pelletier
      Pell_p = p_pell,
      Res_Pell = pass(p_pell),
      
      # Test 5: Frey-McNeil
      FM_p = p_fm,
      Res_FM = pass(p_fm)
    ))
  }
}
# ----------------------------------------------------------
# 2.5) OFFICIAL Basel Traffic Light Test
# Rolling 250-day window, 99% VaR ONLY
# ----------------------------------------------------------

# Basel 红绿灯判定函数（官方定义）
basel_traffic <- function(exceed){
  if(exceed <= 4) return("Green")
  if(exceed <= 9) return("Yellow")
  return("Red")
}

# Rolling Basel 主函数
rolling_basel_test <- function(actual, VaR, window = 250){
  
  stopifnot(length(actual) == length(VaR))
  
  # 违约指示变量
  exceed_ind <- as.numeric(actual < VaR)
  
  # rolling 250 天内违约次数
  roll_exceed <- zoo::rollapply(
    exceed_ind,
    width = window,
    FUN   = sum,
    align = "right",
    fill  = NA
  )
  
  # 每个窗口对应的 Basel 颜色
  roll_color <- sapply(roll_exceed, function(x){
    if(is.na(x)) return(NA)
    basel_traffic(x)
  })
  
  data.frame(
    Exceed_250 = roll_exceed,
    Basel_250  = roll_color
  )
}

# ----------------------------------------------------------
# Apply rolling Basel to each model (1% VaR)
# ----------------------------------------------------------

basel_summary <- data.frame()

for(mod in models_list){
  
  # 只取 99% VaR
  var_1pct <- Risk_TS[valid_idx, paste0(mod, "_VaR_0.01")]
  
  rb <- rolling_basel_test(
    actual = actual_bt,
    VaR    = var_1pct,
    window = 250
  )
  
  # 去掉前 249 个 NA
  rb_valid <- rb[!is.na(rb$Basel_250), ]
  
  # 各颜色占比
  prop_tab <- prop.table(table(rb_valid$Basel_250))
  
  basel_summary <- rbind(basel_summary, data.frame(
    Model = mod,
    
    # 监管最关心：最坏窗口
    Worst_Window = ifelse("Red" %in% rb_valid$Basel_250, "Red",
                          ifelse("Yellow" %in% rb_valid$Basel_250,
                                 "Yellow", "Green")),
    
    Green_Pct  = ifelse("Green"  %in% names(prop_tab), prop_tab["Green"],  0),
    Yellow_Pct = ifelse("Yellow" %in% names(prop_tab), prop_tab["Yellow"], 0),
    Red_Pct    = ifelse("Red"    %in% names(prop_tab), prop_tab["Red"],    0)
  ))
}

# 输出官方 Basel 结果
cat("\n=== Official Basel Traffic Light Test (Rolling 250 days, 99% VaR) ===\n")
print(knitr::kable(
  basel_summary,
  digits = 3,
  align  = "c",
  caption = "Basel Traffic Light Test based on Rolling 250-day Windows"
))

# ----------------------------------------------------------
# 3) 输出最终表格 (Enhanced)
# ----------------------------------------------------------
cat("\n=== (e) Full Backtesting Results ===\n")
cat("H0: Model is Correct (Pass if p > 0.05)\n")

# 计算预期违约数 (Expected Exceedances)
# actual_bt 是回测样本，N_obs = length(actual_bt)
N_obs <- length(actual_bt)
bt_results$Exp_Exc <- N_obs * as.numeric(gsub("%", "", bt_results$Alpha))/100

# 重新排列列顺序，把 Expected 放进去
display_tab <- bt_results[, c("Model", "Alpha", "Exp_Exc", "Exceed", 
                              "Kupiec_p", "Res_Kup", 
                              "Christ1_Ind_p", "Res_Ind", 
                              "Christ2_CC_p", "Res_CC", 
                              "Pell_p", "Res_Pell", 
                              "FM_p", "Res_FM")]


colnames(display_tab) <- c("Model", "Alpha", "Exp", "Act",
                           "Kupiec", "Res",
                           "Ind", "Res",
                           "CC", "Res",
                           "Pell", "Res",
                           "FM", "Res")


print(knitr::kable(display_tab, digits=3, align="c", caption="Backtesting: Expected(Exp) vs Actual(Act) Exceedances"))

# ----------------------------------------------------------
# 3. Final Basel II "Traffic Light" Report
# (Based on the exact table in your image)
# ----------------------------------------------------------

# 1. 定义图中规则的查找函数
get_basel_params <- function(exceedances) {
  # 绿色区域 (Green Zone)
  if (exceedances <= 4) {
    return(list(Zone = "Green", k = 0.00))
  } 
  # 黄色区域 (Yellow Zone) - 精确匹配图中的数值
  else if (exceedances == 5) { return(list(Zone = "Yellow", k = 0.40)) }
  else if (exceedances == 6) { return(list(Zone = "Yellow", k = 0.50)) }
  else if (exceedances == 7) { return(list(Zone = "Yellow", k = 0.65)) }
  else if (exceedances == 8) { return(list(Zone = "Yellow", k = 0.75)) }
  else if (exceedances == 9) { return(list(Zone = "Yellow", k = 0.85)) }
  # 红色区域 (Red Zone)
  else {
    return(list(Zone = "Red", k = 1.00))
  }
}

# 2. 执行计算 (只针对最后 250 天，符合监管标准)
basel_final_report <- data.frame()
window_size <- 250  # 监管要求的窗口大小

cat(sprintf("Generating Basel Report based on the LAST %d days...\n", window_size))

for(mod in models_list){
  
  # 1. 获取 99% VaR 和 实际收益率
  full_var <- Risk_TS[valid_idx, paste0(mod, "_VaR_0.01")]
  full_act <- actual_bt
  
  # 2. 截取最后 250 个观测值 (监管评估窗口)
  # 注意：如果数据不足 250 天，这部分代码会报错或取全部数据
  n_total <- length(full_act)
  if(n_total < window_size) stop("Error: Backtesting sample is smaller than 250 days!")
  
  recent_act <- tail(full_act, window_size)
  recent_var <- tail(full_var, window_size)
  
  # 3. 计算违约次数 (Exceedances)
  n_exceed <- sum(recent_act < recent_var)
  
  # 4. 获取对应的 Zone 和 k (查表)
  params <- get_basel_params(n_exceed)
  
  # 5. 计算累积概率 (cdf) - 对应图中最后一列
  # 这是二项分布概率：在 99% VaR 下，250 天内发生 n_exceed 次或更少违约的概率
  cdf_val <- pbinom(n_exceed, size = window_size, prob = 0.01)
  
  # 6. 组装结果
  basel_final_report <- rbind(basel_final_report, data.frame(
    Model = mod,
    Exceedances = n_exceed,
    Zone = params$Zone,
    Plus_Factor_k = params$k,
    Total_Multiplier = 3 + params$k,  # 监管资本乘数 = 3 (最低) + k
    cdf_percent = paste0(round(cdf_val * 100, 2), "%") # 对应图中的 cdf 列
  ))
}

# ----------------------------------------------------------
# 3. 输出表格
# ----------------------------------------------------------
cat("\n=== Basel Committee 'Traffic Lights' Approach (Last 250 Days) ===\n")
print(knitr::kable(
  basel_final_report, 
  digits = 2, 
  align = "c",
  caption = "Final Regulatory Capital Charge Determination"
))

# ----------------------------------------------------------
# 4) Official Basel Traffic Light Results (Rolling 250 days)
# ----------------------------------------------------------

cat("\n=== Official Basel Traffic Light Test (99% VaR, Rolling 250 days) ===\n")

print(knitr::kable(
  basel_summary,
  digits = 3,
  align  = "c",
  caption = "Basel Traffic Light Test based on Rolling 250-day Windows"
))

library(ggplot2)
library(zoo)
library(reshape2) # 用于数据变形

# ----------------------------------------------------------
# 5) Visualization: Basel Traffic Light Timeline
# ----------------------------------------------------------

# 1. 计算滚动违约次数 (Rolling 250-day Sum of Exceedances)
# ----------------------------------------------------------
window_size <- 250
plot_data <- data.frame(Date = index(rP)[valid_idx]) # 时间轴

# 循环计算每个模型的滚动违约数
for(mod in c("HS", "EWMA", "DCC")){
  # 获取 1% VaR 数据
  var_col <- paste0(mod, "_VaR_0.01")
  var_vec <- Risk_TS[valid_idx, var_col]
  
  # 判断每天是否违约 (0 或 1)
  hit_vec <- as.numeric(actual_bt < var_vec)
  
  # 计算滚动求和 (Rolling Sum)
  # align="right" 表示今天的值代表过去250天的总和
  roll_sum <- zoo::rollsum(hit_vec, k = window_size, fill = NA, align = "right")
  
  plot_data[[mod]] <- roll_sum
}

# 去掉开头的 NA 值 (前249天无法计算)
plot_data <- na.omit(plot_data)

# 2. 转换为长格式 (Long Format) 以便 ggplot 绘图
# ----------------------------------------------------------
plot_df_long <- reshape2::melt(plot_data, id.vars = "Date", 
                               variable.name = "Model", 
                               value.name = "Exceedances")

# 3. 绘图 (这是核心部分)
# ----------------------------------------------------------
p_basel <- ggplot(plot_df_long, aes(x = Date, y = Exceedances)) +
  
  # --- A. 绘制巴塞尔红绿灯背景区域 (Zones) ---
  # 绿色区域 (0-4): 安全
  annotate("rect", xmin = min(plot_data$Date), xmax = max(plot_data$Date), 
           ymin = 0, ymax = 4.5, fill = "darkgreen", alpha = 0.15) +
  # 黄色区域 (5-9): 警告
  annotate("rect", xmin = min(plot_data$Date), xmax = max(plot_data$Date), 
           ymin = 4.5, ymax = 9.5, fill = "gold", alpha = 0.25) +
  # 红色区域 (10+): 违规/罚款
  annotate("rect", xmin = min(plot_data$Date), xmax = max(plot_data$Date), 
           ymin = 9.5, ymax = max(plot_df_long$Exceedances, 20), fill = "darkred", alpha = 0.2) +
  
  # --- B. 添加分界线 ---
  # 把 size 改为 linewidth
  geom_hline(yintercept = 4.5, linetype = "dashed", color = "darkgreen", linewidth=0.3) +
  geom_hline(yintercept = 9.5, linetype = "dashed", color = "red", linewidth=0.3) +
  
  # --- C. 绘制模型曲线 ---
  # 把 size 改为 linewidth
  geom_line(aes(color = Model), linewidth = 1) +
  
  # --- D. 标注与美化 ---
  scale_color_manual(values = c("HS" = "black", "EWMA" = "blue", "DCC" = "red")) +
  labs(title = "Basel II Traffic Light Timeline (Rolling 250-day Window)",
       subtitle = "Green Zone: 0-4 | Yellow Zone: 5-9 | Red Zone: 10+ Exceedances",
       y = "Number of Exceedances (Last 250 days)",
       x = "") +
  theme_minimal() +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()) +
  
  # 强制 Y 轴显示整数刻度
  scale_y_continuous(breaks = seq(0, max(plot_df_long$Exceedances, 15), by = 2))

# 打印图表
print(p_basel)
# ----------------------------------------------------------
# [Recommended] Plot for DCC 1% (Matches Basel Traffic Light)
# ----------------------------------------------------------
# 1. 提取 1% 的数据
var_1pct <- Risk_TS[valid_idx, "DCC_VaR_0.01"] # 确保列名对应你之前生成的Risk_TS

# 2. 绘图设置
par(mfrow=c(1,1), mar=c(4,4,3,1))

# A. 画背景（实际收益率）
plot(index(rP)[valid_idx], actual_bt, 
     type="p", pch=16, cex=0.4, col="gray70", # 使用稍深一点的灰色
     main="Basel Backtesting Visual: DCC-GARCH (VaR 1%)", 
     ylab="Daily Returns", xlab="Time")

# B. 画风险线 (VaR Threshold)
lines(index(rP)[valid_idx], var_1pct, col="blue", lwd=1.5) # 蓝色线表示防线

# C. 标出违约点 (Exceedances) - 即“红灯”时刻
exc_pts <- which(actual_bt < var_1pct)
points(index(rP)[valid_idx][exc_pts], actual_bt[exc_pts], 
       pch=18, col="red", cex=1.2) # pch=18 是实心菱形，看起来更像那张参考图

# D. 添加图例
legend("bottomright", 
       legend=c("Actual Returns", "VaR Limit (1%)", "Violations (Exceedances)"), 
       col=c("gray70", "blue", "red"), 
       pch=c(16, NA, 18), 
       lty=c(NA, 1, NA), 
       pt.cex=c(1, 1, 1.2),
       bty="n", cex=0.8)

# E. (可选) 打印违约的具体日期，增加报告深度
cat("The 6 days where the model failed:\n")
print(index(rP)[valid_idx][exc_pts])


