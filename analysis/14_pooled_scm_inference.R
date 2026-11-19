# Pooled SCM inference for the 7 treated zones (fit-then-average + permutation).
# Implements:
#   (A) Cavallo, Galiani, Noy & Pantano (2013, REStat): average the per-zone
#       effects, build a placebo distribution by averaging random same-size
#       combinations of control zones, rank the real average against it.
#   (B) same, using the RMSPE-ratio of the group-averaged gap series (shares the
#       aggregate's denominator-shrink caveat, reported for comparison).
#   (C) Dube & Zipperer (2015): average the per-zone percentile ranks.
suppressWarnings(suppressMessages({ library(dplyr); library(tidyr); library(CVXR); library(ggplot2) }))

BUILD <- "analysis"
full_data <- readRDS(file.path(BUILD, "ckpt_step1.rds"))$full_data
ids61     <- readRDS(file.path(BUILD, "zones61_ids.rds"))
treat_ids <- c(75,262,263,140,141,236,237)
treat_date <- as.Date("2017-01-01")

md <- full_data %>% filter(locationID %in% ids61) %>%
  mutate(log_obs = log(observations)) %>% arrange(locationID, month)
bal <- md %>% count(locationID) %>% filter(n==max(n))
md  <- md %>% filter(locationID %in% bal$locationID)
Y_wide <- md %>% select(month, locationID, log_obs) %>%
  pivot_wider(names_from=locationID, values_from=log_obs) %>% arrange(month)
months <- Y_wide$month
Y <- as.matrix(Y_wide %>% select(-month)); colnames(Y) <- colnames(Y_wide)[-1]
treat_col_ids <- as.character(treat_ids)
ctrl_col_ids  <- setdiff(colnames(Y), treat_col_ids)
pre_idx  <- which(months <  treat_date); post_idx <- which(months >= treat_date)

fit_gap <- function(y, Yd) {           # return full monthly gap series
  w  <- Variable(ncol(Yd))
  pr <- Problem(Minimize(sum((y[pre_idx] - Yd[pre_idx,] %*% w)^2)),
                list(w>=0, sum(w)==1))
  wh <- as.numeric(solve(pr)$getValue(w))
  y - as.numeric(Yd %*% wh)
}
rmratio <- function(g) sqrt(mean(g[post_idx]^2)) / sqrt(mean(g[pre_idx]^2))
avgpost <- function(g) mean(g[post_idx])

# ---- treated per-zone fits (donor pool = 54 controls) ----
Yctrl <- Y[, ctrl_col_ids]
treat_gap <- sapply(treat_col_ids, function(tz) fit_gap(Y[,tz], Yctrl))  # 24 x 7

# ---- control placebo fits (each control from the other 53) ----
ctrl_gap <- sapply(ctrl_col_ids, function(c0) {
  donors <- setdiff(ctrl_col_ids, c0)
  fit_gap(Y[, c0], Y[, donors])
})                                                                       # 24 x 54
ctrl_ratio   <- apply(ctrl_gap, 2, rmratio)
ctrl_avgpost <- apply(ctrl_gap, 2, avgpost)

# ---- treated group statistics (fit-then-average) ----
grp_gap   <- rowMeans(treat_gap)                 # averaged gap series
t_avgpost <- avgpost(grp_gap)                    # mean post effect (log pts)
t_ratio_series <- rmratio(grp_gap)               # ratio of the averaged series
t_meanpct <- mean(sapply(seq_len(ncol(treat_gap)),
              function(j) mean(ctrl_ratio < rmratio(treat_gap[,j]))))  # DZ statistic

# ---- permutation over random 7-zone control subsets ----
set.seed(42); n_perm <- 5000; k <- length(treat_col_ids); nc <- length(ctrl_col_ids)
perm_avgpost <- numeric(n_perm); perm_ratio <- numeric(n_perm); perm_meanpct <- numeric(n_perm)
for (b in seq_len(n_perm)) {
  s <- sample.int(nc, k)
  gbar <- rowMeans(ctrl_gap[, s, drop=FALSE])
  perm_avgpost[b] <- avgpost(gbar)
  perm_ratio[b]   <- rmratio(gbar)
  # DZ: each picked zone's percentile among the OTHER controls, then average
  perm_meanpct[b] <- mean(sapply(s, function(j)
                      mean(ctrl_ratio[-j] < ctrl_ratio[j])))
}
p_cavallo_gap   <- (sum(perm_avgpost <= t_avgpost) + 1) / (n_perm + 1)   # effect negative
p_cavallo_ratio <- (sum(perm_ratio   >= t_ratio_series) + 1) / (n_perm + 1)
p_dz_rank       <- (sum(perm_meanpct >= t_meanpct) + 1) / (n_perm + 1)

pct <- function(x) (exp(x)-1)*100
cat("\n================ POOLED SCM INFERENCE (7 treated zones) ================\n")
cat(sprintf("Fit-then-average effect:  %+0.4f log pts  (%+0.1f%%)\n", t_avgpost, pct(t_avgpost)))
cat(sprintf("  (compare: aggregate pool-then-fit = -0.0724 / -7.0%%; per-zone avg = -0.1118 / -10.6%%)\n\n"))
cat(sprintf("(A) Cavallo pooled RAW-GAP permutation p   = %0.4f\n", p_cavallo_gap))
cat(sprintf("    treated mean post gap %+0.4f ranks %d of %d averaged placebos\n",
            t_avgpost, sum(perm_avgpost <= t_avgpost)+1, n_perm+1))
cat(sprintf("(B) Cavallo pooled RMSPE-RATIO permutation p = %0.4f  (ratio of averaged series = %0.1f)\n",
            p_cavallo_ratio, t_ratio_series))
cat(sprintf("(C) Dube-Zipperer mean-percentile-rank p     = %0.4f  (treated mean pct = %0.3f)\n",
            p_dz_rank, t_meanpct))
cat("\nPlacebo raw-gap distribution (mean of random 7-control subsets):\n")
print(round(quantile(perm_avgpost, c(.05,.25,.5,.75,.95)), 4))

# ---- figure: the pooled permutation distribution ----
n_beat <- sum(perm_avgpost <= t_avgpost)
gg <- ggplot(data.frame(avg = perm_avgpost), aes(avg)) +
  geom_histogram(bins = 60, fill = "gray78", color = "white") +
  geom_vline(xintercept = t_avgpost, color = "firebrick", size = 1.1) +
  annotate("text", x = t_avgpost, y = Inf,
           label = sprintf("  treated corridor  %+.1f%%  (p = %.3f)", pct(t_avgpost), p_cavallo_gap),
           hjust = 0, vjust = 1.8, color = "firebrick", size = 4) +
  labs(title = "Pooled permutation test: treated corridor vs. random seven-zone control groups",
       subtitle = sprintf(paste0("Each bar counts random groups of 7 control zones by their average post-2017 gap (%s draws).\n",
                                  "Only %d of %s groups fell as far as the real corridor."),
                          format(n_perm, big.mark=","), n_beat, format(n_perm, big.mark=",")),
       x = "Average post-2017 gap  (log points; more negative = demand fell further below synthetic)",
       y = "number of random control groups") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(BUILD, "scm_pooled_placebo_61.png"), gg, width = 9, height = 5.2, dpi = 130)
cat("wrote analysis/scm_pooled_placebo_61.png\n")

write.csv(data.frame(
  stat = c("fit_then_avg_logpts","fit_then_avg_pct",
           "p_cavallo_rawgap","p_cavallo_ratio","p_dz_rank",
           "ratio_avg_series","treated_mean_pct","n_perm"),
  value = c(t_avgpost, pct(t_avgpost),
            p_cavallo_gap, p_cavallo_ratio, p_dz_rank,
            t_ratio_series, t_meanpct, n_perm)),
  file.path(BUILD, "scm_pooled_inference.csv"), row.names = FALSE)
cat("\nwrote analysis/scm_pooled_inference.csv\n")
