# ---------------------------------------------------------------------------
# 04_diffusion.R - value-added analysis: does violence spread to neighbours?
#
# A spatial Markov analysis (Rey 2001) of the quarterly cube. Each
# township-quarter is Quiet or Active; transitions are tabulated conditional on
# how many of the township's neighbours were Active in the same quarter.
# If violence diffuses by contagion, a quiet township's chance of becoming
# active should rise with the activity around it. If it simply relocates, or if
# activity is driven by fixed regional conditions, the neighbourhood should add
# little once region and time are controlled.
#
# Output: data/derived/diffusion.rds
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(sf); library(dplyr); library(tidyr); library(purrr)
})
source("R/00_config.R")
P <- readRDS(file.path(DIR_DERIVED, "prepared.rds"))
L <- readRDS(file.path(DIR_DERIVED, "lmsa.rds"))
tsp <- P$tsp; nb <- L$nb
n <- nrow(tsp); TT <- length(P$quarters)

# counts as an n x T matrix, rows in tsp order
M <- P$cube |> select(TS_PCODE, period, pv) |>
  pivot_wider(names_from = period, values_from = pv) |>
  slice(match(tsp$TS_PCODE, TS_PCODE)) |> select(-TS_PCODE) |> as.matrix()

ctx_lab <- c("No neighbour active", "Some neighbours active", "Most neighbours active")

# Build the transition table for one activity threshold and one neighbour list.
transitions <- function(M, nb, thr) {
  A <- M >= thr
  share <- vapply(seq_len(TT), function(t) vapply(nb, function(j) mean(A[j, t]), 1), numeric(n))
  ctx <- matrix(cut(share, c(-Inf, 0, 0.5 - 1e-9, Inf), labels = ctx_lab), n, TT)
  # own history: active in any of the previous four quarters (NA early on)
  hist4 <- matrix(NA, n, TT)
  for (t in 5:TT) hist4[, t] <- rowSums(A[, (t - 4):(t - 1), drop = FALSE]) > 0
  expand_grid(i = seq_len(n), t = seq_len(TT - 1)) |>
    mutate(from = A[cbind(i, t)], to = A[cbind(i, t + 1)],
           ctx = factor(ctx[cbind(i, t)], levels = ctx_lab),
           share = share[cbind(i, t)], hist4 = hist4[cbind(i, t)],
           ST = tsp$ST[i], TS_PCODE = tsp$TS_PCODE[i], quarter = P$quarters[t])
}

# Conditional transition probabilities and the likelihood-ratio test of
# homogeneity across neighbourhood contexts (Bickenbach and Bode 2003).
markov_tab <- function(tr) {
  tr |> group_by(from, ctx) |>
    summarise(n = n(), p_active_next = mean(to), .groups = "drop") |>
    mutate(from = if_else(from, "Active", "Quiet"))
}
lr_test <- function(tr) {
  G2 <- 0; df <- 0
  for (f in c(FALSE, TRUE)) {
    d <- filter(tr, from == f)
    pooled <- mean(d$to)
    by <- d |> group_by(ctx) |> summarise(k = sum(to), m = n(), .groups = "drop") |> filter(m > 0)
    ll <- function(k, m, p) ifelse(k > 0, k * log(p), 0) + ifelse(m - k > 0, (m - k) * log(1 - p), 0)
    G2 <- G2 + 2 * sum(ll(by$k, by$m, by$k / by$m) - ll(by$k, by$m, pooled))
    df <- df + nrow(by) - 1
  }
  tibble(G2 = G2, df = df, p = pchisq(G2, df, lower.tail = FALSE))
}

tr <- transitions(M, nb, ACTIVE_MIN)
tab <- markov_tab(tr)
lr <- lr_test(tr)
onset <- function(tr, c) mean(tr$to[!tr$from & tr$ctx == c])
gap_obs <- onset(tr, ctx_lab[3]) - onset(tr, ctx_lab[1])

## Spatial permutation test ------------------------------------------------------
# The LR test treats every transition as independent, which they are not: the
# same township appears 18 times and neighbouring townships share conditions.
# This test keeps every township's full time series intact and shuffles which
# location it sits at, so temporal dependence is preserved and only the spatial
# arrangement is destroyed. The statistic is the onset-rate gap between
# townships with most neighbours active and those with none.
set.seed(SEED)
NPERM <- 499L
gap_perm <- replicate(NPERM, {
  Mp <- M[sample(n), , drop = FALSE]
  trp <- transitions(Mp, nb, ACTIVE_MIN)
  onset(trp, ctx_lab[3]) - onset(trp, ctx_lab[1])
})
p_perm <- (sum(gap_perm >= gap_obs) + 1) / (NPERM + 1)

## Controlled model -----------------------------------------------------------------
# Onset among quiet township-quarters, with state/region and quarter fixed
# effects and the township's own recent history. If the neighbour effect were
# only "this is a violent region" or "this was a violent quarter", it would
# shrink towards zero here.
onset_df <- tr |> filter(!from, !is.na(hist4)) |>
  mutate(ctx = relevel(ctx, ctx_lab[1]), qf = factor(quarter))
fit <- glm(to ~ ctx + hist4 + ST + qf, family = binomial, data = onset_df)
fit0 <- glm(to ~ hist4 + ST + qf, family = binomial, data = onset_df)
or_tab <- tibble(term = names(coef(fit)), est = coef(fit)) |>
  filter(grepl("^ctx|hist4", term))
# township-cluster bootstrap for intervals that respect repeated observations
set.seed(SEED)
boot <- replicate(200, {
  ids <- sample(unique(onset_df$TS_PCODE), replace = TRUE)
  d <- bind_rows(lapply(ids, \(id) onset_df[onset_df$TS_PCODE == id, ]))
  f <- try(glm(to ~ ctx + hist4 + ST + qf, family = binomial, data = d), silent = TRUE)
  if (inherits(f, "try-error")) rep(NA, nrow(or_tab)) else coef(f)[or_tab$term]
})
or_tab <- or_tab |> mutate(
  OR = exp(est),
  lo = exp(apply(boot, 1, quantile, 0.025, na.rm = TRUE)),
  hi = exp(apply(boot, 1, quantile, 0.975, na.rm = TRUE)))
lr_ctx <- anova(fit0, fit, test = "LRT")

## Out-of-sample check ------------------------------------------------------------------
# Fit on transitions up to the one ending 2024 Q3; score the four later ones.
# The question is whether knowing the neighbourhood improves on knowing only a
# township's own history, which is the baseline any monitoring team already has.
cut_q <- as.Date("2024-07-01")
trn <- filter(onset_df, quarter < cut_q); tst <- filter(onset_df, quarter >= cut_q)
f_full <- glm(to ~ ctx + hist4 + ST, family = binomial, data = trn)
f_own  <- glm(to ~ hist4 + ST, family = binomial, data = trn)
auc <- function(y, s) {
  r <- rank(s); n1 <- sum(y); n0 <- sum(!y)
  (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
tst$p_full <- predict(f_full, tst, type = "response")
tst$p_own  <- predict(f_own, tst, type = "response")
oos <- tibble(model = c("Own history and state", "Plus neighbourhood"),
              auc = c(auc(tst$to, tst$p_own), auc(tst$to, tst$p_full)),
              n_test = nrow(tst), onsets = sum(tst$to))

## Stability across Operation 1027 ------------------------------------------------------
# Operation 1027 began on 27 Oct 2023, inside 2023 Q4. Transitions that start
# in 2023 Q4 or later are "after". A stationary chain should give similar
# conditional probabilities in both phases.
op_q <- as.Date("2023-10-01")
phase_tab <- tr |>
  mutate(phase = factor(if_else(quarter < op_q, "Before Operation 1027", "From Operation 1027"),
                        levels = c("Before Operation 1027", "From Operation 1027"))) |>
  group_by(phase) |> group_modify(\(d, k) markov_tab(d)) |> ungroup()

## Threshold sensitivity --------------------------------------------------------------------
thr_sens <- map_dfr(c(ACTIVE_SENS[1], ACTIVE_MIN, ACTIVE_SENS[2]), function(th) {
  t2 <- transitions(M, nb, th)
  tibble(threshold = th,
         none = onset(t2, ctx_lab[1]), some = onset(t2, ctx_lab[2]), most = onset(t2, ctx_lab[3]),
         ratio = most / none, share_active = mean(M >= th))
})

## Townships to watch ---------------------------------------------------------------------------
# Quiet in the last observed quarter (2025 Q3) with the model's estimated chance
# of turning active in 2025 Q4. In-sample estimate, not a validated forecast
# beyond the out-of-sample check above.
A <- M >= ACTIVE_MIN
share_last <- vapply(nb, function(j) mean(A[j, TT]), 1)
last <- tibble(i = seq_len(n), from = A[, TT], share = share_last,
               ctx = factor(as.character(cut(share_last, c(-Inf, 0, 0.5 - 1e-9, Inf), labels = ctx_lab)),
                            levels = ctx_lab),
               hist4 = rowSums(A[, (TT - 4):(TT - 1), drop = FALSE]) > 0,
               ST = tsp$ST,
               qf = factor(as.character(max(onset_df$quarter)), levels = levels(onset_df$qf)))
q_new <- filter(last, !from)
q_new$p_onset <- predict(fit, newdata = q_new, type = "response")
watch <- q_new |>
  mutate(TS_PCODE = tsp$TS_PCODE[i], TS = tsp$TS[i], pv_q = M[cbind(i, TT)]) |>
  select(TS_PCODE, TS, ST, ctx, share, hist4, pv_q, p_onset) |>
  arrange(desc(p_onset))

saveRDS(list(tr = tr, tab = tab, lr = lr, gap_obs = gap_obs, gap_perm = gap_perm,
             p_perm = p_perm, or_tab = or_tab, lr_ctx = lr_ctx, oos = oos,
             phase_tab = phase_tab, thr_sens = thr_sens, watch = watch, ctx_lab = ctx_lab),
        file.path(DIR_DERIVED, "diffusion.rds"))

print(tab); print(lr); cat("gap", gap_obs, " perm p", p_perm, " perm mean", mean(gap_perm), "\n")
print(or_tab); print(lr_ctx); print(oos); print(phase_tab, n = 20); print(thr_sens)
print(head(watch, 15)); cat("watch n:", nrow(watch), "\n")
