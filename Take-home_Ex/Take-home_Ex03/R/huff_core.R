# ---------------------------------------------------------------------------
# huff_core.R - the Huff model, written for fitting
#
# The probability that a resident of block i boards at station j is
#     p_ij = A_j * d_ij^lambda / sum_l (A_l * d_il^lambda)      (Huff 1963)
# summed over the stations l in block i's choice set. Attractiveness is
# multiplicative in the station's characteristics, as in an MCI model:
#     A_j = exp(sum_m gamma_m * x_jm),  x_jm = log(1 + count) or an indicator
# so the utility u_ij = lambda * log(d_ij) + sum_m gamma_m * x_jm is linear in
# the parameters: the Huff model is a conditional logit over stations. Every
# variant fitted here (one lambda, a lambda per region, exponential decay) is
# the same function of a pair-level design matrix Z with u = Z %*% theta.
#
# Predicted boardings are the dwelling units allocated to each station, scaled
# by one trip rate k so that predicted and observed totals agree on the log
# scale:  T_hat_j = k * sum_i du_i * p_ij.
#
# MCI::huff.shares() computes exactly these shares (checked in R/04); this
# version exists because the optimiser evaluates it thousands of times, and it
# supplies the analytic gradient of the loss.
# ---------------------------------------------------------------------------

# Pair-level design for the power-decay model: log distance, then the station
# attraction variables repeated for each pair.
huff_Z <- function(pairs, X) {
  Z <- cbind(lambda = pairs$logd, X[pairs$sj, , drop = FALSE])
  colnames(Z) <- c("lambda", colnames(X))
  Z
}

shares_Z <- function(theta, Z, oi) {
  e <- exp(drop(Z %*% theta))
  e / rowsum(e, oi, reorder = TRUE)[oi]
}

alloc_Z <- function(theta, Z, pairs, n_stn) {
  a <- rowsum(pairs$du * shares_Z(theta, Z, pairs$oi), pairs$sj, reorder = TRUE)
  out <- numeric(n_stn); out[as.integer(rownames(a))] <- a[, 1]
  out
}

# Loss: mean squared error of log boardings over the stations in fit_idx, with
# k solved in closed form. Station ridership spans three orders of magnitude,
# so a 20% miss costs the same at a small station as at a large one. Returns
# the gradient as an attribute.
loss_Z <- function(theta, Z, pairs, obs, fit_idx, grad = TRUE) {
  p <- shares_Z(theta, Z, pairs$oi)
  a <- rowsum(pairs$du * p, pairs$sj, reorder = TRUE)
  Tj <- numeric(length(obs)); Tj[as.integer(rownames(a))] <- a[, 1]
  lT <- log(pmax(Tj, 1e-12))
  F <- fit_idx
  r <- (lT[F] - mean(lT[F])) - (log(obs[F]) - mean(log(obs[F])))
  out <- mean(r^2)
  if (grad) {
    # d p_ij / d theta = p_ij (z_ij - sum_l p_il z_il)
    zbar <- rowsum(p * Z, pairs$oi, reorder = TRUE)[pairs$oi, , drop = FALSE]
    G <- rowsum(pairs$du * p * (Z - zbar), pairs$sj, reorder = TRUE)
    dlT <- matrix(0, length(obs), ncol(Z)); dlT[as.integer(rownames(G)), ] <- G / pmax(Tj[as.integer(rownames(G))], 1e-12)
    # sum(r) = 0, so the centring terms drop out of the gradient
    attr(out, "gradient") <- 2 / length(F) * colSums(r * dlT[F, , drop = FALSE])
  }
  out
}

fit_Z <- function(Z, pairs, obs, fit_idx = seq_along(obs), start = NULL) {
  if (is.null(start)) start <- c(-2, rep(0, ncol(Z) - 1))
  # Fit-index weights allow the bootstrap to repeat stations.
  opt <- optim(start, \(th) as.numeric(loss_Z(th, Z, pairs, obs, fit_idx, grad = FALSE)),
               \(th) attr(loss_Z(th, Z, pairs, obs, fit_idx), "gradient"),
               method = "BFGS", control = list(maxit = 1000, reltol = 1e-12))
  a <- alloc_Z(opt$par, Z, pairs, length(obs))
  k <- exp(mean(log(obs[fit_idx]) - log(pmax(a[fit_idx], 1e-12))))
  list(theta = setNames(opt$par, colnames(Z)), k = k, pred = k * a,
       loss = opt$value, convergence = opt$convergence)
}

# Wrappers for the power-decay model with station matrix X -------------------
huff_shares  <- function(theta, pairs, X) shares_Z(theta, huff_Z(pairs, X), pairs$oi)
huff_predict <- function(theta, pairs, X, n_stn) alloc_Z(theta, huff_Z(pairs, X), pairs, n_stn)
huff_fit <- function(pairs, X, obs, fit_idx = seq_along(obs), start = NULL)
  fit_Z(huff_Z(pairs, X), pairs, obs, fit_idx, start)

# Fit statistics on the log scale and the original scale.
fit_stats <- function(obs, pred) {
  tibble::tibble(
    R2_log = 1 - sum((log(obs) - log(pred))^2) / sum((log(obs) - mean(log(obs)))^2),
    RMSLE = sqrt(mean((log(pred) - log(obs))^2)),
    MdAPE_pct = 100 * median(abs(pred - obs) / obs),
    spearman = cor(obs, pred, method = "spearman"))
}
