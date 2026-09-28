chick_sites <- function(K = 12) {
  cw <- as.data.frame(ChickWeight)
  cw$Chick <- factor(as.character(cw$Chick))
  cw$site <- as.integer(cw$Chick) %% K
  cw
}

logistic_sites <- function(K = 10, n_pat = 20, seed = 1) {
  set.seed(seed)
  lapply(seq_len(K), function(k) {
    id <- rep(seq_len(n_pat), each = 3)
    x <- rnorm(length(id))
    u <- rnorm(n_pat, 0, 0.5)[id]
    data.frame(
      pat_id = id, x = x, z = rbinom(length(id), 1, 0.4),
      y = rbinom(length(id), 1, plogis(-0.5 + 0.8 * x + u))
    )
  })
}
