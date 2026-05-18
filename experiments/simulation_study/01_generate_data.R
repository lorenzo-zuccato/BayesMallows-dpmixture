library(BayesMallowsDPMixture)

# ── Helper functions ──────────────────────────────────────────────────────────

make_preferences <- function(rankings, lambda) {
  n_assessors <- nrow(rankings)
  n_items     <- ncol(rankings)
  x           <- rpois(n_assessors, lambda)
  x           <- replace(x, which(x == 0), 1)
  
  preferences <- data.frame(
    assessor    = rep(1:n_assessors, x),
    bottom_item = rep(NA, sum(x)),
    top_item    = rep(NA, sum(x))
  )
  
  count <- 1
  for (i in 1:n_assessors) {
    item1 <- sample(1:n_items, x[i], replace = TRUE)
    repeat {
      item2 <- sample(1:n_items, x[i], replace = TRUE)
      if (!any(item1 == item2)) break
    }
    for (j in 1:x[i]) {
      if (rankings[i, item1[j]] < rankings[i, item2[j]]) {
        preferences$bottom_item[count] <- item2[j]
        preferences$top_item[count]    <- item1[j]
      } else {
        preferences$bottom_item[count] <- item1[j]
        preferences$top_item[count]    <- item2[j]
      }
      count <- count + 1
    }
  }
  return(preferences)
}

make_topk <- function(rankings, lambda) {
  n_assessors <- nrow(rankings)
  topk_ranks  <- rankings
  x           <- rpois(n_assessors, lambda)
  for (i in 1:n_assessors) {
    topk_ranks[i, topk_ranks[i, ] > x[i]] <- NA
  }
  return(topk_ranks)
}

# ── True partition parameters ─────────────────────────────────────────────────

rho1   <- 1:30
rho2   <- 30:1
rho3   <- c(16:30, 1:15)
n1     <- 50;  alpha1 <- 3
n2     <- 25;  alpha2 <- 2
n3     <- 25;  alpha3 <- 5

# ── Generate 30 independent datasets ─────────────────────────────────────────

set.seed(42)

generate_dataset <- function() {
  rankings <- matrix(NA, nrow = 100, ncol = 30)
  rankings[1:n1, ]                              <- sample_mallows(rho1, alpha1, n1)
  rankings[(n1 + 1):(n1 + n2), ]               <- sample_mallows(rho2, alpha2, n2)
  rankings[(n1 + n2 + 1):(n1 + n2 + n3), ]     <- sample_mallows(rho3, alpha3, n3)
  list(
    top8   = make_topk(rankings, 8),
    pref30 = make_preferences(rankings, 30)
  )
}

datasets <- lapply(1:30, function(i) generate_dataset())

# ── Collect into structured lists ─────────────────────────────────────────────

parameters <- list(
  rho1   = rho1,   rho2   = rho2,   rho3   = rho3,
  alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3,
  n1     = n1,     n2     = n2,     n3     = n3
)

data_top8 <- list(
  parameters = parameters,
  datasets   = lapply(datasets, function(d) d$top8)
)

data_pref30 <- list(
  parameters = parameters,
  datasets   = lapply(datasets, function(d) d$pref30)
)

# ── Save ──────────────────────────────────────────────────────────────────────

save(data_top8,   file = "data/data_top8_new.RData")
save(data_pref30, file = "data/data_pref30_new.RData")