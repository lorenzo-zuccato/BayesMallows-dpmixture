# 01_preprocess_data.R
#
# This script preprocesses the MovieLens 1M dataset to obtain the pairwise
# preference data used in the paper. The raw data must be downloaded from:
# https://grouplens.org/datasets/movielens/1m/
# and the files ratings.dat and users.dat placed in the data/ folder.

library(entropy)
library(dplyr)

set.seed(1234)

n_assessors <- 300
n_movies    <- 50

# ── Load raw data ─────────────────────────────────────────────────────────────

ratings_raw <- readLines("data/ratings.dat")
ratings <- do.call(rbind, strsplit(ratings_raw, "::"))
ratings <- as.data.frame(ratings, stringsAsFactors = FALSE)
ratings <- ratings[, 1:3]
colnames(ratings) <- c("user_id", "movie_id", "rating")
ratings$user_id  <- as.integer(ratings$user_id)
ratings$movie_id <- as.integer(ratings$movie_id)
ratings$rating   <- as.integer(ratings$rating)

users_raw <- readLines("data/users.dat")
users <- do.call(rbind, strsplit(users_raw, "::"))
users <- as.data.frame(users, stringsAsFactors = FALSE)
users <- users[, 1:4]
colnames(users) <- c("user_id", "gender", "age", "occupation")
users$user_id     <- as.integer(users$user_id)
users$age         <- as.integer(users$age)
users$occupation  <- as.integer(users$occupation)

# ── Step 1: Keep only ratings for the 200 most popular movies ─────────────────

popular_movies <- as.integer(names(sort(table(ratings$movie_id), decreasing = TRUE))[1:200])
ratings        <- ratings[ratings$movie_id %in% popular_movies, ]

# ── Step 2: Sample n_assessors random assessors ───────────────────────────────

chosen_assessors <- sample(unique(ratings$user_id), n_assessors, replace = FALSE)
ratings          <- ratings[ratings$user_id %in% chosen_assessors, ]

# ── Step 3: Keep the n_movies movies with highest rating entropy ──────────────

calculate_entropy <- function(ratings) {
  probabilities <- table(ratings) / length(ratings)
  entropy::entropy(probabilities, base = 2)
}

entropy_values      <- tapply(ratings$rating, ratings$movie_id, calculate_entropy)
high_entropy_movies <- as.integer(names(sort(entropy_values, decreasing = TRUE))[1:n_movies])
ratings             <- ratings[ratings$movie_id %in% high_entropy_movies, ]

# ── Step 4: Keep assessors with sufficient rating diversity ───────────────────

ratings <- ratings %>%
  group_by(user_id) %>%
  filter(n_distinct(rating) >= 3 |
           (n_distinct(rating) == 2 & !(any(table(rating) == 1)))) %>%
  ungroup()

users <- users[users$user_id %in% unique(ratings$user_id), ]

cat("Assessors remaining after filtering:", length(unique(ratings$user_id)), "\n")

# ── Step 5: Relabel movies and users ─────────────────────────────────────────

ratings$user_id  <- as.numeric(factor(ratings$user_id,  levels = unique(ratings$user_id)))
ratings$movie_id <- as.numeric(factor(ratings$movie_id, levels = unique(ratings$movie_id)))
users$user_id    <- seq(1, length(unique(ratings$user_id)))

cat("Number of assessors:", max(ratings$user_id), "\n")
cat("Number of movies:   ", max(ratings$movie_id), "\n")

# ── Step 6: Sample one preference per assessor to hold out ───────────────────

ratings_pref_temp <- ratings %>%
  inner_join(ratings, by = "user_id", relationship = "many-to-many") %>%
  filter(rating.x < rating.y) %>%
  select(assessor = user_id, bottom_item = movie_id.x, top_item = movie_id.y) %>%
  arrange(assessor, bottom_item, top_item)

deleted_preferences <- ratings_pref_temp %>%
  group_by(assessor) %>%
  group_split() %>%
  lapply(function(df) df %>% slice_sample(n = 1)) %>%
  bind_rows()

# ── Step 7: Remove one rating per assessor involved in held-out preference ────

which_one <- sample(c(2, 3), length(unique(ratings$user_id)), replace = TRUE)

for (i in 1:length(unique(ratings$user_id))) {
  ratings <- ratings %>%
    filter(!(user_id == i &
               movie_id == as.integer(deleted_preferences[i, which_one[i]])))
}

# ── Step 8: Compute final pairwise preferences ────────────────────────────────

ratings_pref <- ratings %>%
  inner_join(ratings, by = "user_id", relationship = "many-to-many") %>%
  filter(rating.x < rating.y) %>%
  select(assessor = user_id, bottom_item = movie_id.x, top_item = movie_id.y) %>%
  arrange(assessor, bottom_item, top_item)

cat("Average preferences per assessor:",
    round(nrow(ratings_pref) / length(unique(ratings_pref$assessor)), 1), "\n")

# ── Save ──────────────────────────────────────────────────────────────────────

save(ratings_pref, deleted_preferences, users,
     file = "data/data_movielens.RData")

cat("Data saved to data/data_movielens.RData\n")
