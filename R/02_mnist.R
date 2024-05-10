#### Regularisation analysis of MNIST dataset ####

#### Metadata ####

# mnist_size <- 5000 # n rows in mnist to (randomly) take; too large otherwise
mnist_size <- 3000
label_sizes <- c(20, 40, 80, 100, 200, 400, 800, 2000) # size of label sets
# label_sizes <- c(10, 20, 40, 60, 80, 100, 200, 300) 
gammas <- c(0.001, 0.01, 0.1, 1) # regularisation parameters
nn <- 6 # nearest neighbours to take for adjacency matrix
seed <- 123
n_rand_split <- 10 # number of random splits to average over

#### Libs ####

source("R/functions.R")
library(igraph)


#### Load Data ###

# load mnistsphere data 
# mnist <- evclass::mnistsphere
mnist <- readRDS(file = "assignment2/data/mnist.RDS")

# join data
mnist <- rbind(
  cbind(mnist$train$images, mnist$train$labels),
  cbind(mnist$test$images, mnist$test$labels)
)

# limit analysis to 7 and 8s (i.e. two-class case)
mnist <- mnist[mnist[, ncol(mnist)] %in% c(7, 8), ]

# randomly take 1000 7s and 1000 8s, data too large otherwise!
set.seed(seed)

mnist <- mnist[
  as.vector(vapply(c(7, 8), \(x) {
    sample(which(mnist[, ncol(mnist)] == x), size = mnist_size, replace = FALSE)
  }, numeric(mnist_size)))
, ]

# Change labels to match -1, 1 structure
mnist[, ncol(mnist)] <- ifelse(mnist[, ncol(mnist)] == 7, -1, 1)

# split again to match ionosphere
mnist <- list(
  "x" = mnist[, -ncol(mnist)], 
  "y" = mnist[, ncol(mnist)]
)


#### Training Errors ####

train_df <- tidyr::crossing(
  "L"           = label_sizes, 
  "gamma"       = c(0, gammas),  # TODO: double check order is right here
  "class_error" = NA,
  "gen_bound"   = NA
)

for (label_size in label_sizes) {
  print(paste("label set of size =", label_size))
  set.seed(seed)
  # misclassification rates and gen bounds to average over for random splits
  misclass_rates <- gen_bounds <- rep(0, (length(gammas) + 1))
  # ignore messages
  sink("/dev/null")
  for (i in seq_len(n_rand_split)) {
   rand_choice <- sample(1:nrow(mnist$x), size = label_size)
   # part of code that must change for (i) random splits and (ii) test data
   # Testing on training data, so make training data the same as test data
   xl <- xu <- mnist$x[rand_choice, ]
   yl <- yu <- mnist$y[rand_choice]
   xugen <- mnist$x[-rand_choice, ]
   # calculate misclassification and gen bounds
   misclass_gen <- calc_tables(xl, xu, xugen, yl, yu, nn, reg = FALSE) 
   # add to sum across all random splits
   misclass_rates <- misclass_rates + misclass_gen$misclass_rates
   gen_bounds <- gen_bounds + misclass_gen$gen_bounds
  }
  # Close connection
  sink() 
  
  # take average across all random splits
  misclass_rates <- misclass_rates / n_rand_split
  gen_bounds <- gen_bounds / n_rand_split
  
  # Input misclassification and gen bounds for this label size
  train_df[train_df$L == label_size, ]$class_error <- misclass_rates
  train_df[train_df$L == label_size, ]$gen_bound <- gen_bounds
}

readr::write_csv(train_df, "mnist_train_df.csv")


#### Test Errors ####

test_df <- tidyr::crossing(
  "L"           = label_sizes, 
  "gamma"       = c(0, gammas),  # TODO: double check order is right here
  "class_error" = NA,
  "gen_bound"   = NA
)

for (label_size in label_sizes) {
  print(paste("label set of size =", label_size))
  set.seed(seed)
  # misclassification rates and gen bounds to average over for random splits
  misclass_rates <- gen_bounds <- rep(0, (length(gammas) + 1))
  # ignore messages
  sink("/dev/null")
  for (i in seq_len(n_rand_split)) {
   rand_choice <- sample(1:nrow(mnist$x), size = label_size)
   # part of code that must change for (i) random splits and (ii) test data
   # Testing on training data, so make training data the same as test data
   # xl <- xu <- mnist$x[rand_choice, ]
   xl <- mnist$x[rand_choice, ]
   xu <- mnist$x[-rand_choice, ]
   # yl <- yu <- mnist$y[rand_choice]
   yl <- mnist$y[rand_choice]
   yu <- mnist$y[-rand_choice]
   # calculate misclassification and gen bounds
   misclass_gen <- calc_tables(xl, xu, yl, yu, nn) 
   # add to sum across all random splits
   misclass_rates <- misclass_rates + misclass_gen$misclass_rates
   gen_bounds <- gen_bounds + misclass_gen$gen_bounds
  }
  # Close connection
  sink() 
  
  # take average across all random splits
  misclass_rates <- misclass_rates / n_rand_split
  gen_bounds <- gen_bounds / n_rand_split
  
  # Input misclassification and gen bounds for this label size
  train_df[train_df$L == label_size, ]$class_error <- misclass_rates
  train_df[train_df$L == label_size, ]$gen_bound <- gen_bounds
}
