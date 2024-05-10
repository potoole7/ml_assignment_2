#### Regularisation analysis of Ionosphere dataset ####

# Connectivity of graphs
# https://cedar.wwu.edu/cgi/viewcontent.cgi?article=1082&context=math_facpubs
# TODO: Change to for loop (done)
# TODO: Add each value to dataframe (done)
# TODO: Functionalise everything to share with testing (done)
# TODO: Average over 10 random splits (done)
# TODO: Allow use of exponential smooth function (done, but doesn't seem to work!)
# TODO: Calculate for test data (won't work for just 6 neighbours I'd imagine!)
# TODO: Check lambda1 for S, should be the same as in paper (34.9907) (far smaller! What if I use 34.9907?)
# Even max isn't that large!
# TODO: Generalisation bound calculation wrong, need to fix somehow
# TODO: Also test with exponential smooth function!


#### Metadata ####

label_sizes <- c(10, 20, 40, 60, 80, 100, 200, 300) # size of label sets
gammas <- c(0.001, 0.01, 0.1, 1) # regularisation parameters
nn <- 6 # nearest neighbours to take for adjacency matrix
seed <- 123
n_rand_split <- 10 # number of random splits to average over

#### Libs ####

source("R/functions.R")
library(igraph)


#### Load Data ###

# load ionosphere data 
iono <- evclass::ionosphere

# rescale labels
iono$y <- ifelse(iono$y == 1, 1, -1)  


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
   rand_choice <- sample(1:nrow(iono$x), size = label_size)
   # part of code that must change for (i) random splits and (ii) test data
   # Testing on training data, so make training data the same as test data
   xl <- xu <- iono$x[rand_choice, ]
   yl <- yu <- iono$y[rand_choice]
   # xu for generalisation bounds
   xugen <- iono$x[-rand_choice, ]
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

readr::write_csv(train_df, "iono_train_df.csv")


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
   rand_choice <- sample(1:nrow(iono$x), size = label_size)
   # part of code that must change for (i) random splits and (ii) test data
   # Testing on training data, so make training data the same as test data
   # xl <- xu <- iono$x[rand_choice, ]
   xl <- iono$x[rand_choice, ]
   xu <- iono$x[-rand_choice, ]
   # yl <- yu <- iono$y[rand_choice]
   yl <- iono$y[rand_choice]
   yu <- iono$y[-rand_choice]
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
