#### Functions for interpolated regularisation ####

#' Regression on graphs (from https://github.com/cran/SSL/blob/5122dffd155a26912e0c484c17f591a31417a155/R/SSL.R)
#' @description \code{sslRegress} develops a regularization framework on graphs.It supports many
#' kinds of distance measurements and graph representations. However, it only supports binary classifications.
#' @param xl a n * p matrix or data.frame of labeled data.
#' @param yl a n * 1 binary labels(1 or -1).
#' @param xu a m * p matrix or data.frame of unlabeled data.
#' @param graph.type character string; which type of graph should be created? Options
#' include\code{tanh} and \code{exp}.
#' \itemize{\item \code{tanh}:tanh-weighted graphs.  \code{w(i,j) = (tanh(alpha1(d(i,j) - alpha2)) + 1)/2}.where \code{d(i,j)} denotes the distance between point i and j. Hyperparameters \code{alpha1} and \code{alpha2} control the slope and cutoff value respectively.
#' \item \code{exp} :exp-weighted graphs.\code{w(i,j) = exp(-d(i,j)^2/alpha^2)},where \code{d(i,j)} denotes the distance between point i and j. Hyperparameter \code{alpha} controls the decay rate.}
#' @param dist.type character string; this parameter controls the type of distance measurement.(see \code{\link{dist}} or \code{\link{pr_DB}}).
#' @param alpha numeric parameter needed when \code{graph.type = exp}
#' @param alpha1 numeric parameter needed when \code{graph.type = tanh}
#' @param alpha2 numeric parameter needed when \code{graph.type = tanh}
#' @param p an integer parameter controls the power of Laplacian for regularization.
#' @param method character string; this parameter choose two possible algorithms:"Tikhonov" means  Tikhonov regularization;"Interpolated" means Interpolated regularization.
#' @param gamma a parameter of Tikhonov regularization.
#' @param knn number of k-nearest neighbours to include in graph, if desired.
#' @return  a m * 1 integer vector representing the predicted labels  of  unlabeled data(1 or -1).
#' @author Junxiang Wang
#' @export
#' @importFrom proxy dist
#' @examples
#' data(iris)
#' xl <- iris[c(1:20, 51:70), -5]
#' xu <- iris[c(21:50, 71:100), -5]
#' yl <- rep(c(1, -1), each = 20)
#' # Tikhonov regularization
#' yu1 <- sslRegress(xl, yl, xu, graph.type = "tanh", alpha1 = -2, alpha2 = 1)
#' yu2 <- sslRegress(xl, yl, xu, graph.type = "exp", alpha = 1)
#' # Interpolated regularization
#' yu3 <- sslRegress(xl, yl, xu, graph.type = "tanh", alpha1 = -2, alpha2 = 1, method = "Interpolated")
#' yu4 <- sslRegress(xl, yl, xu, graph.type = "exp", alpha = 1, method = "Interpolated")
#' @seealso \code{\link{pr_DB}} \code{\link{dist}}
#' @references Belkin, M., Matveeva, I., & Niyogi, P. (2004a). Regularization and semisupervised learning on large graphs. COLT
sslregress <- function(
    xl, yl, xu, nn = 6, p = 2, method = "Tikhonov", gamma = 1
  ) {
  # all data, both for labelled and unlabelled
  x <- rbind(xl, xu)
  all.Obs <- nrow(x) # number of observations

  # Calculate weights matrix W (equal to adjacency matrix) (nn may be L)
  if (!is.matrix(nn) && is.numeric(nn)) {
    print("Calculating adjacency matrix...")
    L = calculate_laplacian(x, nn)
  } else {
    L = nn
  }
  
  # find closest matrix to L which is SPD, if L is not
  if (!isSymmetric(L)) {
    L <- as.matrix(Matrix::nearPD(L)$mat)
  }
  
  k <- length(yl) # k = number of known obs
  meanY <- mean(yl) # for mean subtraction (section 2.2)
  # join known labels with unknown (labelled with 0), matching xl and xu in x
  y1 <- c(yl, rep(0, all.Obs - k))
  # Calculate power of Laplacian to use as smoothness functional
  s <- rep(1, all.Obs)
  S <- diag(s)
  print(paste0("Calculating S = L^", p, "..."))
  for (i in seq_along(p)) {
    S <- S %*% L
  } 

  rm(L)
  if (method == "Tikhonov") {
    temp <- k * gamma * S + diag(c(rep(1, k), rep(0, all.Obs - k)))
    temp <- solve(temp)
    mu <- as.numeric((-1) * (s %*% (temp %*% y1)) / (s %*% (temp %*% s)))
    f <- temp %*% (y1 + mu)
    yu <- ifelse(f[-(1:k)] > 0, 1, -1)
  }
  if (method == "Interpolated") {
    # partition S as in paper
    S2 <- S[1:k, -(1:k)]
    S3 <- S[-(1:k), -(1:k)]
    t <- matrix(s[-(1:k)], nrow = 1) # ?
    temp <- -solve(S3) %*% t(S2) # S3^{-1} * S_3^T
    y2 <- as.matrix(yl - meanY) # mean subtract
    mu <- as.numeric((-1) * (t %*% (temp %*% y2)) / (t %*% (temp %*% s[1:k])))
    f <- temp %*% (y2 + mu) # f_tilde, as before
    f <- as.vector(f)
    yu <- ifelse(f > 0, 1, -1) # perform binary classification
  }
  return(yu)
}

sslregress_exp <- function(
    xl, yl, xu, nn = 6, t = 0.5, method = "Tikhonov", gamma = 1
  ) {
  # all data, both for labelled and unlabelled
  x <- rbind(xl, xu)
  all.Obs <- nrow(x) # number of observations
  # print(x)

  # Calculate weights matrix W (equal to adjacency matrix) (nn may be L)
  if (!is.matrix(nn) && is.numeric(nn)) {
    print("Calculating adjacency matrix...")
    L = calculate_laplacian(x, nn)
  } else {
    L = nn
  }
  
  # find closest matrix to L which is SPD, if L is not
  if (!isSymmetric(L)) {
    L <- as.matrix(Matrix::nearPD(L)$mat)
  }
  
  k <- length(yl) # k = number of known obs
  meanY <- mean(yl) # for mean subtraction (section 2.2)
  # join known labels with unknown (labelled with 0), matching xl and xu in x
  y1 <- c(yl, rep(0, all.Obs - k))
  # Calculate power of Laplacian to use as smoothness functional
  s <- rep(1, all.Obs)
  S <- expm::expm(-t * L)

  rm(L)
  if (method == "Tikhonov") {
    I <- diag(c(rep(1, k), rep(0, all.Obs - k)))
    temp1 <- (1 / k) * I + gamma * S
    temp1 <- solve(temp1)
    lambda <- as.numeric((2 * s %*% temp1 %*% y1) / (s %*% temp1 %*% s))
    f <- 1 / ((1/k) * gamma * S) %*% ((1/k) * y1 - (lambda/2) * s)
    # perform binary classification
    yu <- ifelse(f[-(1:k)] > 0, 1, -1)
  }
  if (method == "Interpolated") {
    stop("exp(-tL) only implemented for Tikhonov")
  }
  return(yu)
}


#### Calculate classification tables ####

# take training + test data, perform + evaluate semi-supervised predictions
# May not be a pure function!
calc_tables <- \(xl, xu, xu_gen, yl, yu, nn, reg = TRUE, gen_bound = TRUE) {
  misclass_rates <- rep(NA, (length(gammas) + 1))
  if (reg == TRUE) {
    # calculate 6-nearest neighbours graph
    L <- calculate_laplacian(rbind(xl, xu), nn)
    
    # calculate tikhonov solution for various gamma values
    yl_tiks <- lapply(gammas, \(gamma) {
      sslregress(xl = xl, yl = yl, xu = xu, nn = L, p = 2, gamma = gamma)
    })
    
    # interpolated regularization (note condition number smaller than for tik)
    yl_int <- sslregress(xl = xl, yl = yl, xu = xu, nn = L, method = "Interpolated")
    
    # forecasts 
    yls <- c(list(yl_int), yl_tiks) # match order in df
    
    # calculate misclassification rate
    yu_fact <- factor(yu, levels = c(-1, 1))
    yls <- lapply(yls, factor, levels = c(-1, 1))
    # for(x in yls) {
    for(i in seq_along(yls)) {
      # print(norm(yu - x, type = "2"))
      conf_mat <- caret::confusionMatrix(data = yls[[i]], reference = yu_fact)
      misclass_rates[[i]] <- 1 - conf_mat$overall[["Accuracy"]]
    }
  }
  
  gen_bounds <- rep(NA, (length(gammas) + 1))
  if (gen_bound == TRUE) {
    # calculate generalisation error
    for (i in seq_along(gammas)) {
      gen_bound <- calc_gen_bound(xl, xu_gen, 3, gammas[i], 0.01)
      gen_bounds[i+1] <- gen_bound # i = 1 needs to stay NA, for interpolated 
    }   
  }
  
  return(list("misclass_rates" = misclass_rates, "gen_bounds" = gen_bounds))
}

#### Calculate generalisation bound ####

# Qs for Sangeetha:
# - What to do for imaginary lambda? Also, what is "nontrivial" in this context?
# - t = max vertices in training data, right?
# - Dominated by choice of gamma, makes sense right?
# - M always 1 here, right?
# - Should NN be done for training data, or just take nn of training data within
# all data?
# - Generalisation bound should generally get smaller for larger k, not bigger!
# - Isn't generalisation bound NA for interpolated regression, since gamma = 0?

# Calculate generalisation bound, as described in theorem 1
calc_gen_bound <- \(xl, xu, nn, gamma, delta = 0.1) {
  # nn <- 3
  x <- rbind(xl, xu)
  knn <- dbscan::kNN(x, nn, approx = TRUE) 
  w <- knn_to_neighbour(knn) 
  
  # Create an adjacency graph from the adjacency matrix
  graph <- igraph::graph_from_adjacency_matrix(w, mode = "max")
  
  # Calculate Laplacian
  d <- diag(rowSums(w)) # degree matrix
  L <- d - w
  
  # find closest SPD matrix to L
  if (!isSymmetric(L)) {
    L <- as.matrix(Matrix::nearPD(L)$mat)
  }
  
  # Calculate smooth function (L^2)
  S <- diag(rep(1, nrow(x)))
  p <- 2
  for (i in seq_along(p)) {
    S <- S %*% L
  }
  
  k <- nrow(xl) # number of vertices in training data (higher => lower bound)
  M = 1 # bound on max of vertex values (always 1?)
  t <- max(degree(graph))
  lambda1 <- eigen(S)$values
  lambda1 <- min(lambda1[lambda1 != 0]) 
  print(paste("lambda1 =", lambda1))
  
  # diameter of graph (maximum shortest path length between any pair of vertices)
  D <- diameter(graph)
  # Maximum entry of gamma-regularisation problem (prop 2)
  K = M * sqrt(D/gamma)
  
  beta = (3 * M * sqrt(t * k)) / ((k * gamma * lambda1 - t)^2) + 
    (4 * M) / (k * gamma * lambda1 - t)
  
  gen_bound <- beta + sqrt((2 * log(2 / delta)) / k) * (k * beta + (K + M)^2)
  
  return(gen_bound)
}


#### Calculate Laplacian Matrix ####

# calculate Laplacian L = D - W, D = Degree matrix, W = adjacency matrix
calculate_laplacian <- function(x, nn) { # p) {
# calculate_laplacian <- function(xl, xu, nn) { # p) {
  # x <- rbind(xl, xu)
  # all.Obs <- nrow(x)
  knn <- dbscan::kNN(x, nn, approx = TRUE) # "exact" wayyyy too slow here
  # knn <- train_test_neighbours(xl, xu, nn)
  w <- knn_to_neighbour(knn) 
  
  # compute Laplacian as in section 2.1 of paper
  d <- diag(rowSums(w)) # degree matrix
  L <- d - w
  return(L)
}

# Find closest train neighbours in test
train_test_neighbours <- function(xl, xu, nn = 6) {
  
  x <- rbind(xl, xu)
  nl <- nrow(xl)
  nu <- nrow(xu)
  # if looking at training error, just return network as normal
  if (nl == nu && all(xl == xu)) {
    return( dbscan::kNN(x, nn, approx = TRUE))
  }
  train_set <- seq_len(nl)
  test_set <- seq_len(nu) + nl
  max_nn <- nl + nu - 1
  
  # take large number of nearest neighbours
  nn_large <- min(100, max_nn)
  knn <- dbscan::kNN(x, nn_large, approx = TRUE)
  
  # needs to be nn neighbours between each test and training point, if not do
  # above again
  test_fun <- function(knn) {
    apply(knn$id[test_set, ], 1, \(y) {
      length(y[y %in% train_set])
    })
  }
  # increase neighbourhood, try again!
  while (any(test_fun(knn) < nn)) {
    print("Not enough neighbours, trying again!")
    nn_large <- min(nn_large + 25,max_nn - 1)
    if (nn_large == nl + nu) stop("Can't increase any further!")
    knn <- dbscan::kNN(x, nn_large, approx = TRUE)
  }
  
  knn_final <- knn
  # take first nn neighbours for training set 
  knn_final$id <- knn$id[train_set, 1:nn]
  # only take neighbours in test set that are in training set
  test_train_nn <- t(apply(knn$id[test_set, ], 1, \(y) {
    train <- y[y %in% train_set]
    if (length(train) == 0) stop("No nn")
    # if point doesn't have nn neighbours in train, join to other test values
    if (length(train) < nn) {
      print("train not long enough, adding test points in")
      n_remain <- nn - length(train)
      train <- c(y[y %in% test_set][1:n_remain], train)
      train <- train[order(as.numeric(names(train)))]
    } else {
      train <- train[1:nn]
    }
  }))
  knn_final$id <- rbind(knn_final$id, test_train_nn)
  return(knn_final)
}


# calculate adjacency matrix from knn object
knn_to_neighbour <- function(knn) {
  n <- nrow(knn$id)
  N <- diag(0, n)
  for (i in seq_len(n)) {
    N[i, knn$id[i, , drop = TRUE]] <- 1
  }
  return(N)
}

