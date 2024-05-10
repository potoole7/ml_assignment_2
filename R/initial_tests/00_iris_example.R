#### Semi-supervised learning using Tikhonov & Interpolated Regularisation ####
# Tested on iris data

library(igraph)

source("R/functions.R")

data(iris)
# labelled data
xl<-iris[c(1:20,51:70),-5] # just take setosa and versicolor, remove label
# unlabelled data
xu<-iris[c(21:50,71:100),-5]
# labels for labelled data
yl<-rep(c(1,-1),each=20) 
# true labels
yu <- rep(c(1,-1),each=30) 


#### Training Error ####

nn_train <- 2
Ll <- calculate_laplacian(rbind(xl, xl), nn_train)
# Ll <- calculate_laplacian(xl, xl, 6)

# Tikhonov regularization (on training data)
yl1 <- sslregress(xl, yl, xl, nn = Ll, p = 2, gamma = 1)
# try for exponential smooth function
yl1_exp <- sslregress_exp(xl, yl, xl, nn = Ll, t = 1, gamma = 1)

# Interpolated regularization
yl2 <- sslregress(xl, yl, xl, nn = Ll, p = 2, method = "Interpolated")

# compare true labels to regularisation solution (i.e. compute test error)
yls <- list(yl1, yl1_exp, yl2)
for(x in yls) {
  print(norm(yl - x, type = "2")) # /length(yl))
  # compute classification error rate on training data
  conf_mat <- caret::confusionMatrix(factor(x), factor(yl))
  print(paste0("Misclassification Rate = ", 1 - conf_mat$overall[["Accuracy"]]))
}


#### Test Error ####

nn_test <- 3

Lu <- calculate_laplacian(rbind(xl, xu), nn_test)

# Tikhonov regularization (on test)
yu1 <- sslregress(xl, yl, xu, nn = Lu, p = 2, gamma = 1)
yu1_exp <- sslregress_exp(xl, yl, xu, nn = Lu, t = 1, gamma = 1)


# Interpolated regularization
yu2 <- sslregress(xl, yl, xu, nn = Lu, p = 2, method = "Interpolated")

# compare true labels to regularisation solution (i.e. compute test error)
yus <- list(yu1, yu1_exp, yu2)
for(x in yus) {
  print(norm(yu - x, type = "2"))
  conf_mat <- caret::confusionMatrix(factor(x), factor(yu))
  print(paste0("Misclassification Rate = ", 1 - conf_mat$overall[["Accuracy"]]))
}

# how to compute training error here? Need to return this as well!

#### Calculate Generalisation Error Bound ####

lapply(c(1, 0.1, 0.01, 0.001), \(x) {
  print(paste(
    "Generalisation bound for gamma =", 
    x, 
    "is",
    round(calc_gen_bound(xl, xu, 3, x, 0.1), 3)
  ))
})

# calculate nearest-neighbours adjacency matrix between training and test data
nn <- 3
x <- rbind(xl, xu)
# TODO: Allow this to be input to calculate_laplacian
knn <- dbscan::kNN(x, nn, approx = TRUE) # "exact" wayyyy too slow here
# knn <- train_test_neighbours(xl, xu, nn)
w <- knn_to_neighbour(knn) 

# Create an adjacency graph from the adjacency matrix
graph <- graph.adjacency(w, mode = "undirected")
graph_train <- graph.adjacency(
  w[seq_len(nrow(xl)), seq_len(nrow(xl))], 
  mode = "undirected"
)
plot(graph) # notice that it is unconnected!

# Calculate Laplacian
d <- diag(rowSums(w)) # degree matrix
L <- d - w

# find closest SPD matrix to L
if (!isSymmetric(L)) {
  L <- as.matrix(Matrix::nearPD(L)$mat)
}

# Calculate smooth function (L^2)
S <- L %*% L

delta <- 0.1 # 1 - delta = confidence level
gamma <- 0.1 # regularisation parameter
k <- nrow(xu) # number of vertices in training data
M = 1 # bound on max of vertex values (always 1?)
t <- max(degree(graph_train))
lambda1 <- Re(eigen(S)$values)  # also take real component, in case of Im (?)
lambda1 <- min(lambda1[lambda1 != 0]) 

# diameter of graph (maximum shortest path length between any pair of vertices)
D <- diameter(graph)
# Maximum entry of gamma-regularisation problem (prop 2)
K = M * sqrt(D/gamma)

beta = (3 * M * sqrt(t * k)) / ((k * gamma * lambda1 - t)^2) + 
  (4 * M) / (k * gamma * lambda1 - t)

gen_bound <- beta + sqrt((2 * log(2 / delta)) / k) * (k * beta + (K + M)^2)
print(paste("gen_bound =", round(gen_bound, 3)))
