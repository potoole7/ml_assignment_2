#### Semi-supervised learning using Tikhonov & Interpolated Regularisation ####
# Tested on iris data


# TODO: Investigate why S3 singular for small nn for test but not train
# TODO: Plot adjacency matrices for S3s, see if training data not connected to test
# Maybe Sangeetha/Bill might know? I've tried everything!
# TODO: Calculate classification error rate (test and train) (done)
# Presume this just means 1 - accuracy?
# TODO: Calculate generalisation bound

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
# yl1<-sslRegress(xl,yl,xl,graph.type="tanh",alpha1=-2,alpha2=1, nn = 2)
# yl2<-sslRegress(xl,yl,xl,graph.type="exp",alpha = 1, nn = 2)
yl1 <- sslregress(xl, yl, xl, nn = Ll, p = 2, gamma = 1)
# try for exponential smooth function
yl1_exp <- sslregress_exp(xl, yl, xl, nn = Ll, t = 1, gamma = 1)

# Interpolated regularization
yl2 <- sslregress(xl, yl, xl, nn = Ll, p = 2, method = "Interpolated")

# compare true labels to regularisation solution (i.e. compute test error)
# yls <- list(yl1, yl2, yl3, yl4)
# yls <- list(yl1, yl2)
yls <- list(yl1, yl1_exp, yl2)
for(x in yls) {
  print(norm(yl - x, type = "2")) # /length(yl))
  # compute classification error rate on training data
  conf_mat <- caret::confusionMatrix(factor(x), factor(yl))
  print(paste0("Misclassification Rate = ", 1 - conf_mat$overall[["Accuracy"]]))
}


#### Test Error ####

# For nn_test = 2, algorithm fails, but for nn_test = 3, works fine
# For both, graph is unconnected, so this is not a condition for fitting
# However, for both algorithms matrix to be inverted is singular for nn_test = 2
# Why???
# For nn_test = 2, det(Lu) = 0, while for nn_test = 3, det(Lu) > 0 (but v small)
# Therefore, Lu cannot be singular, or else the algorithm fails!
# Is there anything on how to have a non-singular Laplacian?
# However, for train data, can have up to 5 nn and still have det(Ll) = 0 
# but inverted matrix will still be non-singular! Strange
# So for training, constant * S + I is non-singular for singular S, but for 
# test it is? interesting!
# TODO: Need to find some papers with theoretical guarantees on this kind of 
# thing
# Rowsums are 0 for both Lu, Ll, colsums both != 0
# Both not symmetric
# What's the difference?! Look at interpolated and S3 more closely:
# Again, det(S3_test) = 0, not the case for det(S3_train)! But why?
# Recall, this is with det(Ll) = det(Lu) = 0
# L is diagonally dominant (the magnitude of the diagonal entry is greater than or equal to the sum of the magnitudes of the off-diagonal entries in its row)
# Where equal, L is singular, but this doesn't seem to matter here
# But how is S3_test different to S3_train?!?! Other than having 0 determinant
# last two eigenvalues of S3_train are 1, 1, while for S3_test they're
# 0.3819660112501051529854 -0.0000000000000001982447, so ill-conditioned! 
# (Similar before that)

nn_test <- 3

# not working, fix!
# Still causes error!
Lu <- calculate_laplacian(rbind(xl, xu), nn_test)
# Lu <- calculate_laplacian(xl, xu, 6)

# Tikhonov regularization (on test)
# yu1<-sslRegress(xl,yl,xu,graph.type="tanh",alpha1=-2,alpha2=1)
# yu2<-sslRegress(xl,yl,xu,graph.type="exp",alpha = 1)
yu1 <- sslregress(xl, yl, xu, nn = Lu, p = 2, gamma = 1)
yu1_exp <- sslregress_exp(xl, yl, xu, nn = Lu, t = 1, gamma = 1)


# Interpolated regularization
# yu3<-sslRegress(xl,yl,xu,graph.type="tanh",alpha1=-2,alpha2=1,method="Interpolated")
# yu4<-sslRegress(xl,yl,xu,graph.type="exp",alpha = 1,method="Interpolated")
yu2 <- sslregress(xl, yl, xu, nn = Lu, p = 2, method = "Interpolated")

# compare true labels to regularisation solution (i.e. compute test error)
# yus <- list(yu1, yu2, yu3, yu4)
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

# Calculate smooth function (L^2)
S <- L %*% L

delta <- 0.1 # 1 - delta = confidence level
gamma <- 0.1 # regularisation parameter
k <- nrow(xu) # number of vertices in training data
# k <- 1000 # higher k => lower generalisation bound
M = 1 # bound on max of vertex values (always 1?)
# max multiplicity of vertices (always nn?) (colSums?) 
# May also be in just training set, not all data!
# t <- max(degree(graph)) # lower for this
t <- max(degree(graph_train))
# lambda1 # (second?) smallest non-trivial eigenvalue of S (may be imaginary!)
# Does non-trivial mean != 0, or greater than a very small value?
lambda1 <- Re(eigen(S)$values)  # also take real component, in case of Im (?)
# lambda1 <- min(lambda1[lambda1 > 1e-5]) # doesn't have much effect
lambda1 <- min(lambda1[lambda1 != 0]) 
# lambda1 <- 2 # for larger lambda, bound increases (Should this be the case?)

# diameter of graph (maximum shortest path length between any pair of vertices)
D <- diameter(graph)
# Maximum entry of gamma-regularisation problem (prop 2)
K = M * sqrt(D/gamma)

beta = (3 * M * sqrt(t * k)) / ((k * gamma * lambda1 - t)^2) + 
  (4 * M) / (k * gamma * lambda1 - t)

gen_bound <- beta + sqrt((2 * log(2 / delta)) / k) * (k * beta + (K + M)^2)
print(paste("gen_bound =", round(gen_bound, 3)))
# definitely vacuous, probably too large at ~3000, ~2000 largest seen for 
# ionosphere in paper
# A lot smaller for larger gamma! but classification error worse for higher gamma?!
