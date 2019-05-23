# Example

data = read.csv("seeds.txt", sep = "", header = FALSE)

head(data)

source("mcmc_ess.r")

T = 1000

# Run the MCMC

out = post_ess_graphs(Y = data, m = 14, T = T, verbose = TRUE)

# Compute (estimated) posterior probabilities of edge inclusion

q = ncol(data)

burn = 200

probs  = matrix((rowMeans(out[,(burn + 1):(T)])), q, q)
probs
