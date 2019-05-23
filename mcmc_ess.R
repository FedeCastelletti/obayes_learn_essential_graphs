# MCMC scheme

library(pcalg)
library(igraph)
library(gRbase)
library(graph)
library(ggm)

source("markov_chain_move.r")
source("marg_like_ess.r")


post_ess_graphs = function(Y, m, T, verbose = FALSE){
  
  # Y       : the (n,q) data matrix
  # m       : maximum number of edges (undirected or directed) in the EG space
  # T       : number of MCMC iterations
  # verbose : TRUE/FALSE to see or not plot of posterior density of visited EGs
  
  # Output:
  
  # A_chain : a (q*q,T) matrix collecting the (vectorized) adjacency matrices of the EGs visited by the chain
  
  q = ncol(Y)
  
  # initial graph (empty graph)
  
  A_0 = matrix(0,q,q)
  colnames(A_0) = rownames(A_0) = 1:q
  
  # matrix collecting the output of the MCMC (each column is the by-column-vectorized adjacency matrix of the accepted graph)
  A_chain = matrix(0, nrow = q*q, ncol = T+1)
  A_chain[,1] = A_0
  
  # store space for logposterior
  logpost = rep(-Inf, T+1)
  
  A_old = A_0
  m_old = marg_like(A_old, Y)
  G_und = G_undirected(A_old)
  
  colnames(A_0) = rownames(A_0) = 1:q
  
  # generate randoms for MCMC
  # runifs = log(runif(T,0,0.3))

  # MCMC iterations
  for(t in 2:T){
    
    # proposed move
    prop = move(A_old,m,G_und=G_und,q)
    A_new = prop[[1]]
    m_new = marg_like(A_new,Y)
    
    # log prior ratio evaluation
    logprior.new = lgamma(n.edge(A_new)+1) + 
      lgamma(q*(q-1)/2-n.edge(A_new)+(2*q-2)/3-1)
    logprior.old = lgamma(n.edge(A_old)+1) + 
      lgamma(q*(q-1)/2-n.edge(A_old)+(2*q-2)/3-1)
    logprior = logprior.new - logprior.old
       
    # acceptance ratio
    ratio = min(0,m_new-m_old+logprior)
    # ratio = min(0,m_new-m_old+logprior)
    
    if(log(runif(1,0,1)) < ratio){ # accept move
      
      A_old = A_new
      m_old = m_new
      G_und = G_undirected(A_old)
      logprior.old = logprior.new
      }
    
    # store chain value and posterior density
    A_chain[,t] = A_old
    logpost[t]  = m_old + logprior.old
    # print(ratio)
    
    # show plot of posterior if verbose
    if(verbose){
      if(t%%10==0 & t > 150) plot(logpost[100:t],type="l")
    }
 
  }
  
  return(A_chain)
  
}
