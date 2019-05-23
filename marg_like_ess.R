chain_comp = function(ess){
  
  # Given an EG with adjacency matrix "ess" returns its set of chain components Tau
  
  ess = as(ess,"igraph")
  
  amat = as.matrix(get.adjacency(ess))  # if the argument is ess.obs (the essential graph, a graph object)
  wmat = matrix(as.integer((amat + t(amat)) > 1), nrow = nrow(amat))
  wg = graph.adjacency(wmat, mode = "undirected")
  cc = clusters(wg)
  neworder  =  order(cc$membership)
  a = matrix(0, nrow = length(cc$csize), ncol = length(cc$csize))
  b = cumsum(cc$csize)
  wmat = amat[neworder, neworder]
  
  for(i in 1: length(cc$csize)){
    for(j in 1: length(cc$csize)){
      if(j != i){
        a[i,j] = as.integer(sum(wmat[(max(b[i-1],0)+1):b[i],
                                      (max(b[j-1],0)+1):b[j]]) > 0)
      }
    }
  }
  
  rownames(a) = colnames(a) = as.character(1:length(b))
  
  chainorder = topOrder(a)
  vertorder = c()
  chainsize = c()
  
  for(k in 1:length(b)){
    vertorder = c(vertorder, which(cc$membership == chainorder[k]))
    chainsize = c(chainsize, cc$csize[chainorder[k]])
  }
  
  q = list(vert.order=vertorder,chain.size=chainsize)
  order = q$vert.order
  size = c(1,q$chain.size)
  
  Tau = list(0)
  for(i in 2:(length(size))){
    Tau[[i-1]] = order[sum(size[1:(i-1)]):(sum(size[2:i]))]
  }
  
  return(Tau)
  
}


marg_like = function(ess,Y){
  
  # Returns the FBF (Fractional Bayes Factor) marginal likelihood of an EG as in Castelletti et al. (2018)  
  
  # ess : (q,q) adjacency matrix of an EG
  # Y   : (n,q) data matrix
  
  n = nrow(Y)
  A = ess
  ess = as(ess,"graphNEL")
  
  # need to specify first ess.obs (the essential graph)
  
  Tau = chain_comp(ess)
  return(sum(unlist(lapply(Tau,m_tau,ess = ess, Y = Y))))
  
}

Gamma = function(s,a){
  (.25*s*(s-1))*log(pi)+sum(lgamma(a+.5*(1-1:s)))
}

Ehat = function(y,x){
  beta = solve(t(x)%*%x)%*%t(x)%*%y
  y-x%*%beta
}

m_J = function(J,pa_tau,tau_c,pa_tau_c,a_D,n_0,Y){
  
  n = nrow(Y)
  
  # J is, from time to time, a Clique or a Separator
  
  if(identical(unlist(J), character(0))){
    m = 0
  }
  
  else{
    
    J = as.numeric(J)
    Y_J = as.matrix(Y[,J])
    J_c = length(J)
    
    X_J = as.matrix(cbind(1,Y[,pa_tau]))
    
    m = (-.5*(n-n_0)*J_c)*log(pi)+
      Gamma(J_c,0.5*(a_D+n-pa_tau_c-1-tau_c+J_c))-Gamma(J_c,0.5*(a_D+n_0-pa_tau_c-1-tau_c+J_c))+
      .5*(J_c*(a_D+n_0-tau_c+J_c))*log(n_0/n)+
      (-.5*(n-n_0))*log(det(t(Ehat(Y_J,X_J))%*%Ehat(Y_J,X_J)))
    
  }
  
  m
  
}

m_tau = function(tau,ess,Y){
  
  # Compute the marginal likelihood for a chain component tau
  
  tau_c = length(tau)
  pa_tau = as.numeric(parents(tau,ess))
  pa_tau_c = length(pa_tau)
  n_0 = ncol(Y)+1
  a_D = tau_c-1
  
  if(tau_c <= 2){
    
    m_J(tau,pa_tau,tau_c,pa_tau_c,a_D,n_0,Y)
    
  } else{
    
    ess = as(ess,"graphNEL")
    taug = subGraph(as.character(tau),ess)
    
    C = mpd(taug)$cliques
    S = mpd(taug)$separators
    
    sum(unlist(lapply(C,m_J,pa_tau,tau_c,pa_tau_c,a_D,n_0,Y)))-
      sum(unlist(lapply(S,m_J,pa_tau,tau_c,pa_tau_c,a_D,n_0,Y)))
    
  }
} 