##### Monte Carlo simulation

MC = function(n1, n2, a, rn, K, sigma) {
  
  M = 300  # Number of replications
  
  ## Storage of the results: nxn = n1 * n2
  results <- matrix(
    0,
    nrow = M,
    ncol = nxn
  )
  
  ## Replications
  for (m in 1:M) {
    
    results[m, ] <- erreur_prd(
      n1,
      n2,
      a,
      rn,
      K,
      sigma
    )$ep
  }
  
  return(results)
}


i = n1 * n2 * M


##### MSE

MSE = sum(
  (
    erreur_prd(
      n1,
      n2,
      a,
      rn,
      K,
      sigma
    )$et -
      t(MC(
        n1,
        n2,
        a,
        rn,
        K,
        sigma
      ))
  )^2
) / i
