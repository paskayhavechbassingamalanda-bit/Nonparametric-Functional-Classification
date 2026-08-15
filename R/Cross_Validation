
vgc1 <- function(n1, n2, a, sigma) {
  
  don <- as.matrix(data2(n1, n2, a, sigma))
  nxn <- n1 * n2
  
  argvals <- seq(0, 1, 1/364)
  p <- length(argvals)
  
  X <- as.matrix(don[1:nxn, 1:p])
  Y <- don[1:nxn, p+1]   # Response variable
  delta <- don[1:nxn, p+2]
  Y_G <- don[1:nxn, p+3] # True response variable
  
  n_folds <- 10
  rn <- seq(0.5, 2.5, length.out = 10)
  Kmax <- min(50, floor(nxn / 3))
  K <- seq(5, Kmax, by = 3)
  
  
  ## Auxiliary function: distance
  dista = function(x, y) {
    
    x = as.matrix(x)
    y = as.matrix(y)
    
    m1 = nrow(x)
    m2 = nrow(y)
    p = ncol(x)
    
    D = matrix(0, m1, m2)
    
    for (i in 1:m1) {
      for (j in 1:m2) {
        D[i, j] = sqrt(
          t(x[i, ] - y[j, ]) %*%
            (x[i, ] - y[j, ]) / p
        )
      }
    }
    
    D
  }
  
  
  ## Closed ball
  ind2 = function(A, B, rn) {
    
    A = as.matrix(A)
    B = as.matrix(B)
    rn = as.numeric(rn)
    
    m1 = nrow(A)
    m2 = nrow(B)
    
    G = array(0, c(m1, m2, length(rn)))
    
    for (i in 1:m1) {
      for (j in 1:m2) {
        for (k in 1:length(rn)) {
          
          if (dista(A, B)[i, j] > rn[k]) {
            G[i, j, k] = 0
          } else {
            G[i, j, k] = 1
          }
        }
      }
    }
    
    G
  }
  
  
  ## Binary transformation
  ind1 <- function(A) {
    ifelse(abs(as.numeric(A)) > 3/2, 1, 0)
  }
  
  
  ## Epanechnikov kernel
  K2 <- function(A) {
    ifelse(
      abs(A) <= 1,
      (3/4) * (1 - A^2),
      0
    )
  }
  
  
  ## U_sub function
  U_sub <- function(Y_G, X, delta, rn, h) {
    
    v <- which(delta == 1)      # Uncensored observations
    vtilde <- which(delta == 0) # Censored observations
    
    nv <- length(v)
    nc <- length(vtilde)
    
    U <- numeric(nc)
    M <- matrix(0, nc, nv)
    
    for (j in 1:nc) {
      
      for (l in 1:nv) {
        
        in_ball <- ind2(
          matrix(X[v[l], ], 1),
          matrix(X[vtilde[j], ], 1),
          rn
        )
        
        M[j, l] <- in_ball *
          K2((Y_G[v[l]] - Y_G[vtilde[j]]) / h)
        
        U[j] <- U[j] +
          ind1(Y_G[v[l]]) * M[j, l]
      }
      
      sum_M <- sum(M[j, ])
      
      if (sum_M > 0) {
        U[j] <- U[j] / sum_M
      } else {
        U[j] <- 0
      }
    }
    
    return(U)
  }
  
  
  ## Cross-validation
  set.seed(123)
  
  data <- cbind(X, Y, Y_G, delta)
  shuffled_data <- data[sample(nrow(data)), ]
  
  
  # Extraction of the shuffled variables
  X <- as.matrix(shuffled_data[, 1:p])
  Y <- shuffled_data[, p+1]
  Y_G <- shuffled_data[, p+2]
  delta <- shuffled_data[, p+3]
  
  
  # Splitting the data into folds
  folds <- cut(
    seq(1, nrow(shuffled_data)),
    breaks = n_folds,
    labels = FALSE
  )
  
  
  # Error array
  error <- array(
    0,
    dim = c(length(rn), length(K), n_folds)
  )
  
  
  for (r_idx in 1:length(rn)) {
    
    rn_val <- rn[r_idx]
    
    for (k_idx in 1:length(K)) {
      
      K_val <- K[k_idx]
      
      for (l in 1:n_folds) {
        
        test_indices <- which(
          folds == l,
          arr.ind = TRUE
        )
        
        test_data <- shuffled_data[test_indices, ]
        train_data <- shuffled_data[-test_indices, ]
        
        
        # Training data
        X_train <- train_data[, 1:p]
        Y_train <- train_data[, p+1]
        Y_G_train <- train_data[, p+2]
        delta_train <- train_data[, p+3]
        
        
        # Test data
        X_test <- test_data[, 1:p]
        Y_test <- test_data[, p+1]
        Y_G_test <- test_data[, p+2]
        delta_test <- test_data[, p+3]
        
        n_train <- nrow(X_train)
        n_test <- nrow(X_test)
        
        
        # Estimation of U
        Un_hat <- U_sub(
          Y_G_train,
          X_train,
          delta_train,
          rn = rn_val,
          h = n_train^(-0.2)
        )
        
        
        # Prediction
        pred <- numeric(n_test)
        
        idx_V <- which(delta_train == 1)
        idx_Vtilde <- which(delta_train == 0)
        
        
        for (i in 1:n_test) {
          
          xi <- matrix(X_test[i, ], 1)
          
          rv <- 0
          rc <- 0
          
          
          if (length(idx_V) > 0) {
            
            for (a in idx_V) {
              
              rv <- rv +
                ind1(Y_G_train[a]) *
                ind2(
                  matrix(X_train[a, ], 1),
                  xi,
                  rn_val
                )
            }
          }
          
          
          if (length(idx_Vtilde) > 0) {
            
            for (j in seq_along(idx_Vtilde)) {
              
              rc <- rc +
                Un_hat[j] *
                ind2(
                  matrix(X_train[idx_Vtilde[j], ], 1),
                  xi,
                  rn_val
                )
            }
          }
          
          
          term1 <- rv / length(idx_V) *
            sqrt(K_val / n_train)
          
          term2 <- rc / length(idx_Vtilde) *
            sqrt(K_val / n_train)
          
          pred[i] <- term1 + term2
        }
        
        
        # Store the error in the appropriate position
        error[r_idx, k_idx, l] <-
          mean(
            (pred - ind1(Y_G_test))^2
          )
      }
    }
  }
  
  
  # Mean cross-validation error
  mean_error <- apply(
    error,
    c(1, 2),
    mean
  )
  
  
  # Best parameter pair
  best_idx <- which(
    mean_error == min(mean_error),
    arr.ind = TRUE
  )
  
  best_rn <- rn[best_idx[1]]
  best_K <- K[best_idx[2]]
  best_error <- mean_error[best_idx]
  
  
  return(
    list(
      cv = best_error,
      best_rn = best_rn,
      best_K = best_K
    )
  )
}
