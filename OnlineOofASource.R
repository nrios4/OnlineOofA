# OnlineOofASource.R

rm(list = ls())

library(OofAExp)
library(gtools)


# helper function
# takes a matrix D where each row is a permutation of 1:m
# returns the CP model matrix with baseline coding, excluding the intercept
get_CP_model_matrix <- function(D) {
  
  P <- gen.cp(D)
  
  P <- as.matrix(P)
  n <- nrow(P)
  m <- ncol(P)
  
  # Position of each component in each permutation
  pos <- matrix(NA, n, m)
  
  for (i in seq_len(n)) {
    pos[i, ] <- match(seq_len(m), P[i, ])
  }
  
  # Convert positions to factors with SAS contrasts
  facs <- vector("list", m)
  
  for (j in seq_len(m)) {
    facs[[j]] <- factor(pos[, j], levels = seq_len(m))
    contrasts(facs[[j]]) <- contr.SAS(m)
  }
  
  dat <- as.data.frame(facs)
  names(dat) <- paste0("cp", seq_len(m))
  
  model.matrix(~ ., data = dat)[,-1]
}

# fit_all_models
# This function fits all of the OofA models in "models" to the data
# Inputs:
# ---- models: a vector of model names, options are: "pwo" (pairwise order),
#     "cp" (component-position), "lcpf" (first order flexible position),
#      or "pwo1", "pwo2", "pwo3" (but these are only used for Section 4 example).
# ---- D: a (n by m) matrix where each row is a permutation of 1:m
# ---- y: a vector of length n containing the responses
fit_all_models <- function(models, D, y){
  
  M <- length(models)
  model_list <- vector("list", M)
  # m = ncol(D)
  
  for(j in 1:M){
    
    if(models[j] == "cp"){
      X <- get_CP_model_matrix(D)
    }
    
    if(models[j] == "lcpf"){
      X <- gen.pf(D, "lcpf")
      X <- X[, -ncol(D)]
    }
    
    if(models[j] == "pwo"){
      X <- gen.pf(D, "pwof")
    }
    
    if(models[j] == "pwo1"){
      X <- gen.pf(D, "pwof")
      X <- X[, seq(1, ncol(X), by = 3), drop = FALSE]
    }
    
    if(models[j] == "pwo2"){
      X <- gen.pf(D, "pwof")
      X <- X[, seq(2, ncol(X), by = 3), drop = FALSE]
    }
    
    if(models[j] == "pwo3"){
      X <- gen.pf(D, "pwof")
      X <- X[, seq(3, ncol(X), by = 3), drop = FALSE]
    }
    
    
    if(models[j] == "lcpf+qcpf"){
      X1 = gen.pf(D,"lcpf")
      X2 = gen.pf(D,"qcpf")
      X = cbind(X1,X2)
    }
    
    df_tmp <- data.frame(cbind(X, y))
    
    fit_tmp = lm(y ~ ., data = df_tmp)
    
    model_list[[j]] <- list(
      type = models[j],
      fit  = fit_tmp,
      loglik = as.numeric(logLik(fit_tmp))
    )
  }
  
  model_list
}

predict_all_models <- function(model_list, new_order){
  
  Dnew <- matrix(new_order, nrow = 1)
  m = ncol(D)
  
  preds <- numeric(length(model_list))
  
  for(j in seq_along(model_list)){
    
    model_type <- model_list[[j]]$type
    fit <- model_list[[j]]$fit
    
    if(model_type == "cp"){
      Xnew <- t(get_CP_model_matrix(Dnew))
    }
    
    if(model_type == "lcpf"){
      Xnew <- gen.pf(Dnew, "lcpf")
      Xnew <- Xnew[, -ncol(Dnew), drop = FALSE]
    }
    
    if(model_type == "pwo"){
      Xnew <- gen.pf(Dnew, "pwof")
    }
    
    if(model_type == "pwo1"){
      Xnew <- gen.pf(Dnew, "pwof")
      Xnew <- Xnew[, seq(1, ncol(Xnew), by = 3), drop = FALSE]
    }
    
    if(model_type == "pwo2"){
      Xnew <- gen.pf(Dnew, "pwof")
      Xnew <- Xnew[, seq(2, ncol(Xnew), by = 3), drop = FALSE]
    }
    
    if(model_type == "pwo3"){
      Xnew <- gen.pf(Dnew, "pwof")
      Xnew <- Xnew[, seq(3, ncol(Xnew), by = 3), drop = FALSE]
    }
    
    if(model_type == "lcpf+qcpf"){
      X1 = gen.pf(Dnew,"lcpf")
      X2 = gen.pf(Dnew,"qcpf")
      Xnew = cbind(X1, X2)
    }
    
    preds[j] <- predict(fit,newdata = data.frame(Xnew))
  }
  
  names(preds) <- sapply(model_list, `[[`, "type")
  
  preds
}


# generate_response
# Helper function to online_oofa_sampling
# generates the true response under a specified true model
# Inputs:
# ---- sampled_order: a permutation of 1:m
# ---- true_model: one of "pwo", "cp", "lcpf", or "lcpf+qcpf"
# ---- sigma: error standard deviation
generate_response <- function(sampled_order, true_model = "lcpf+qcpf",
                              sigma = 0.05){
  
  m = length(sampled_order)
  if(true_model == "pwo"){
    x = gen.pf(sampled_order, "pwof")
    p = choose(m,2)
    beta_vec = 1:p/p
    y = as.numeric(x%*%beta_vec  + rnorm(1,mean = 0, sd = sigma))
    return(y)
  }
  if(true_model == "cp"){
    x = get_CP_model_matrix(sampled_order)
    p = length(x)
    # beta_vec = 1:p/p
    opt = m:1 # we want this to be the optimal permutation
    pvec <- match(1:m, opt)
    B <- outer(1:m,1:(m-1),function(i,j) (j - pvec[i])^2
    )
    
    beta_vec <- as.vector(t(B))
    y = as.numeric(x%*%beta_vec  + rnorm(1,mean = 0, sd = sigma)) + 1 # add 1 to avoid dividing by 0 in relative optimality
    return(y)
  }
  if(true_model == "lcpf"){
    x = gen.pf(sampled_order,"lcpf")
    beta1 = 1:m/m
    y = as.numeric(x%*%beta1  + rnorm(1,mean = 0, sd = sigma))
    return(y)
  }
  if(true_model == "lcpf+qcpf"){
    x1 = gen.pf(sampled_order,"lcpf")[-m]
    x2 = gen.pf(sampled_order,"qcpf")[-c(m-1,m)]
    beta1 = 1:(m-1)/m
    beta2 = seq(from = 0, to = 1, length.out = m-2)
    y = as.numeric(x1%*%beta1 + x2%*%beta2 + rnorm(1,mean = 0, sd = sigma))
    return(y)
  }
  
  
}

# get_neighboring_orders
# Helper function that finds the set of orders that "neighbor" the current best order
# Inputs:
# ---- best_order: a permutation of 1:m
# Returns: a matrix where row is a permutation of 1:m that neighbors best_order
get_neighboring_orders <- function(best_order) {
  
  m <- length(best_order)
  
  neighbors <- matrix(NA, nrow = m - 1, ncol = m)
  
  for (i in 1:(m-1)){
    x <- best_order
    x[c(i, i + 1)] <- x[c(i + 1, i)]
    neighbors[i, ] <- x
  }
  
  neighbors
}

# currently, this function just finds the optimal order for given true_models
# and it just uses brute force search
# here, optimal means MINIMUM
get_true_optimal_order <- function(true_model, m){
  
  all_permutations = permutations(n=m,r=m)
  true_ys = apply(all_permutations, 1, generate_response,
                  true_model = true_model, sigma = 0) # sigma = 0 gives no error
  min_ind = which.min(true_ys)
  return(list(opt_order = all_permutations[min_ind,],
              opt_y = true_ys[min_ind]))
  
}

get_design_row <- function(order, model_type){
  
  Dnew <- matrix(order, nrow = 1)
  
  if(model_type == "cp"){
    X <- get_CP_model_matrix(Dnew)
  }
  
  if(model_type == "lcpf"){
    X <- gen.pf(Dnew, "lcpf")
    X <- X[, -ncol(Dnew), drop = FALSE]
  }
  
  if(model_type == "pwo"){
    X <- gen.pf(Dnew, "pwof")
  }
  
  if(model_type == "pwo1"){
    X <- gen.pf(Dnew, "pwof")
    X <- X[, seq(1, ncol(X), by = 3), drop = FALSE]
  }
  
  if(model_type == "pwo2"){
    X <- gen.pf(Dnew, "pwof")
    X <- X[, seq(2, ncol(X), by = 3), drop = FALSE]
  }
  
  if(model_type == "pwo3"){
    X <- gen.pf(Dnew, "pwof")
    X <- X[, seq(3, ncol(X), by = 3), drop = FALSE]
  }
  
  
  if(model_type == "lcpf+qcpf"){
    X1 <- gen.pf(Dnew, "lcpf")
    X2 <- gen.pf(Dnew, "qcpf")
    X  <- cbind(X1, X2)
  }
  
  X
}

model_ei <- function(fit, model_type,
                     neighbor_order,
                     incumbent_order){
  
  Xa <- get_design_row(neighbor_order, model_type)
  Xb <- get_design_row(incumbent_order, model_type)
  
  beta_hat <- coef(fit)
  Vbeta    <- vcov(fit)
  
  # Keep estimable coefficients only
  keep <- !is.na(beta_hat)
  
  beta_hat <- beta_hat[keep]
  Vbeta    <- Vbeta[keep, keep, drop = FALSE]
  
  # Difference vector (intercept cancels)
  d <- c(0, as.numeric(Xa - Xb))
  d <- d[keep]
  
  gain_hat <- -sum(d * beta_hat) # add negative for minimizing
  
  var_gain <- as.numeric(t(d) %*% Vbeta %*% d)
  var_gain <- max(var_gain, 1e-12)
  
  se_gain <- sqrt(var_gain)
  
  z <- gain_hat / se_gain
  
  ei <- se_gain * dnorm(z) +
    gain_hat * pnorm(z)
  
  return(ei)
}

weighted_ei <- function(model_list,
                        weights,
                        neighbor_order,
                        incumbent_order){
  
  eis <- numeric(length(model_list))
  
  for(j in seq_along(model_list)){
    
    eis[j] <- model_ei(
      fit             = model_list[[j]]$fit,
      model_type      = model_list[[j]]$type,
      neighbor_order  = neighbor_order,
      incumbent_order = incumbent_order
    )
    
  }
  
  sum(weights * eis)
}

# Two-Phase Online OofA Experiment via Mixture-of-Experts
# This implements Algorithm 1 with simulated responses
# the goal is to identify an order that MINIMIZES the response
# Inputs:
# ---- n1: phase I sample size
# ---- n_max: maximum total sample size
# ---- m: number of components
# ---- models:  a vector of model names, options are: "pwo" (pairwise order),
#     "cp" (component-position), "lcpf" (first order flexible position)
#      default is models = c("pwo","cp","lcpf")
# ---- learning_rate: learning rate (positive number), default is 0.05
# ---- true_model: the true model used to generate the response. Default is "lcpf+qcpf"
# ---- sigma: error standard deviation (must be positive)
# ---- tol: convergence criteria to stop early (default is tol = 0.05). The method stops 
#      in Phase II when the weighted expected improvement is less than tol.
# ---- n_check: If the method is stuck at the same best order for this many
#      iterations, we stop sampling. Default is n_check = 3.
online_oofa_sampling <- function(n1, n_max, m, models = c("pwo","cp","lcpf", "lcpf+qcpf"),
                                 learning_rate = 0.05, true_model = "lcpf+qcpf", 
                                 sigma = 0.05, tol = 0.05, n_check = 3){
  
  # For now, Phase I uses random sampling
  # in the future we will use an OofA OA
  sampled_orders = matrix(NA, nrow = n_max, ncol = m)
  y = numeric(n_max)
  # if(n1 %% 12 != 0 && m < 9){
  #   stop("n1 needs to be a multiple of 12")
  # }
  if(m < 9 || n1 %% 12 != 0){ # if we have m < 9 and can find an OofA-OA, do it
    OofA_OA = as.matrix(des.cat(m,n1,"OofA-OA")$dem)    
  } else{
    # for m >= 9 or if OofA-OA does not exist, use algorithmic search
    D1 = matrix(NA, nrow = n1, ncol = m)
    for(j in 1:nrow(D1)){
      D1[j,] = sample(1:m)
    }
    # check that initial design is not singular
    found_init_design = FALSE
    while(!found_init_design){
      X_PWO = gen.pf(D1, "pwof")
      # found_init_design = is.finite(kappa(t(X_PWO)%*%X_PWO/n1))
      found_init_design = ncol(X_PWO) == qr(X_PWO)$rank
    }
    OofA_OA = find_opt_design_SA(init_runs = D1, n_iter = 100)
    
  }
  

  for(i in 1:n1){
    
    sampled_order = OofA_OA[i,]
    sampled_orders[i,] = sampled_order
    # sampled_order = sample(1:m)
    # sampled_orders[i,] = sampled_order
    y[i] = generate_response(sampled_order, true_model = true_model, sigma = sigma)
  }
  
  # find the best order so far
  best_index = which.min(y[1:n1])
  best_order = sampled_orders[best_index,]
  best_y = y[best_index]
  
  M = length(models)
  weights = rep(1/M, M)
  
  # Phase II
  sample_size = n1
  counter = 0 # count the number of times we stay at prev_best_order
  for(i in (n1+1):n_max){
    
    fitted_models = fit_all_models(models, D = sampled_orders[1:(i-1),], y = y[1:(i-1)])
    neighbors = get_neighboring_orders(best_order)

    n_neighbors = nrow(neighbors)
    yhats = numeric(n_neighbors)
    for(k in 1:n_neighbors){
      yhat_models = predict_all_models(model_list = fitted_models, new_order = neighbors[k,])
      yhats[k] = sum(weights*yhat_models)
    }
    
    eis <- numeric(n_neighbors)
    
    for(k in 1:n_neighbors){
      
      eis[k] <- weighted_ei(
        model_list      = fitted_models,
        weights         = weights,
        neighbor_order  = neighbors[k,],
        incumbent_order = best_order
      )
      
    }
    
    delta <- max(eis)
    
    cat("Max weighted EI =", delta, "\n")
    
    if(delta < tol){
      break
    }
    
    best_index = which.max(eis)
    sampled_order = neighbors[best_index,]
    sampled_orders[i,] = sampled_order
    
    cat("Sampled order:", sampled_order, "\n")
    
    y[i] = generate_response(sampled_order = sampled_order, true_model = true_model,
                             sigma = sigma)
    
    
    # find the loss
    yhats_models = predict_all_models(model_list = fitted_models, new_order = sampled_order)
    losses = sapply(yhats_models, function(x) (y[i] - x)^2  )
    
    # update the weights
    weights = weights*exp(-learning_rate*losses)
    
    # normalize the weights
    weights = weights/sum(weights)

    
    # save the best order from the last iteration
    prev_best_order = best_order
    
    # update the best order
    best_index = which.min(y[1:i])
    best_order = sampled_orders[best_index,]
    best_y = y[best_index]
    
    # increase sample size counter
    sample_size = sample_size + 1
   
    # check if we changed orders
    if(all(prev_best_order == best_order)){
      counter = counter + 1
    } else{
      counter = 0
    }
    
    # if we have been "stuck" for more than n_check iterations, stop
    if(counter >= n_check){
      break
    }
 
    
  }
  
  best_est_y = generate_response(sampled_order = best_order,
                                 true_model = true_model, sigma = 0)
  
  return(list(best_order = best_order,
              best_y = best_est_y,
              sample_size = sample_size))
  
  
}

# Recursive partitioning step for Lin & Chen (2023)
recursive_partition <- function(active_components,
                                current_order,
                                crit_val,
                                true_model,
                                sigma,
                                sample_size){
  
  k <- length(active_components)
  
  # base case
  if(k <= 1){
    return(list(
      order = current_order,
      sample_size = sample_size
    ))
  }
  
  pivot <- active_components[1]
  
  left_set   <- c()
  middle_set <- pivot
  right_set  <- c()
  
  # Construct the k orders used in Algorithm 2
  orders <- matrix(rep(current_order, k),
                   nrow = k,
                   byrow = TRUE)
  
  y <- numeric(k)
  
  cat("Order:", paste(current_order, collapse = " "), "\n")
  
  # Step 1: baseline order
  y[1] <- generate_response(
    current_order,
    true_model = true_model,
    sigma = sigma
  )
  
  sample_size <- sample_size + 1
  
  # Step 2: swap pivot with each component
  for(t in 2:k){
    
    comp <- active_components[t]
    
    pivot_pos <- which(orders[t, ] == pivot)
    comp_pos  <- which(orders[t, ] == comp)
    
    orders[t, c(pivot_pos, comp_pos)] <-
      orders[t, c(comp_pos, pivot_pos)]
    
    cat("Order:", paste(orders[t,], collapse = " "), "\n")
    y[t] <- generate_response(
      orders[t, ],
      true_model = true_model,
      sigma = sigma
    )
    
    sample_size <- sample_size + 1
  }
  
  # Step 3: Partition using y_{t-1} - y_t
  for(t in 2:k){
    
    comp <- active_components[t]
    
    diff <- y[t - 1] - y[t]
    
    if(diff > crit_val){
      
      left_set <- c(left_set, comp)
      
    } else if(diff < -crit_val){
      
      right_set <- c(right_set, comp)
      
    } else{
      
      middle_set <- c(middle_set, comp)
    }
  }
  
  # Update current ordering
  # new_suborder <- c(left_set, middle_set, right_set)
  # 
  # idx <- match(active_components, current_order)
  # 
  # current_order[idx] <- new_suborder
  
  # Randomly arrange components in the middle set
  if(length(middle_set) > 1){
    middle_set <- sample(middle_set)
  }
  
  # Update current ordering
  new_suborder <- c(left_set, middle_set, right_set)
  
  idx <- match(active_components, current_order)
  
  current_order[idx] <- new_suborder
  
  cat("Left :", paste(left_set, collapse = " "), "\n")
  cat("Middle:", paste(middle_set, collapse = " "), "\n")
  cat("Right :", paste(right_set, collapse = " "), "\n\n")
  
  # Recursive call on left subset
  if(length(left_set) > 1){
    
    left_result <- recursive_partition(
      active_components = left_set,
      current_order     = current_order,
      crit_val          = crit_val,
      true_model        = true_model,
      sigma             = sigma,
      sample_size       = sample_size
    )
    
    current_order <- left_result$order
    sample_size   <- left_result$sample_size
  }
  
  # Recursive call on right subset
  if(length(right_set) > 1){
    
    right_result <- recursive_partition(
      active_components = right_set,
      current_order     = current_order,
      crit_val          = crit_val,
      true_model        = true_model,
      sigma             = sigma,
      sample_size       = sample_size
    )
    
    current_order <- right_result$order
    sample_size   <- right_result$sample_size
  }
  
  return(list(
    order = current_order,
    sample_size = sample_size
  ))
}

# quicksort_algorithm
quicksort_algorithm <- function(n1,
                                m,
                                true_model = "lcpf+qcpf",
                                sigma = 0.05,
                                pivot = 1){
  
  initial_order <- c(
    pivot,
    sample(setdiff(1:m, pivot))
  )
  
  # Estimate variance once
  starting_ys <- numeric(n1)
  
  for(i in 1:n1){
    starting_ys[i] <- generate_response(
      initial_order,
      true_model = true_model,
      sigma = sigma
    )
  }
  
  sigma2_hat <- var(starting_ys)
  
  crit_val <- qnorm(0.975) * sqrt(2 * sigma2_hat)
  
  result <- recursive_partition(
    active_components = initial_order,
    current_order     = initial_order,
    crit_val          = crit_val,
    true_model        = true_model,
    sigma             = sigma,
    sample_size       = n1
  )
  
  estimated_order <- result$order
  
  estimated_y <- generate_response(
    estimated_order,
    true_model = true_model,
    sigma = 0
  )
  
  return(list(
    best_order = estimated_order,
    best_y     = estimated_y,
    sample_size     = result$sample_size
  ))
}

hamming_distance <- function(x,y){
  sum(x != y)
}


# Takes a model matrix X
# returns D_efficiency
D_efficiency <- function(X){
  k <- ncol(X)
  N <- nrow(X)
  return( (1/N)*det(t(X)%*%X)^(1/k) )
}

## find_opt_design_SA
## input: 
## init_runs: a matrix of initial runs of the desired size
## n_iter: number of iterations (default 1000)
## trace: if TRUE, prints iteration number (default TRUE)
## show_plot: if TRUE, plots the D-efficiency by iteration (default FALSE)
## criterion: criterion used for efficiency of design, by default, this is D
##### Each iteration selects a random row of the design matrix
##### And randomly swaps two adjacent elements 
##### Compares the d-efficiency of current with new one 
find_opt_design_SA <- function(init_runs, n_iter = 1000, trace = TRUE,
                                 show_plot = FALSE, criterion = "D",
                                 temp_init = 1){
  m <- ncol(init_runs)
  current_design <- init_runs
  size <- nrow(init_runs)
  Z_current <- gen.pf(init_runs, "pwof")
  C_current <- solve(t(Z_current)%*%Z_current)
  if(criterion == "D"){
    det_current <- det(t(Z_current)%*%Z_current)
    d_current <- D_efficiency(Z_current)
  }
  if(criterion == "A"){
    d_current <- sum(diag(C_current))
  }
  d_best <- d_current
  Z_best <- Z_current
  best_design <- current_design
  ds <- numeric(n_iter)
  for(t in 1:n_iter){
    ds[t] <- d_current
    if(trace){
      print(paste("Beginning iteration ",t,sep=""))
    }
    new_design <- current_design
    temp <- temp_init/(t+1)
    # for each iteration, randomly select a row (1...size)
    # and a column in the row (col from 1...m-1)
    i <- sample(1:size,1)  # choose random row
    j <- sample(1:(m-1),1) # choose random column
    
    tmp_design <- current_design
    tmp_design[i,j] <- current_design[i,j+1]
    tmp_design[i,j+1] <- current_design[i,j]
    Z_new <- gen.pf(tmp_design, "pwof")
    F1 <- cbind(Z_new[i,], -Z_current[i,])
    F2 <- cbind(Z_new[i,], Z_current[i,])
    tmp0 <- diag(2) + t(F2)%*%C_current%*%F1
    if(rcond(tmp0) < 1e-14){
      next
    }
    tmp <- solve(tmp0)
    C_new <- C_current - C_current%*%F1%*%tmp%*%t(F2)%*%C_current
    if(criterion == "D"){
      f1 <- F2[,2]
      f2 <- F2[,1]
      v12 <- t(f1)%*%C_current%*%f2
      v1 <- t(f1)%*%C_current%*%f1
      v2 <- t(f2)%*%C_current%*%f2
      delta <- 1 + v2 - v1 + v12^2 - v2*v1
      det_new <- det_current*delta
      d_new <-  det_new^(1/ncol(Z_new))/nrow(Z_new)
      if(!is.na(d_new) & !is.na(d_current)){
        aprob <- exp(-(d_current - d_new)/temp)
      }
      else{
        aprob <- -1
      }
      
    }
    if(criterion == "A"){
      d_new <- sum(diag(C_new))
      if(!is.na(d_new) & !is.na(d_current)){
        aprob <- exp(-(d_new - d_current)/temp)
      }
      else{
        aprob <- -1
      }
    }
    
    print(paste("Acceptance probability: ", aprob, sep = ""))
    
    if( aprob > runif(1)){
      print("==========ACCEPTED NEW DESIGN===============")
      current_design <- new_design
      Z_current <- Z_new
      d_current <- d_new
      C_current <- C_new
      if(criterion == "D"){
        det_current <- det_new
      }
      
      if((d_new > d_best) & criterion == "D"){
        best_design <- new_design
        Z_best <- Z_new
        d_best <- d_new
      }
      if((d_new < d_best) & criterion == "A"){
        best_design <- new_design
        Z_best <- Z_new
        d_best <- d_new
      }
    }
    
    # if(aprob < 1e-20){
    #   break
    # }
    
    
  }
  
  if(show_plot){
    if(criterion == "D"){
      plot(1:n_iter,ds,type="l", xlab = "iteration", ylab = "D-efficiency")
    }
    if(criterion == "A"){
      plot(1:n_iter,ds,type="l", xlab = "iteration", ylab = "A-efficiency")
    }
  }
  
  return(best_design)
  
}

################################################################################
######################### Section 4 Functions ##################################
################################################################################

## find_opt_design_SA_reduced
## input: 
## init_runs: a matrix of initial runs of the desired size
## pwo_cols: column indices of the PWO coefficients included in reduced model
## n_iter: number of iterations (default 1000)
## trace: if TRUE, prints iteration number (default TRUE)
## show_plot: if TRUE, plots the D-efficiency by iteration (default FALSE)
## criterion: criterion used for efficiency of design, by default D
find_opt_design_SA_reduced <- function(init_runs, pwo_cols,
                                       n_iter = 1000, trace = TRUE,
                                       show_plot = FALSE, criterion = "D",
                                       temp_init = 1){
  
  m <- ncol(init_runs)
  current_design <- init_runs
  size <- nrow(init_runs)
  
  # Reduced PWO matrix
  Z_current <- gen.pf(init_runs, "pwof")
  Z_current <- Z_current[, pwo_cols, drop = FALSE]
  
  # Check estimability of initial design
  if(qr(cbind(1, Z_current))$rank < ncol(Z_current) + 1){
    stop("Initial design is not estimable under the reduced PWO model.")
  }
  
  C_current <- solve(t(Z_current) %*% Z_current)
  
  if(criterion == "D"){
    det_current <- det(t(Z_current) %*% Z_current)
    d_current <- D_efficiency(Z_current)
  }
  
  if(criterion == "A"){
    d_current <- sum(diag(C_current))
  }
  
  d_best <- d_current
  Z_best <- Z_current
  best_design <- current_design
  ds <- numeric(n_iter)
  
  for(t in 1:n_iter){
    
    ds[t] <- d_current
    
    if(trace){
      print(paste("Beginning iteration ", t, sep = ""))
    }
    
    temp <- temp_init/(t + 1)
    
    # randomly select row and adjacent positions
    i <- sample(1:size, 1)
    j <- sample(1:(m-1), 1)
    
    # generate proposed design
    new_design <- current_design
    new_design[i,j] <- current_design[i,j+1]
    new_design[i,j+1] <- current_design[i,j]
    
    # Reduced PWO matrix for proposed design
    Z_new <- gen.pf(new_design, "pwof")
    Z_new <- Z_new[, pwo_cols, drop = FALSE]
    
    F1 <- cbind(Z_new[i,], -Z_current[i,])
    F2 <- cbind(Z_new[i,], Z_current[i,])
    
    tmp0 <- diag(2) + t(F2)%*%C_current%*%F1
    if(rcond(tmp0) < 1e-14){
      next
    }
    tmp <- solve(tmp0)
    
    C_new <- C_current -
      C_current %*% F1 %*% tmp %*% t(F2) %*% C_current
    
    if(criterion == "D"){
      
      f1 <- F2[,2]
      f2 <- F2[,1]
      
      v12 <- t(f1) %*% C_current %*% f2
      v1 <- t(f1) %*% C_current %*% f1
      v2 <- t(f2) %*% C_current %*% f2
      
      delta <- 1 + v2 - v1 + v12^2 - v2*v1
      
      det_new <- det_current * delta
      
      d_new <- det_new^(1/ncol(Z_new))/nrow(Z_new)
      
      if(!is.na(d_new) & !is.na(d_current)){
        aprob <- exp(-(d_current - d_new)/temp)
      } else{
        aprob <- -1
      }
    }
    
    if(criterion == "A"){
      
      d_new <- sum(diag(C_new))
      
      if(!is.na(d_new) & !is.na(d_current)){
        aprob <- exp(-(d_new - d_current)/temp)
      } else{
        aprob <- -1
      }
    }
    
    if(trace){
      print(paste("Acceptance probability: ", aprob, sep = ""))
    }
    
    if(aprob > runif(1)){
      
      if(trace){
        print("==========ACCEPTED NEW DESIGN===============")
      }
      
      current_design <- new_design
      Z_current <- Z_new
      d_current <- d_new
      C_current <- C_new
      
      if(criterion == "D"){
        det_current <- det_new
      }
      
      if((d_new > d_best) & criterion == "D"){
        best_design <- new_design
        Z_best <- Z_new
        d_best <- d_new
      }
      
      if((d_new < d_best) & criterion == "A"){
        best_design <- new_design
        Z_best <- Z_new
        d_best <- d_new
      }
    }
  }
  
  if(show_plot){
    
    if(criterion == "D"){
      plot(1:n_iter, ds, type = "l",
           xlab = "iteration",
           ylab = "D-efficiency")
    }
    
    if(criterion == "A"){
      plot(1:n_iter, ds, type = "l",
           xlab = "iteration",
           ylab = "A-efficiency")
    }
  }
  
  return(best_design)
}

# helper function to get the cost of an order in the job scheduling problem (Section 4)
# sigma controls the noise of the error (sigma = 0 gives deterministic responses)
get_order_cost <- function(order, p_vec, omega_vec, sigma = 0.5){
  
  m = length(order)
  proc_time = cumsum(p_vec[order])
  order_cost = sum(omega_vec * (proc_time)^2) + rnorm(1, mean = 0, sd = sigma)
  return(order_cost)
}

# Two-Phase Online OofA Experiment via Mixture-of-Experts
# This implements Algorithm 1 for the Job Scheduling Example (Section 4)
# the goal is to identify an order that MINIMIZES the response
# Inputs:
# ---- initial_design: Phase I design (n1 rows, to be found by SA algorithm)
# ---- n_max: maximum total sample size
# ---- m: number of components
# ---- models:  a vector of model names, options are: "pwo" (full pairwise order),
#     "pwo1", "pwo2", "pwo3" (reduced models defined in Section 4), 
#     "cp" (component-position), "lcpf" (first order flexible position)
#      default is models = c("pwo1","pwo2","pwo3","lcpf", "lcpf+qcpf")
# ---- p_vec: vector of processing times for each job (must be non-negative)
# ---- omega_vec: vector of weights for each job (must be non-negative)
# ---- learning_rate: learning rate (positive number), default is 0.05
# ---- sigma: error standard deviation (must be positive)
# ---- tol: convergence criteria to stop early (default is tol = 0.05). The method stops 
#      in Phase II when the weighted expected improvement is less than tol.
# ---- n_check: If the method is stuck at the same best order for this many
#      iterations, we stop sampling. Default is n_check = 3.
online_oofa_sampling_Section4 <- function(initial_design, n_max, m, models = c("pwo1","pwo2","pwo3","lcpf", "lcpf+qcpf"),
                                 p_vec, omega_vec, learning_rate = 0.05, sigma = 0.05, tol = 0.05, n_check = 5){
  
  # For now, Phase I uses random sampling
  # in the future we will use an OofA OA
  sampled_orders = matrix(NA, nrow = n_max, ncol = m)
  y = numeric(n_max)
  OofA_OA = initial_design
  
  for(i in 1:n1){
    
    sampled_order = OofA_OA[i,]
    sampled_orders[i,] = sampled_order
    # sampled_order = sample(1:m)
    # sampled_orders[i,] = sampled_order
    y[i] = get_order_cost(sampled_order, p_vec = p_vec, omega_vec = omega_vec, 
                          sigma = sigma)
  }
  
  # find the best order so far
  best_index = which.min(y[1:n1])
  best_order = sampled_orders[best_index,]
  best_y = y[best_index]
  
  M = length(models)
  weights = rep(1/M, M)
  
  # Phase II
  sample_size = n1
  counter = 0 # count the number of times we stay at prev_best_order
  for(i in (n1+1):n_max){
    
    fitted_models = fit_all_models(models, D = sampled_orders[1:(i-1),], y = y[1:(i-1)])
    neighbors = get_neighboring_orders(best_order)
    
    n_neighbors = nrow(neighbors)
    yhats = numeric(n_neighbors)
    for(k in 1:n_neighbors){
      yhat_models = predict_all_models(model_list = fitted_models, new_order = neighbors[k,])
      yhats[k] = sum(weights*yhat_models)
    }
    
    eis <- numeric(n_neighbors)
    
    for(k in 1:n_neighbors){
      
      eis[k] <- weighted_ei(
        model_list      = fitted_models,
        weights         = weights,
        neighbor_order  = neighbors[k,],
        incumbent_order = best_order
      )
      
    }
    
    delta <- max(eis)
    
    cat("Max weighted EI =", delta, "\n")
    
    if(delta < tol){
      break
    }
    
    best_index = which.max(eis)
    sampled_order = neighbors[best_index,]
    sampled_orders[i,] = sampled_order
    
    cat("Sampled order:", sampled_order, "\n")
    
    y[i] = get_order_cost(sampled_order, p_vec = p_vec, omega_vec = omega_vec, 
                           sigma = sigma)
    
    
    # find the loss
    yhats_models = predict_all_models(model_list = fitted_models, new_order = sampled_order)
    losses = sapply(yhats_models, function(x) (y[i] - x)^2  )
    
    # update the weights
    weights = weights*exp(-learning_rate*losses)
    
    # normalize the weights
    weights = weights/sum(weights)
    
    
    # save the best order from the last iteration
    prev_best_order = best_order
    
    # update the best order
    best_index = which.min(y[1:i])
    best_order = sampled_orders[best_index,]
    best_y = y[best_index]
    
    # increase sample size counter
    sample_size = sample_size + 1
    
    # check if we changed orders
    if(all(prev_best_order == best_order)){
      counter = counter + 1
    } else{
      counter = 0
    }
    
    # if we have been "stuck" for more than n_check iterations, stop
    if(counter >= n_check){
      break
    }
    
    
  }
  
  best_est_y = get_order_cost(sampled_order, p_vec = p_vec, omega_vec = omega_vec, 
                               sigma = 0)
    
  
  return(list(best_order = best_order,
              best_y = best_est_y,
              sample_size = sample_size))
  
  
}



# Recursive partitioning step for Lin & Chen (2023), adapted to Section 4 
# job scheduling example
recursive_partition_Section4 <- function(active_components,
                                current_order,
                                crit_val,
                                p_vec, omega_vec,
                                sigma,
                                sample_size){
  
  k <- length(active_components)
  
  # base case
  if(k <= 1){
    return(list(
      order = current_order,
      sample_size = sample_size
    ))
  }
  
  pivot <- active_components[1]
  
  left_set   <- c()
  middle_set <- pivot
  right_set  <- c()
  
  # Construct the k orders used in Algorithm 2
  orders <- matrix(rep(current_order, k),
                   nrow = k,
                   byrow = TRUE)
  
  y <- numeric(k)
  
  cat("Order:", paste(current_order, collapse = " "), "\n")
  
  # Step 1: baseline order
  y[1] <- get_order_cost(order = current_order, 
                         p_vec = p_vec, omega_vec = omega_vec, sigma = sigma)

  sample_size <- sample_size + 1
  
  # Step 2: swap pivot with each component
  for(t in 2:k){
    
    comp <- active_components[t]
    
    pivot_pos <- which(orders[t, ] == pivot)
    comp_pos  <- which(orders[t, ] == comp)
    
    orders[t, c(pivot_pos, comp_pos)] <-
      orders[t, c(comp_pos, pivot_pos)]
    
    cat("Order:", paste(orders[t,], collapse = " "), "\n")
    y[t] <- get_order_cost(order = orders[t,], 
                           p_vec = p_vec, omega_vec = omega_vec, sigma = sigma)
      
    
    sample_size <- sample_size + 1
  }
  
  # Step 3: Partition using y_{t-1} - y_t
  for(t in 2:k){
    
    comp <- active_components[t]
    
    diff <- y[t - 1] - y[t]
    
    if(diff > crit_val){
      
      left_set <- c(left_set, comp)
      
    } else if(diff < -crit_val){
      
      right_set <- c(right_set, comp)
      
    } else{
      
      middle_set <- c(middle_set, comp)
    }
  }
  
  
  # Randomly arrange components in the middle set
  if(length(middle_set) > 1){
    middle_set <- sample(middle_set)
  }
  
  # Update current ordering
  new_suborder <- c(left_set, middle_set, right_set)
  
  idx <- match(active_components, current_order)
  
  current_order[idx] <- new_suborder
  
  cat("Left :", paste(left_set, collapse = " "), "\n")
  cat("Middle:", paste(middle_set, collapse = " "), "\n")
  cat("Right :", paste(right_set, collapse = " "), "\n\n")
  
  # Recursive call on left subset
  if(length(left_set) > 1){
    
    left_result <- recursive_partition_Section4(
      active_components = left_set,
      current_order     = current_order,
      crit_val          = crit_val,
      p_vec = p_vec, omega_vec = omega_vec,
      sigma             = sigma,
      sample_size       = sample_size
    )
    
    current_order <- left_result$order
    sample_size   <- left_result$sample_size
  }
  
  # Recursive call on right subset
  if(length(right_set) > 1){
    
    right_result <- recursive_partition_Section4(
      active_components = right_set,
      current_order     = current_order,
      crit_val          = crit_val,
      p_vec = p_vec, omega_vec = omega_vec,
      sigma             = sigma,
      sample_size       = sample_size
    )
    
    current_order <- right_result$order
    sample_size   <- right_result$sample_size
  }
  
  return(list(
    order = current_order,
    sample_size = sample_size
  ))
}

# quicksort_algorithm, for the job scheduling problem in Section 4
quicksort_algorithm_Section4 <- function(n1,
                                m,
                                p_vec, omega_vec,
                                sigma = 0.5,
                                pivot = 1){
  
  initial_order <- c(
    pivot,
    sample(setdiff(1:m, pivot))
  )
  
  # Estimate variance once
  starting_ys <- numeric(n1)
  
  for(i in 1:n1){
    starting_ys[i] <- get_order_cost(initial_order,
                                     p_vec = p_vec, omega_vec = omega_vec,
                                     sigma = sigma)
  }
  
  sigma2_hat <- var(starting_ys)
  
  crit_val <- qnorm(0.975) * sqrt(2 * sigma2_hat)
  
  result <- recursive_partition_Section4(
    active_components = initial_order,
    current_order     = initial_order,
    crit_val          = crit_val,
    p_vec = p_vec, omega_vec = omega_vec,
    sigma             = sigma,
    sample_size       = n1
  )
  
  estimated_order <- result$order
  
  estimated_y <- get_order_cost(estimated_order,
                                p_vec = p_vec, omega_vec = omega_vec,
                                sigma = 0)
  
  return(list(
    best_order = estimated_order,
    best_y     = estimated_y,
    sample_size     = result$sample_size
  ))
}
