# TITLE: hybrid_wasserstein.R
# AUTHOR: Kevin O'Connor
# DATE CREATED: 3/28/19
# DATE MODIFIED: 3/29/19

wasserstein_normal_mu_sigma <- function(mu_x, sigma_x, mu_y, sigma_y){
  return(t(mu_x - mu_y) %*% (mu_x - mu_y) + tr(sigma_x) + tr(sigma_y)
         - 2 * tr(sqrtm(sqrtm(sigma_x)$B %*% sigma_y %*% sqrtm(sigma_x)$B)$B))
}


wasserstein_normal_x_y <- function(x, y){
  # Computes the wasserstein-2^2 distance between z_x and z_y where z_x is the 
  #  Gaussian distribution with mean(z_x) = mean(x) and cov(z_x) = cov(x).
  #
  # Args:
  #  x: p x n dimensional data matrix.
  #  y: p x n dimensional data matrix.
  #
  # Returns:
  #  Wasserstein distance between z_x and z_y.
  require(pracma)
  require(psych)
  
  mu_x <- rowMeans(x)
  mu_y <- rowMeans(y)
  sigma_x <- cov(t(x))
  sigma_y <- cov(t(y))
  
  return(wasserstein_normal_mu_sigma(mu_x, sigma_x, mu_y, sigma_y))
}


wasserstein_normal <- function(dist_df){
  # Computes the matrix of wasserstein-2^2 distances between moment-matched
  #  Gaussians of distributions in dist_df.
  #
  # Args:
  #  dist_df: Dataframe with "class" column indicating the distribution and at 
  #   least one additional column with numeric observations from each dist.
  #
  # Returns:
  #  Matrix of normal wasserstein distances between distributions in dist_df.
  
  wass_norm_dist_x_y <- function(x, y){
    if(x == y){0}
    else{
      dist_x <- filter(dist_df, class == x) %>% dplyr::select(-class) %>% 
        as.matrix() %>% t()
      dist_y <- filter(dist_df, class == y) %>% dplyr::select(-class) %>% 
        as.matrix() %>% t()
      wasserstein_normal_x_y(dist_x, dist_y) %>% as.numeric()
    }
  }
  wass_norm_dist_x_y <- Vectorize(wass_norm_dist_x_y)
  outer(unique(dist_df$class), unique(dist_df$class), FUN = wass_norm_dist_x_y)
}


get_reference_measure <- function(dist_df){
  # Computes reference measure given a dataframe of distributions.
  #
  # Args:
  #  dist_df: Dataframe with "class" column indicating the distribution and at 
  #   least one additional column with numeric observations from each dist.
  #
  # Returns:
  #  Reference distribution constructed from input distributions.
  require(MASS)
  require(pracma)
  require(tidyverse)
  require(ks)
  select <- dplyr::select
  
  ns <- group_by(dist_df, class) %>% summarize(count = n()) %>% 
    select(count) %>% unlist() %>% as.vector()
  pis <- ns / sum(ns)
  
  # Center and standardize distributions.
  std_dist <- c()
  for(c in unique(dist_df$class)){
    dist <- filter(dist_df, class == c) %>% select(-class) %>% 
      as.matrix() %>% t()  # p x n
    dist <- dist - rowMeans(dist)
    dist <- sqrtm(ginv(cov(t(dist))))$B %*% dist
    std_dist <- rbind(std_dist, t(dist))
  }
  # Get kernel density estimate.
  return(kde(std_dist))
}


get_tangent_approx_mat <- function(dist_df, m = "None"){
  # Computes matrix of tangent approximation distances between distributions in
  #  dist_df.
  #
  # Args:
  # dist_df: Dataframe with column "class" and at least one other column 
  #  corresponding to the points generating the distributions.
  # m: Scalar giving the subsample size. If "None", choose as minimum 
  #  distribution size. 
  #
  # Returns:
  #  Matrix giving the tangent approximation distance between each distribution.
  
  require(clue)
  
  # Get reference measure and sample from it.
  if(m == "None"){
    m <- group_by(dist_df, class) %>% summarize(count = n()) %>% 
      dplyr::select(count) %>% as.vector() %>% min()
  }
  R <- get_reference_measure(dist_df)
  R_samp <- rkde(m, R) %>% t()  # p x m
  
  # Get transport maps.
  T_list <- list()
  T_Us_list <- list()
  for(i in 1:length(unique(dist_df$class))){
    c <- unique(dist_df$class)[i]
    dist <- filter(dist_df, class == c) %>% dplyr::select(-class) %>% 
      as.matrix() %>% t()  # p x n
    dist <- dist[, sample.int(ncol(dist), m)]  # p x m
    # Solve for distance matrix.
    l2_dist_R_samp <- function(x, y){
      (dist[, x] - R_samp[, y])^2 %>% sum()
    }
    l2_dist_R_samp <- Vectorize(l2_dist_R_samp)
    distances <- outer(1:m, 1:m, FUN = l2_dist_R_samp)
    # Solve for optimal transport map.
    T_list[[i]] <- solve_LSAP(distances)
    # Compute T(U_s) for s = 1, ..., m.
    ## Compute 1 nearest neighbors of U_s to distribution.
    U_s_nn <- apply(R_samp, 2, function(r){
      apply(dist, 2, function(d){sum((r - d)^2)}) %>% which.min()
    })
    ## Get T(U_s) via 1 nearest neighbors.
    T_Us_list[[i]] <- R_samp[, T_list[[i]][U_s_nn]]
  }
  
  # Compute distances between T_j(U_s) and T_k(U_s).
  l2_dist_Tj_Tk <- function(x, y){sum((x - y)^2)}
  l2_dist_Tj_Tk <- Vectorize(l2_dist_Tj_Tk)
  outer(T_Us_list, T_Us_list, FUN = l2_dist_Tj_Tk)
}


hybrid_wasserstein <- function(dist_df){
  # Compute hybrid wasserstein distance between distributions in dist_df.
  #
  # Args:
  #  dist_df: Dataframe with a column called "class" and at least one other ...
  
  return(wasserstein_normal(dist_df) + get_tangent_approx_mat(dist_df))
}


get_psis <- function(dist_df, R, R_samp){
  require(ks)
  require(clue)
  
  # Get transport maps.
  T_list <- list()
  T_Us_list <- list()
  for(c in unique(dist_df$class)){
    dist <- filter(dist_df, class == c) %>% dplyr::select(-class) %>% 
      as.matrix() %>% t()  # p x n
    dist <- dist[, sample.int(ncol(dist), m)]  # p x m
    # Solve for distance matrix.
    l2_dist_R_samp <- function(x, y){
      (dist[, x] - R_samp[, y])^2 %>% sum()
    }
    l2_dist_R_samp <- Vectorize(l2_dist_R_samp)
    distances <- outer(1:m, 1:m, FUN = l2_dist_R_samp)
    # Solve for optimal transport map.
    T_list[[c]] <- solve_LSAP(distances)
  }
  # Create psi functions.  
  psis <- lapply(unique(dist_df$class), function(c){
    T_j <- T_list[[c]]
    psi <- function(z){
      # Compute T(U_s) for s = 1, ..., m.
      ## Compute 1 nearest neighbors of U_s to distribution.
      dist <- filter(dist_df, class == c) %>% dplyr::select(-class) %>% 
        as.matrix() %>% t()
      nn <- apply(dist, 2, function(d){sum((z - d)^2)}) %>% which.min()
      ## Get T(U_s) via 1 nearest neighbors.
      return((R_samp[, T_j[nn]] - z) * sqrt(dkde(z, R)))
    }
    return(psi)
  })
  return(psis)
}


get_hybrid_barycenter <- function(dist_df, m = "None", n_iter = 10){
  require(MASS)
  require(ks)
  require(pracma)
  
  ## Get lambdas.
  ns <- group_by(dist_df, class) %>% summarize(count = n()) %>% 
    pull(count)
  lambdas <- ns / sum(ns)
  
  # Get mu_bar.
  mus <- group_by(dist_df, class) %>% summarise(mu1 = mean(dim1), mu2 = mean(dim2)) %>% dplyr::select(mu1, mu2)
  mu_bar <- c(mus$mu1 %*% lambdas, mus$mu2 %*% lambdas)
  
  # Get psi_bar.
  ## Get reference measure and sample from it.
  if(m == "None"){
    m <- group_by(dist_df, class) %>% summarize(count = n()) %>% 
      dplyr::select(count) %>% as.vector() %>% min()
  }
  R <- get_reference_measure(dist_df)
  R_samp <- rkde(m, R) %>% t()  # p x m
  ## Compute psi_bar.
  psi_bar <- function(z){
    psis_list <- get_psis(dist_df, R, R_samp)
    lam_psi_list <- lapply(1:length(psis_list), function(i){
      lambda_psi <- function(z){psis_list[[i]](z) * lambdas[i]}
      return(lambda_psi)
    })
    return(rowSums(sapply(lam_psi_list, function(lam_psi){lam_psi(z)})))
  }
  
  # Get sigma_bar.
  sigma <- dist_df %>% filter(class == 1) %>% dplyr::select(-class) %>% 
    as.matrix() %>% cov()
  for(i in 1:n_iter){
    term1 <- lapply(unique(dist_df$class), function(c){
      dist_cov <- dist_df %>% filter(class == c) %>% dplyr::select(-class) %>% 
        as.matrix() %>% cov()
      lambdas[c] * sqrtm(sqrtm(sigma)$B %*% dist_cov %*% sqrtm(sigma)$B)$B
    }) %>% 
      Reduce(f = `+`) %>%
      (function(x){x %*% x})

    sigma <- sqrtm(ginv(sigma))$B %*% term1 %*% sqrtm(ginv(sigma))$B
  }
  
  return(list("mu_bar" = mu_bar, "sigma_bar" = sigma, "psi_bar" = psi_bar))
}


get_inv_tmap <- function(dist_df, mu_bar, sigma_bar, psi_bar, m = "None"){
  
  # Get reference measure and sample from it.
  if(m == "None"){
    m <- group_by(dist_df, class) %>% summarize(count = n()) %>% 
      dplyr::select(count) %>% as.vector() %>% min()
  }
  R <- get_reference_measure(dist_df)
  R_samp <- rkde(m, R) %>% t()  # p x m
  
  T_bar <- function(z){
    psi_bar(z) / dkde(z, R) + z
  }
  
  return(T_bar)
  
  
}


#hybrid_wasserstein <- function(phi1, phi2){}