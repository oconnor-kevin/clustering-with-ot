# TITLE: hybrid_wasserstein.R
# AUTHOR: Kevin O'Connor
# DATE CREATED: 3/28/19
# DATE MODIFIED: 3/29/19

wasserstein_normal <- function(x, y){
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
  
  return(t(mu_x - mu_y) %*% (mu_x - mu_y) + tr(sigma_x) + tr(sigma_y)
         - 2 * tr(sqrtm(sqrtm(sigma_x)$B %*% sigma_y %*% sqrtm(sigma_x)$B)$B))
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


get_tangent_approx_mat <- function(dist_df){
  # Computes tangent approximation of wasserstein distance.
  
  require(clue)
  
  # Get reference measure and sample from it.
  m <- group_by(dist_df, class) %>% summarize(count = n()) %>% 
    dplyr::select(count) %>% as.vector() %>% min()
  R <- get_reference_measure(dist_df)
  R_samp <- rkde(m, R) %>% t()  # p x m
  
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
    distances <- outer(1:m, 1:m, FUN = l2_dist)
    # Solve for optimal transport map.
    T_list[[c]] <- solve_LSAP(distances)
    # Compute T(U_s) for s = 1, ..., m.
    ## Compute 1 nearest neighbors of U_s to distribution.
    U_s_nn <- apply(R_samp, 2, function(r){
      apply(dist, 2, function(d){sum((r - d)^2)}) %>% which.min()
    })
    ## Get T(U_s) via 1 nearest neighbors.
    T_Us_list[[c]] <- R_samp[, T_list[[c]][U_s_nn]]
  }
  
  # Compute distances between T_j(U_s) and T_k(U_s).
  l2_dist_Tj_Tk <- function(x, y){sum((x - y)^2)}
  l2_dist_Tj_Tk <- Vectorize(l2_dist_Tj_Tk)
  outer(T_Us_list, T_Us_list, FUN = l2_dist_Tj_Tk)
}


hybrid_wasserstein <- function(x, y){}