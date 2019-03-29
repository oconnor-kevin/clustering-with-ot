# TITLE: hybrid_wasserstein.R
# AUTHOR: Kevin O'Connor
# DATE CREATED: 3/28/19

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


get_tangent_approx <- function(dist_df){
  # Computes tangent approximation of wasserstein distance.
  
  m <- group_by(dist_df, class) %>% summarize(count = n()) %>% 
    dplyr::select(count) %>% as.vector() %>% min()
  R <- get_reference_measure(dist_df)
  
  # Get transport maps.
  for(c in unique(dist_df$class)){
    dist <- filter(dist_df, class == c) %>% dplyr::select(-class) %>% 
      as.matrix() %>% t()  # p x n
    
  }
  rkde(m, R)
}


hybrid_wasserstein <- function(x, y){}