# This file was generated by Rcpp::compileAttributes
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

get_rho_delta <- function(data, similarity = "euclidean", method = "neighbors", percent = 0.01, threads = 4L, SNN_percent = 0.1) {
    .Call('dpCluster_get_rho_delta', PACKAGE = 'dpCluster', data, similarity, method, percent, threads, SNN_percent)
}

dpCluster_cpp <- function(parameters, peaks, use_halo = TRUE) {
    .Call('dpCluster_dpCluster_cpp', PACKAGE = 'dpCluster', parameters, peaks, use_halo)
}

SNN <- function(dist, size, percent) {
    .Call('dpCluster_SNN', PACKAGE = 'dpCluster', dist, size, percent)
}

