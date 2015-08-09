#' @name dpCluster
#'
#' @title dpCluster
#'
#' @param x data.frame or matrix
#'
#' @param percent method `gaussian` and `withinDc` use the nearest `percent` of points to estimate `dc`;
#'
#'        method `neighbors` use the nearest `percent` of points to estimate density directly.
#'
#'        percent is suggested between 0.01 and 0.02
#'
#' @param thres.rho threshold of rho to detect peaks, which is often get from "Decision Graph".
#'
#' @param thres.delta threshold of delta to detect peaks, which is often get from "Decision Graph".
#'
#' @param method there are three methods to estimate density.
#'
#'        `gaussian` : gaussian kernel density estimation $e^{(\frac{-1}{2}(\frac{\delta(x_{i}, x_{j})}{2dc})^{2})}$
#'
#'        `withinDc` : the density of $x_{i}$ is defined as the number of points within dc distance of $x_{i}$
#'
#'        `neighbors` : the density of $x_{i}$ is defined as the reciprocal of the mean distance of $x_{i}$'s neighbors
#'
#' @param threads number of threads to use for task scheduling
#'
#' @param halo.detection logical, whether to remove the potential noise point with low density
#'
#' @export
#'
dpCluster <- function(x,
                      percent = 0.01,
                      thres.rho = NULL,
                      thres.delta = NULL,
                      similarity = c("euclidean", "SNN"),
                      method = c("gaussian", "withinDc", "neighbors"),
                      threads = 4,
                      halo.detection = TRUE)
{
  if(!is.data.frame(x) && !is.matrix(x)) stop("x is not a data.frame or a matrix!\n")
  if(is.data.frame(x)) x <- as.matrix(x)
  method <- match.arg(method)
  similarity <- match.arg(similarity)
  threads <- min(threads, defaultNumThreads())

  params <- get_rho_delta(data = x, similarity = similarity, method = method, percent = percent, threads = threads)

  while(is.null(thres.rho) || is.null(thres.delta))
  {
    plot(params$rho, params$delta, xlab = "rho", ylab = "delta", main = "The Decision Graph")
    cat("Please click on plot to select threashold\n")
    click <- locator(1)
    if(is.null(thres.rho)) thres.rho <- click$x
    if(is.null(thres.delta)) thres.delta <- click$y
  }

  peaks <- which(params$rho >= thres.rho & params$delta >= thres.delta)
  plot(params$rho, params$delta, xlab = "rho", ylab = "delta", main = "The Decision Graph")
  points(params$rho[peaks], params$delta[peaks], col = 2:(1+length(peaks)), pch = 19)

  if(length(peaks) == 1) stop("there's only one cluster!\n")
  res <- dpCluster_cpp(parameters = params, peaks = peaks - 1, use_halo = halo.detection)

  res$peaks <- res$peaks + 1
  class(res) <- c("partition", "dpCluster")
  res
}
