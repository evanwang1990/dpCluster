#include "dpCluster.h"
#include <utility>
#include <Rcpp.h>
//[[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace std;
using namespace Rcpp;
using namespace RcppParallel;

NumericVector Dist(const NumericMatrix &data);
pair<size_t, int> get_grainSize(const size_t size, double percent, const int threads);
double getDc(const NumericVector &dist, double percent, const int threads);
void getRho_withinDc(const NumericVector &dist, const size_t size, const double dc, NumericVector &rho);
void getRho_neighbors(const NumericVector &dist, const size_t size, const double percent, NumericVector &rho);
void getRho_gaussian(const NumericVector &dist, const size_t size, const double dc, NumericVector &rho);
void getDelta(const NumericVector &dist, const NumericVector &rho, size_t size, NumericVector &delta);

//[[Rcpp::export]]
List get_rho_delta(NumericMatrix data, String method = "neighbors", double percent = 0.01, int threads = 4)
{
  size_t size = data.nrow();
  NumericVector rho(size, 0.0);
  NumericVector delta(size, R_PosInf);
  double dc = 0.0;

  //calculate the distance
  NumericVector dist = Dist(data);

  //calculate rho and delta
  if(method == "neighbors")
  {
    getRho_neighbors(dist, size, percent, rho);
  }
  else
  {
    dc = getDc(dist, percent, threads);
    if(method == "withinDc")
    {
      getRho_withinDc(dist, size, dc, rho);
    }
    else
    {
      getRho_gaussian(dist, size, dc, rho);
    }
  }

  getDelta(dist, rho, size, delta);
  return List::create(
    Named("dist") = dist,
    Named("rho") = rho,
    Named("delta") = delta,
    Named("method") = method,
    Named("dc") = dc,
    Named("size") = size
  );
}

//[[Rcpp::export]]
List dpCluster_cpp(List parameters, IntegerVector peaks, bool use_halo = true)
{
  size_t size = parameters["size"];
  vector<int> cluster(size);
  NumericVector dist = parameters["dist"];
  double dc = parameters["dc"];
  NumericVector rho = parameters["rho"];
  NumericVector delta = parameters["delta"];
  String method = parameters["method"];

  if(!use_halo)
  {
    double dist_temp;
    for(size_t i = 0; i < size; ++i)
    {
      double min_dist = R_PosInf;
      for(size_t j = 0; j < peaks.size(); ++j)
      {

        dist_temp = peaks[j] == i ? 0:dist[locate(i, peaks[j], size)];
        if(min_dist > dist_temp)
        {
          min_dist = dist_temp;
          cluster[i] = j + 1;
        }
      }
    }

    return List::create(
      Named("peaks") = peaks,
      Named("rho") = rho,
      Named("delta") = delta,
      Named("dc") = dc,
      Named("method") = method,
      Named("clusters") = cluster,
      Named("halo") = NA_REAL
    );
  }
  else
  {
    vector<NodeDist> to_peak_dist(peaks.size());
    vector<double> halo_max_rho(peaks.size(), R_NegInf);//max rho of halo in each cluster
    int loc;

    for(size_t i = 0; i < size; ++i)
    {
      if(method == "neighbors") dc = 1/rho[i];
      for(size_t j = 0; j < peaks.size(); ++j)
      {
        loc = locate(i, peaks[j], size);
        to_peak_dist[j].dist = loc == -1 ? 0:dist[loc];
        to_peak_dist[j].node = j;
      }
      nth_element(to_peak_dist.begin(), to_peak_dist.begin() + 1, to_peak_dist.end(), Less);
      cluster[i] = to_peak_dist[0].node;
      if(to_peak_dist[1].dist - to_peak_dist[0].dist < dc)
        halo_max_rho[cluster[i]] = max(halo_max_rho[cluster[i]], rho[i]);
    }

    for(size_t i = 0; i < size; ++i)
    {
      if(rho[i] <= halo_max_rho[cluster[i]]) cluster[i] = -1;
      else cluster[i] ++;
    }

    return List::create(
      Named("peaks") = peaks,
      Named("rho") = rho,
      Named("delta") = delta,
      Named("dc") = dc,
      Named("method") = method,
      Named("clusters") = cluster,
      Named("halo") = halo_max_rho
    );
  }
}

NumericVector Dist(const NumericMatrix &data)
{
  size_t nrow = data.nrow();
  NumericVector dist((nrow - 1) * nrow / 2);
  parallelDist parallelDist(data, dist);
  parallelFor(1, nrow, parallelDist, 2e4);
  return(dist);
}

pair<size_t, int> get_grainSize(const size_t size, double percent, const int threads)
{
  int iter = 16;
  double lower = 0.1;
  double upper = 0.2;

  for(; iter > 1; iter/=2)
  {
    if(iter == 16 && iter*threads*percent < lower) break; //percent is too small use max iter(16)
    if(iter*threads*percent >= lower and iter*threads*percent <= upper) break;
  }

  if(iter == 1 and iter*threads*percent > upper)
    return(make_pair(0,0)); //data is too small for each worker to work

  size_t grainSize = size / (iter * threads) + 1;
  return(make_pair(grainSize, iter));
}

double getDc(const NumericVector &dist, double percent, const int threads)
{
  size_t k = percent * dist.size();
  NumericVector x = clone(dist);
  pair<size_t, int> grainSize;
  while(true)
  {
      if(grainSize.first == 0)
        break;
      NumericVector res(grainSize.second * threads * k);

      parallelgetDc parallelgetDc(x, grainSize.first, k, res);
      parallelFor(0, x.size(), parallelgetDc, grainSize.first);
      x = res;

    //for loop
    percent = 1.0/(grainSize.second * threads);
    grainSize = get_grainSize(x.size(), percent, threads);
  }

  nth_element(x.begin(), x.begin()+k, x.end(), less<double>());
  return(x[k-1]);
}

void getRho_withinDc(const NumericVector &dist, const size_t size, const double dc, NumericVector &rho)
{
  for(size_t i = 1; i < size; ++i)
  {
    for(size_t j = 0; j < i; ++j)
      if(dist[locate(i, j, size)] <= dc)
      {
        rho[i] ++;
        rho[j] ++;
      }
  }
}

void getRho_neighbors(const NumericVector &dist, const size_t size, const double percent, NumericVector &rho)
{
  size_t neighbors = size * percent;
  parallelRho parallelRho(dist, size, neighbors, rho);
  parallelFor(0, size, parallelRho, 2e4);
}

void getRho_gaussian(const NumericVector &dist, const size_t size, const double dc, NumericVector &rho)
{
  double temp_rho;
  for(size_t i = 0; i < size; ++i)
  {
    for(size_t j = i + 1; j < size; ++j)
    {

      temp_rho = exp(-pow(dist[locate(j, i, size)] / (2 * dc), 2) / 2);
      rho[i] += temp_rho;
      rho[j] += temp_rho;
    }
  }
}

void getDelta(const NumericVector &dist, const NumericVector &rho, size_t size, NumericVector &delta)
{
  double max_rho = *max_element(rho.begin(), rho.end());
  for(size_t i = 0; i < size; ++i)
  {
    for(size_t j = i + 1; j < size; ++j)
    {
      if(rho[i] == max_rho)
      {
        if(delta[i] == R_PosInf) delta[i] = dist[locate(j, i, size)];
        else delta[i] = max(delta[i], dist[locate(j, i, size)]);
        delta[j] = min(delta[j], dist[locate(j, i, size)]);
      }
      else if(rho[j] == max_rho)
      {
        if(delta[j] == R_PosInf) delta[j] = dist[locate(j, i, size)];
        else delta[j] = max(delta[j], dist[locate(j, i, size)]);
        delta[i] = min(delta[i], dist[locate(j, i, size)]);
      }
      else if(rho[i] > rho[j])
      {
        delta[j] = min(delta[j], dist[locate(j, i, size)]);
      }
      else
      {
        delta[i] = min(delta[i], dist[locate(j, i, size)]);
      }
    }
  }
}

