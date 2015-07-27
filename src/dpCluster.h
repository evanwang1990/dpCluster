#ifndef dpCluster_h
#define dpCluster_h

#include <vector>
#include <Rcpp.h>
//[[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace std;
using namespace Rcpp;
using namespace RcppParallel;

struct NodeDist
{
  int node;
  double dist;
};

const bool Less(NodeDist & a, NodeDist & b)
{
  return(a.dist < b.dist);
}

const inline int locate(size_t i, size_t j, size_t n)
{
  if(i > j) return((2*n-j-1)*j/2+i-j-1);
  else if(i < j) return((2*n-i-1)*i/2-i+j-1);
  else return(-1);
}

struct parallelDist: public Worker
{
  const RMatrix<double> data;
  RVector<double> dist;
  parallelDist(NumericMatrix data, NumericVector dist)
    : data(data), dist(dist) {}

  void operator()(size_t begin, size_t end)
  {
    size_t nrow = data.nrow();
    for(size_t i = begin; i < end; ++i)
    {
      RMatrix<double>::Row node0 = data.row(i);
      for(size_t j = 0; j < i; ++j)
      {
        RMatrix<double>::Row node1 = data.row(j);
        double dist_ = 0;
        for(size_t len = 0; len < node1.length(); ++len)
        {
          dist_ += pow(node1[len] - node0[len], 2);
        }
        dist[locate(i, j, nrow)] = sqrt(dist_);
      }
    }
  }
};

struct parallelgetDc: public Worker
{
  RVector<double> dist;
  const size_t grainSize;
  const size_t k;
  RVector<double> kmins;

  parallelgetDc(NumericVector dist, const size_t grainSize, const size_t k, NumericVector kmins)
    : dist(dist), grainSize(grainSize), k(k), kmins(kmins) {}

  void operator()(size_t begin, size_t end)
  {
    int iter = floor(double(begin)/ (grainSize - 1));
    nth_element(dist.begin()+begin, dist.begin()+begin+k, dist.begin()+end, less<double>());
    copy(dist.begin() + begin, dist.begin()+begin+k, kmins.begin()+iter*k);
  }
};

struct parallelRho: public Worker
{
  const RVector<double> dist;
  const int size;
  const int neighbors;
  RVector<double> rho;

  parallelRho(const NumericVector dist, const int size, const int neighbors, NumericVector rho)
    : dist(dist), size(size), neighbors(neighbors), rho(rho) {}

  void operator()(size_t begin, size_t end)
  {
    for(size_t i = begin; i < end; ++i)
    {
      vector<double> dist_(size - 1, 0);
      for(size_t j = 0; j < size; ++j)
      {
        if(i != j) dist_[j] = dist[locate(i, j, size)];
      }

      nth_element(dist_.begin(), dist_.begin() + neighbors + 1, dist_.end(), less<double>());
      double rho_ = 0;
      for(vector<double>::iterator pt = dist_.begin(); pt != dist_.begin() + neighbors + 1; pt++)
      {
        rho_ += *pt;
      }
      rho[i] = neighbors / rho_;
    }
  }
};

#endif
