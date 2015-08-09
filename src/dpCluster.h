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

struct parallelSNN: public Worker
{
  const RVector<double> dist;
  const int size;
  const int k;
  RMatrix<int> k_neighbors;

  parallelSNN(const NumericVector dist, const int size, const int k, IntegerMatrix k_neighbors)
    : dist(dist), size(size), k(k), k_neighbors(k_neighbors) {};

  void operator()(size_t begin, size_t end)
  {
    vector<NodeDist> dist_(size - 1);
    for(size_t i = begin; i < end; ++i)
    {
      vector<NodeDist>::iterator pt = dist_.begin();
      for(size_t j = 0; j < size; ++j)
      {
        if(i != j)
        {
          pt->node = j;
          pt->dist = dist[locate(i, j, size)];
          ++pt;
        }
      }

      nth_element(dist_.begin(), dist_.begin() + k - 1, dist_.end(), Less);

      for(size_t n = 0; n < k; ++n)
      {
        k_neighbors(i, n) = dist_[n].node;
      }
    }
  }
};

struct parallelgetDc: public Worker
{
  RVector<double> dist;
  const size_t grainSize;
  const size_t k;
  const bool Less;
  RVector<double> kmins;

  parallelgetDc(NumericVector dist, const size_t grainSize, const size_t k, const bool Less, NumericVector kmins)
    : dist(dist), grainSize(grainSize), k(k), Less(Less), kmins(kmins) {}

  void operator()(size_t begin, size_t end)
  {
    int iter = floor(double(begin)/ (grainSize - 1));
    if(Less)
    {
      nth_element(dist.begin()+begin, dist.begin()+begin+k, dist.begin()+end, less<double>());
    }
    else
    {
      nth_element(dist.begin()+begin, dist.begin()+begin+k, dist.begin()+end, greater<double>());
    }

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
