// This code is part of the project "Parallel Local Clustering
// Algorithms".  Copyright (c) 2016 Julian Shun
//
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights (to use, copy,
// modify, merge, publish, distribute, sublicense, and/or sell copies
// of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
// BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
// ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
// CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

// This is a serial implementation of PR-Nibble using the optimized
// update rule and uses a priority queue.  Currently only works with
// uncompressed graphs, and not with compressed graphs.
#include "ligra.h"
#include <unordered_map>
#include <set>
#include <queue>
#include "sweep.h"
using namespace std;

typedef pair<double, double> pairDouble;
typedef pair<uintE, double> pairIF;
struct pq_compare {
  bool operator()(pairIF a, pairIF b) {
    return a.second > b.second;
  }
};

typedef std::pair<double, uintE> WN;
std::priority_queue<WN, std::vector<WN>, std::greater<WN>> part2process;
int schedule_partition() {
  int pi = -1;
  if (!part2process.empty()) {
    pi = part2process.top().second;
    while (!part2process.empty() && pi == part2process.top().second) {
      part2process.pop();
    }
  }
  return pi;
}


//#define WORK
#ifdef WORK
long _N, _M;
#endif

#define YIELD1 15

timer t1;
template<class vertex>
void Compute(graph<vertex> &GA, graph<vertex> &GB, int NP) {
  int yield_h1 = YIELD1;
  long k = start_l.size();
//  const double alpha = P.getOptionDoubleValue("-a", 0.15);
  const double alpha = 0.15;
//  const double epsilon = P.getOptionDoubleValue("-e", 0.000000001);
  const double epsilon = 0.000000001;
  const intE n = GA.n;
  const double twoAlphaOverOnePlusAlpha = 2 * alpha / (1 + alpha);
  const double oneMinusAlphaOverOnePlusAlpha = (1 - alpha) / (1 + alpha);

  int count[NP];
  double closest[NP];
  {
    parallel_for (int i = 0; i < NP; i++) {
      count[i] = 0;
      closest[i] = 0;
    }
  }

  unordered_map<uintE, pairDouble> pr_l[k];
  multiset<pairIF, pq_compare> q_l[k];
  std::vector<pairIF> buffer_bucket_l[k][NP];
  for (int i = 0; i < k; i++) {
    for (int j = 0; j < NP; j++) {
      buffer_bucket_l[i][j].reserve(GA.n / 100);
    }
    int pid = id2part[start_l[i]];
    pr_l[i][start_l[i]] = make_pair(0.0, 1.0); //starting vertex
    buffer_bucket_l[i][pid].emplace_back(start_l[i], 1.0);
    count[pid]++;
    closest[pid] = 1.0;
  }

  t1.start();
  int pi = -1;
  while (pi = schedule_partition(), pi > 0) {
    count[pi] = 0;
    closest[pi] = 0;
    int nprocs = omp_get_max_threads();
    int nthreads1 = 1;
    auto f = [&](int qi,
                 std::vector<pairIF> *buffer_buckets) {
      auto &Frontier = buffer_buckets[pi];
      auto &pr = pr_l[qi];
      auto q = q_l[qi];
      int nthreads1 = 1;
#pragma omp parallel num_threads(nthreads1)
      {
        std::vector<pairIF> local_buffer;
#pragma omp for nowait schedule(dynamic)
        for (int i = 0; i < Frontier.size(); i++) {
          q.insert(Frontier[i]);
        }
        long totalPushes = 0;
        while (!q.empty()) {
          totalPushes++;
          if (totalPushes > yield_h1) {
            break;
          }
          pairIF top = *q.begin();
          uintE v = top.first;
          q.erase(q.begin());
          pairDouble v_pr = pr[v];
          uintE d = GA.V[v].getOutDegree();
          v_pr.first += twoAlphaOverOnePlusAlpha * v_pr.second;
          double old_vr = v_pr.second;
          local_buffer.emplace_back(make_pair(v, old_vr));
          v_pr.second = 0;
          pr[v] = v_pr;
#ifdef WORK
          #pragma omp critical
        _M+=d;
#endif
          for (long i = 0; i < d; i++) {
            uintE ngh = GA.V[v].getOutNeighbor(i);
            pairDouble ngh_pr = pr[ngh]; //creates default entry if non-existent in table
            double oldRes = ngh_pr.second;
            uintE ngh_d = GA.V[ngh].getOutDegree() + GB.V[ngh].getOutDegree();
            ngh_pr.second += old_vr * oneMinusAlphaOverOnePlusAlpha / (d + GB.V[v].getOutDegree());
            pr[ngh] = ngh_pr;
            double epsd = epsilon * ngh_d;
            if (ngh_pr.second > epsd && oldRes < epsd) {
              //if previous residual is small, that means it is not in q, so insert it
              q.insert(make_pair(ngh, ngh_pr.second / ngh_d));
            }
          }
        }
#ifdef WORK
        _N += totalPushes;
#endif
        Frontier.clear();

        for (auto p : local_buffer) {
          auto v = p.first;
          double old_vr = p.second;
          uintE d = GB.V[v].getOutDegree();
          for (long i = 0; i < d; i++) {
            uintE ngh = GB.V[v].getOutNeighbor(i);
            int pid = id2part[ngh];
            pairDouble ngh_pr;
            ngh_pr = pr[ngh]; //creates default entry if non-existent in table
            double oldRes = ngh_pr.second;
            ngh_pr.second += old_vr * oneMinusAlphaOverOnePlusAlpha / (d + GA.V[v].getOutDegree());
            pr[ngh] = ngh_pr;
            uintE ngh_d = GB.V[ngh].getOutDegree() + GA.V[ngh].getOutDegree();
#ifdef WORK
            #pragma omp critical
          _M+=GB.V[ngh].getOutDegree();
#endif
            double epsd = epsilon * ngh_d;
            if (ngh_pr.second > epsd && oldRes < epsd) {
              //if previous residual is small, that means it is not in q, so insert it
              auto pr_value = ngh_pr.second / ngh_d;
              if (pr_value < epsd) {
                continue;
              }
              buffer_buckets[pid].emplace_back(ngh, pr_value);
              part2process.push(std::make_pair(pr_value, pid));
              count[pid]++;
            }
          }
        }
      }
    };

#pragma omp parallel for num_threads(nprocs) schedule(dynamic)
    for (int qi = 0; qi < k; qi++) {
      if (!buffer_bucket_l[qi][pi].empty()) {
        f(qi, buffer_bucket_l[qi]);
      }
    }
  }

  t1.stop();
  t1.reportTotal("computation time");
#ifdef WORK
  std::cout << _N << " " << _M << std::endl;
#endif
}
