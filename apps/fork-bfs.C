// This code is part of the project "Ligra: A Lightweight Graph Processing
// Framework for Shared Memory", presented at Principles and Practice of 
// Parallel Programming, 2013.
// Copyright (c) 2013 Julian Shun and Guy Blelloch
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#define WEIGHTED 1
#include "ligra.h"
#include <vector>
#include <queue>
#include "time.hpp"

//#define STATS
//#define COUNTSINGLE
//#define COUNTALL
//#define WORK

#ifdef WORK
long _N, _M;
#endif

#define YIELD1 15
#define YIELD2 3

// no need to define prioirty functor. can use the default one
typedef std::pair<intE, uintE> WN;
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

template<class vertex>
void Compute(graph<vertex> &GA, graph<vertex> &GB, int NP) {
  int yield_h1 = YIELD1;
  int yield_h2 = YIELD2;

#ifdef COUNTSINGLE
  int totalVisitedPartitionCount[start_l.size()];
  int totalValidVisitedPartitionCount[start_l.size()];
  for (int _i = 0; _i < start_l.size(); _i++) {
    totalVisitedPartitionCount[_i] = 0;
    totalValidVisitedPartitionCount[_i] = 0;
  }
#endif
#ifdef COUNTALL
  int _totalVisitedPartitionCount = 0;
#endif

  long n = GA.n;
  std::vector<long> s_l = start_l;
  long k = s_l.size();
  int count[NP];
  int closest[NP];
  std::priority_queue<WN, std::vector<WN>, std::greater<WN>> lq_l[k][NP];

#ifdef STATS
  long touch_vertex_c = 0;
  long touch_edge_c = 0;
#endif
  {
    parallel_for (int i = 0; i < NP; i++) {
      count[i] = 0;
      closest[i] = INT_MAX;
    }
  }

  //initialize BfsLevel to "infinity"
  std::vector<intE *> BfsLevel_l;
  std::vector<int *> Visited_l;
  std::vector<long> fq_l[k];
  std::vector<long> buffer_buckets_l[k][NP];
  for (int i = 0; i < k; i++) {
    { parallel_for (int j = 0; j < NP; j++) buffer_buckets_l[i][j].reserve(100000); }
  }
  for (int i = 0; i < k; i++) {
    intE *BfsLevel = newA(intE, n);
    int *Visited = newA(int, n);
    { parallel_for (long j = 0; j < n; j++) BfsLevel[j] = INT_MAX / 2; }
    { parallel_for (long j = 0; j < n; j++) Visited[j] = 0; }
    fq_l[i].reserve(GA.n);

    BfsLevel[s_l[i]] = 0;
    BfsLevel_l.push_back(BfsLevel);
    Visited_l.push_back(Visited);
    int pid = id2part[s_l[i]];
    buffer_buckets_l[i][pid].push_back(s_l[i]);
    count[pid]++;
    closest[pid] = 0;
    part2process.push(std::make_pair(0, pid));
  }

  Timer t;
  t.tic();
  int pi = -1;
  while (pi = schedule_partition(), pi > 0) {
#ifdef COUNTALL
    _totalVisitedPartitionCount+=1;
#endif
    count[pi] = 0;
    closest[pi] = INT_MAX;
    auto f = [&](int qi,
                 std::vector<long> *buffer_buckets,
                 std::priority_queue<WN, std::vector<WN>, std::greater<WN>> *local_frontier) {
      auto &BfsLevel = BfsLevel_l[qi];
      auto &Visited = Visited_l[qi];
      auto &bucket = buffer_buckets[pi];
      int touched_vertices = 0;
      {
        auto &pq = local_frontier[pi];
        std::vector<long> local_buffer;
        local_buffer.reserve(GA.n);
        // get local frontiers from buffer
        for (int i = 0; i < bucket.size(); i++) {
          pq.push(std::make_pair(BfsLevel[bucket[i]], bucket[i]));
          Visited[bucket[i]] = 1;
        }
        int counter = pq.size() / yield_h1;
        int current_level = pq.top().first;
        while (!pq.empty()) {
          uintE vi = pq.top().second;
          intE td = pq.top().first;
          if (counter-- <= 0) {
            if (part2process.top().first > td) {
              counter = pq.size() / yield_h1 * 2;
              continue;
            }
            closest[pi] = td;
#pragma omp critical
            part2process.push(std::make_pair(td, pi));
            break;
          }
          if (td > current_level + yield_h2) {
            closest[pi] = td;
#pragma omp critical
            part2process.push(std::make_pair(td, pi));
            break;
          }
          pq.pop();
          //sent from other partitions
          if (td == BfsLevel[vi]) {
#ifdef WORK
#pragma omp critical
            _N++;
#endif
            touched_vertices++;
            local_buffer.push_back(vi);
#ifdef STATS
            writeAdd(&touch_vertex_c, (long) 1);
#endif
            vertex &v = GA.V[vi];
            Visited[vi] = 1;
            auto degree = v.getOutDegree();
#ifdef STATS
            writeAdd(&touch_edge_c, (long) degree);
#endif
#ifdef WORK
#pragma omp critical
            _M += degree;
#endif
            for (int d = 0; d < degree; d++) {
              uintE ngh = v.getOutNeighbor(d);
              if (BfsLevel[ngh] > BfsLevel[vi] + 1) {
                BfsLevel[ngh] = BfsLevel[vi] + 1;
                pq.push(std::make_pair(BfsLevel[ngh], ngh));
              }
            }
          }
        }
        for (int i = 0; i < local_buffer.size(); i++) {
          auto vi = local_buffer[i];
          if (!Visited[vi]) {
            continue;
          }
#ifdef STATS
          writeAdd(&touch_vertex_c, (long) 1);
#endif
          Visited[vi] = 0;
          vertex &v = GB.V[vi];
          auto degree = v.getOutDegree();
#ifdef STATS
          writeAdd(&touch_edge_c, (long) degree);
#endif
#ifdef WORK
#pragma omp critical
          _M += degree;
#endif
          for (int d = 0; d < degree; d++) {
            uintE ngb = v.getOutNeighbor(d);
            if (BfsLevel[ngb] > BfsLevel[vi] + 1) {
              BfsLevel[ngb] = BfsLevel[vi] + 1;
              int pid = id2part[ngb];
              buffer_buckets[pid].push_back((long) ngb);
              if (BfsLevel[ngb] < closest[pid]) {
                closest[pid] = BfsLevel[ngb];
#pragma omp critical
                part2process.push(std::make_pair(BfsLevel[ngb], pid));
              }
              count[pid]++;
            }
          }
        }
#ifdef COUNTSINGLE
        if (touched_vertices * 5 > part_left[pi] - part_right[pi]) {
          totalValidVisitedPartitionCount[qi]++;
        }
#endif
        local_buffer.clear();
      }
      bucket.clear();
    };

#pragma omp parallel for schedule(dynamic)
    for (int qi = 0; qi < k; qi++) {
      if (!buffer_buckets_l[qi][pi].empty() || !lq_l[qi][pi].empty()) {
#ifdef COUNTSINGLE
        totalVisitedPartitionCount[qi]++;
#endif
        // sequential algorithm
        f(qi, buffer_buckets_l[qi], lq_l[qi]);
      }
    }
  }
  t.toc();
  t.print_ms("exuection time");

#ifdef STATS
  std::cout << "Touched Vertices: " << touch_vertex_c << std::endl;
  std::cout << "Touched Edges: " << touch_edge_c << std::endl;
#endif

#ifdef COUNTSINGLE
  for (int i = 0; i < start_l.size(); i++) {
    std::cout << totalVisitedPartitionCount[i] << " ";
  }
  std::cout << std::endl;
  for (int i = 0; i < start_l.size(); i++) {
    std::cout << totalValidVisitedPartitionCount[i] << " ";
  }
  std::cout << std::endl;
#endif

#ifdef COUNTALL
  std::cout << _totalVisitedPartitionCount << std::endl;
#endif
#ifdef WORK
  std::cout << _N << " " << _M << std::endl;
#endif
}
