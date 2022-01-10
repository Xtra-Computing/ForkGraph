# ForkGraph: Cache-Efficient Fork-Processing Patterns on Large Graphs


Organization
--------

This repository contains code for our SIGMOD paper "Cache-Efficient Fork-Processing Patterns on Large Graphs". It is implemented based on Julian Shun's **Ligra** https://github.com/jshun/ligra


Compilation
--------

Compiler:
* g++ >= 7.5.0


Build system:
* CMake >= 3.12


To build:
```sh
$ cd ForkGraph/ # go to the source folder
$ mkdir build
$ cd build && cmake ..
```


Input Formats
-----------
We support the adjacency graph format used by the [Problem Based Benchmark suite](http://www.cs.cmu.edu/~pbbs/benchmarks/graphIO.html) and [Ligra](https://github.com/jshun/ligra).

The adjacency graph format starts with a sequence of offsets one for each vertex, followed by a sequence of directed edges ordered by their source vertex. The offset for a vertex i refers to the location of the start of a contiguous block of out edges for vertex i in the sequence of edges. The block continues until the offset of the next vertex, or the end if i is the last vertex. All vertices and offsets are 0 based and represented in decimal. The specific format is as follows:

```
AdjacencyGraph
<n>
<m>
<o0>
<o1>
...
<o(n-1)>
<e0>
<e1>
...
<e(m-1)>
```

This file is represented as plain text.

Weighted graphs are represented in the weighted adjacency graph format. The file should start with the string "WeightedAdjacencyGraph". The m edge weights should be stored after all of the edge targets in the .adj file.

**Using SNAP graphs**

Graphs from the [SNAP dataset collection](https://snap.stanford.edu/data/index.html) are commonly used for graph algorithm benchmarks. Pleae use the tool that converts the most common SNAP graph format to the adjacency graph format that ForkGraph accepts. The tool can be found in [GBBS: Graph Based Benchmark Suite](https://github.com/ParAlg/gbbs).

## How to cite ForkGraph
If you use ForkGraph in your paper, please cite our work.

```
@inproceedings{lu2021cache,
  title={Cache-Efficient Fork-Processing Patterns on Large Graphs},
  author={Lu, Shengliang and Sun, Shixuan and Paul, Johns and Li, Yuchen and He, Bingsheng},
  booktitle={Proceedings of the 2021 International Conference on Management of Data},
  pages={1208--1221},
  year={2021}
}
```
