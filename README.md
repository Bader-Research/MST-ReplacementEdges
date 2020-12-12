# MST-ReplacementEdges: Find Minimum Spanning Tree Replacement Edges

Given an undirected, weighted graph, the minimum spanning tree (MST)
is a tree that connects all of the vertices of the graph with minimum
sum of edge weights. In real world applications, network designers
often seek to quickly find a replacement edge for each edge in the
MST. This code finds the lowest cost replacement edge for each edge of
the MST, based on this paper:

David A. Bader and Paul Burkhardt, "[A Linear Time Algorithm for
Finding Minimum Spanning Tree Replacement
Edges](https://arxiv.org/abs/1908.03473)", arXiv:1908.03473v3, 2020.

- MST-ReplacementEdges-Tarjan: This subdirectory contains a superlinear
  implementation using the Tarjan approach for disjoint set unions.

- MST-ReplacementEdges-Gabow-Tarjan: This subdirectory implements the
  linear-time Bader-Burkhardt algorithm using the Gabow-Tarjan
  approach for disjoint set unions when the union tree is known in
  advance.

This code can also find the most vital edge -- the tree edge whose
removal causes the highest cost -- in linear time.

