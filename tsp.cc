#include <algorithm>
#include <cfloat>
#include <iostream>
#include <stack>
#include <vector>


using namespace std;


struct Edge {
  int head, tail;
  float cost;

  Edge(int h, int t, float c) : head(h), tail(t), cost(c) { }

  bool operator<(const Edge& e) const
  {
    return cost < e.cost;
  }
};

struct UnionFind {
  vector<int> parent;
  vector<int> rank;

  UnionFind(int n) : parent(n), rank(n)
  {
    for (int i = 0; i < n; ++i)
      make_set(i);
  }

  void make_set(int i)
  {
    parent[i] = i;
    rank[i] = 0;
  }

  inline int find_set(int i)
  {
    if (parent[i] != i)
      parent[i] = find_set(parent[i]);
    return parent[i];
  }

  void join(int i, int j)
  {
    i = find_set(i);
    j = find_set(j);
    if (rank[i] > rank[j])
      parent[j] = i;
    else {
      parent[i] = j;
      if (rank[i] == rank[j])
        rank[j] += 1;
    }
  }
};


// Kruskal's minimum spanning tree algorithm
// expects symmetric costs (i.e. an undirected graph)
vector<Edge> mst(int n, vector<vector<float> >& cost)
{
  vector<Edge> edges;

  for (int i = 0; i < n; ++i)
    for (int j = i + 1; j < n; ++j)
      edges.push_back(Edge(i, j, cost[i][j]));

  sort(edges.begin(), edges.end());

  vector<Edge> tree;
  UnionFind forest(n);

  int size = 0;
  for (int i = 0; size < n - 1; ++i) {
    int u = edges[i].head;
    int v = edges[i].tail;
    if (forest.find_set(u) != forest.find_set(v)) {
      tree.push_back(edges[i]);
      forest.join(u, v);
      ++size;
    }
  }

  return tree;
}

// 2-approximation for metric travelling salesman problem
// expects symmetric costs that respect the triangle inequality
vector<int> metric_tsp(int n, vector< vector<float> >& cost)
{
  vector<Edge> tree = mst(n, cost);

  // adjacency lists relative to the tree subgraph
  vector<int> adjacency_lists[n];

  for (int i = 0; i < n - 1; ++i) {
    int u = tree[i].head;
    int v = tree[i].tail;
    adjacency_lists[u].push_back(v);
    adjacency_lists[v].push_back(u);
  }

  vector<int> tour;

  // traverse the Eulerian circuit implicitly defined by the tree, skipping
  // already visited vertices
  bool visited[n];
  fill(visited, visited + n, false);

  stack<int> open;
  open.push(0);

  while (tour.size() < n) {
    int u = open.top();
    open.pop();
    if (!visited[u]) {
      tour.push_back(u);
      visited[u] = true;
      for (vector<int>::iterator v = adjacency_lists[u].begin(); 
           v != adjacency_lists[u].end();
           ++v)
        open.push(*v);
    }
  }
  tour.push_back(0);

  return tour;
}

// exact exponential-time algorithm using dynamic programming
vector<int> tsp(int n, vector< vector<float> >& cost)
{
  long nsub = 1 << n;
  vector< vector<float> > opt(nsub, vector<float>(n));

  for (long s = 1; s < nsub; s += 2)
    for (int i = 1; i < n; ++i) {
      vector<int> subset;
      for (int u = 0; u < n; ++u)
        if (s & (1 << u))
          subset.push_back(u);

      if (subset.size() == 2)
        opt[s][i] = cost[0][i];

      else if (subset.size() > 2) {
        float min_subpath = FLT_MAX;
        long t = s & ~(1 << i);
        for (vector<int>::iterator j = subset.begin(); j != subset.end(); ++j)
          if (*j != i && opt[t][*j] + cost[*j][i] < min_subpath)
            min_subpath = opt[t][*j] + cost[*j][i];
        opt[s][i] = min_subpath;
      }
    }
  
  vector<int> tour;
  tour.push_back(0);

  bool selected[n];
  fill(selected, selected + n, false);
  selected[0] = true;

  long s = nsub - 1;

  for (int i = 0; i < n - 1; ++i) {
    int j = tour.back();
    float min_subpath = FLT_MAX;
    int best_k;
    for (int k = 0; k < n; ++k)
      if (!selected[k] && opt[s][k] + cost[k][j] < min_subpath) {
        min_subpath = opt[s][k] + cost[k][j];
        best_k = k;
      }
    tour.push_back(best_k);
    selected[best_k] = true;
    s -= 1 << best_k;
  }
  tour.push_back(0);

  return tour;
}



int main()
{
  int n = 4;

  vector< vector<float> > c(n, vector<float>(n));

  c[0][1] = c[1][0] = 1;
  c[0][2] = c[2][0] = 3;
  c[0][3] = c[3][0] = 2;
  c[1][2] = c[2][1] = 2;
  c[1][3] = c[3][1] = 4;
  c[2][3] = c[3][2] = 3;

  cout << endl << "Metric TSP approximation:" << endl << "Path: ";
  vector<int> tour = metric_tsp(n, c);
  for (int i = 0; i <= n; i++)
    cout << tour[i % n] << " ";
  cout << endl;
  cout << "Cost: ";
  float total_cost = 0;
  for (int i = 1; i <= n; ++i)
    total_cost += c[tour[i - 1]][tour[i % n]];
  cout << total_cost << endl << endl;

  cout << "Exact exponential solution:" << endl << "Path: ";
  tour = tsp(n, c);
  for (int i = 0; i <= n; i++)
    cout << tour[i % n] << " ";
  cout << endl;
  cout << "Cost: ";
  total_cost = 0;
  for (int i = 1; i <= n; ++i)
    total_cost += c[tour[i - 1]][tour[i % n]];
  cout << total_cost << endl << endl;

  return 0;
}



