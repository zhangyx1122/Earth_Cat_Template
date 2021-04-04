#include <bits/stdc++.h>
using namespace std;
const int MAXN = 10005;
const int INF = 1000000000;
struct edge {
    int to, length;
    edge() {}
    edge(int a, int b) : to(a), length(b) {}
};


vector<edge> g[MAXN];

bool centroid[MAXN];
int subtree_size[MAXN];

int ans;

//计算子树大小
int compute_subtree_size(int v, int p) {
    int c = 1;
    for (int i = 0; i < g[v].size(); i++) {
        int w = g[v][i].to;
        if (w == p || centroid[w]) continue;
        c += compute_subtree_size(w, v);
    }
    subtree_size[v] = c;
    return c;
}

//查找重心，t为连通分量大小
// pair（最大子树顶点数，顶点编号）
pair<int, int> search_centroid(int v, int p, int t) {
    pair<int, int> res = pair<int, int>(INF, -1);
    int s = 1, m = 0;
    for (int i = 0; i < g[v].size(); i++) {
        int w = g[v][i].to;
        if (w == p || centroid[w]) continue;
        res = min(res, search_centroid(w, v, t));
        m = max(m, subtree_size[w]);
        s += subtree_size[w];
    }
    m = max(m, t - s);
    res = min(res, pair<int, int>(m, v));
    return res;
}

void init(int n) {
    memset(centroid, 0, sizeof(centroid));
    memset(subtree_size, 0, sizeof(subtree_size));
    for (int i = 0; i <= n; i++) g[i].clear();
    ans = 0;
}

int solve(int u) {
    compute_subtree_size(u, -1);
    int s = search_centroid(u, -1, subtree_size[u]).second;
    centroid[s] = 1;
    for (int i = 0; i < g[s].size(); i++) {
        int v = g[s][i].to;
        if (centroid[v]) continue;
        /*solve()*/
    }
    /*do something*/
    centroid[s] = 0;
    return ans;
}
