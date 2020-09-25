#include <bits/stdc++.h>
using namespace std;
const int MAXN = 100000;

int c[MAXN];
int topo[MAXN], t, V;
vector<int> g[MAXN];

bool dfs(int u) {
    c[u] = -1;
    for (int i = 0; i < g[u].size(); i++) {
        int v = g[u][i];
        if (c[v] < 0)
            return false;
        else if (!c[v] && !dfs(v))
            return false;
    }
    c[u] = 1;
    topo[t--] = u;
    return true;
}

bool toposort(int n) {
    V = n;
    t = n;
    memset(c, 0, sizeof(c));
    for (int u = 1; u <= V; u++)
        if (!c[u] && !dfs(u)) return false;
    return true;
}
