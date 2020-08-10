#include <bits/stdc++.h>
using namespace std;
struct SCC {
    static const int MAXV = 100000;
    int V;
    vector<int> g[MAXV], rg[MAXV], vs;
    bool used[MAXV];
    int cmp[MAXV];

    void add_edge(int from, int to) {
        g[from].push_back(to);
        rg[to].push_back(from);
    }

    void dfs(int v) {
        used[v] = 1;
        for (int i = 0; i < g[v].size(); i++) {
            if (!used[g[v][i]]) dfs(g[v][i]);
        }
        vs.push_back(v);
    }

    void rdfs(int v, int k) {
        used[v] = 1;
        cmp[v] = k;
        for (int i = 0; i < rg[v].size(); i++) {
            if (!used[rg[v][i]]) rdfs(rg[v][i], k);
        }
    }

    int solve() {
        memset(used, 0, sizeof(used));
        vs.clear();
        for (int v = 1; v <= V; v++) {
            if (!used[v]) dfs(v);
        }
        memset(used, 0, sizeof(used));
        int k = 0;
        for (int i = (int)vs.size() - 1; i >= 0; i--) {
            if (!used[vs[i]]) rdfs(vs[i], ++k);
        }
        return k;
    }

    void init(int n) {
        V = n;
        vs.clear();
        for (int i = 0; i < MAXV; i++) {
            g[i].clear();
            rg[i].clear();
            used[i] = 0;
            cmp[i] = 0;
        }
    }

} scc;

//记得调用init()
