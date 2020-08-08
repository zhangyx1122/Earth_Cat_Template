#include <bits/stdc++.h>
using namespace std;

struct SCC {
    static const int MAXN = 100000;
    vector<int> g[MAXN];
    int dfn[MAXN], lowlink[MAXN], sccno[MAXN], dfs_clock, scc_cnt;
    stack<int> S;

    void dfs(int u) {
        dfn[u] = lowlink[u] = ++dfs_clock;
        S.push(u);
        for (int i = 0; i < g[u].size(); i++) {
            int v = g[u][i];
            if (!dfn[v]) {
                dfs(v);
                lowlink[u] = min(lowlink[u], lowlink[v]);
            } else if (!sccno[v]) {
                lowlink[u] = min(lowlink[u], dfn[v]);
            }
        }
        if (lowlink[u] == dfn[u]) {
            ++scc_cnt;
            for (;;) {
                int x = S.top();
                S.pop();
                sccno[x] = scc_cnt;
                if (x == u) break;
            }
        }
    }

    void solve(int n) {
        dfs_clock = scc_cnt = 0;
        memset(sccno, 0, sizeof(sccno));
        memset(dfn, 0, sizeof(dfn));
        memset(lowlink, 0, sizeof(lowlink));
        for (int i = 1; i <= n; i++) {
            if (!dfn[i]) dfs(i);
        }
    }
} scc;

// scc_cnt为SCC计数器，sccno[i]为i所在SCC的编号
// vector<int> g[MAXN]中加边
//之后再补充init()
