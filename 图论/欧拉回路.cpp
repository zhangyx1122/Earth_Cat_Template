#include <bits/stdc++.h>

using namespace std;
typedef long long ll;
const int N = 1e6 + 10;


int stk[N], top;
struct edge {
    int to, idx;
};

vector<edge> g[N];

namespace Euler1 {  //无向图欧拉回路
    bool vis[N];
    int cur[N];

    void dfs(int u, const int &w) {
        vis[abs(w)] = true;
        for (int &i = cur[u]; i < g[u].size();) {
            int idx = g[u][i].idx, v = g[u][i].to;
            i++;
            if (!vis[abs(idx)]) dfs(v, idx);
        }
        stk[++top] = w;
    }

    bool solve(int n) {
        // init();
        for (int i = 0; i <= n; i++) cur[i] = 0;
        for (int i = 0; i <= n; i++) vis[i] = 0;
        // calculate degree
        for (int i = 1; i <= n; i++) {
            if (g[i].size() & 1) return false;
        }
        // Hierholzer
        for (int i = 1; i <= n; i++)
            if (!g[i].empty()) {
                dfs(i, 0);
                break;
            }
        return true;
    }
}  // namespace Euler1

namespace Euler2 {  // 有向图欧拉回路
    int deg[N], cur[N];

    void dfs(int u, const int &w) {
        for (int &i = cur[u]; i < g[u].size();) {
            int idx = g[u][i].idx, v = g[u][i].to;
            i++;
            dfs(v, idx);
        }
        stk[++top] = w;
    }

    bool solve(int n) {
        // init
        for (int i = 0; i <= n; i++) deg[i] = 0;
        for (int i = 0; i <= n; i++) cur[i] = 0;
        // calculate degree
        for (int i = 1; i <= n; ++i) {
            for (auto x: g[i]) deg[i]++, deg[x.to]--;
        }
        for (int i = 1; i <= n; ++i)
            if (deg[i]) return false;
        // Hierholzer
        for (int i = 1; i <= n; ++i)
            if (!g[i].empty()) {
                dfs(i, 0);
                break;
            }
        return true;
    }
}  // namespace Euler2

int main() {
    int t, n, m;
    cin >> t >> n >> m;
    for (int u, v, i = 1; i <= m; i++) {
        cin >> u >> v;
        g[u].push_back({v, i});
        if (t == 1) g[v].push_back({u, -i});
    }
    // solve
    bool flag = t == 1 ? Euler1::solve(n) : Euler2::solve(n);
    // output
    if (!flag || (m > 0 && top - 1 < m))
        puts("NO");
    else {
        puts("YES");
        for (int i = top - 1; i > 0; --i) printf("%d%c", stk[i], " \n"[i == 1]);
    }
    return 0;
}
