#include <bits/stdc++.h>
using namespace std;
int n, m, root;
const int MAX_N = 500005;
const int MAX = 1 << 20;
vector<int> g[MAX_N];
vector<int> vs;
pair<int, int> tree[MAX * 2 + 10];
int fir[MAX_N];
int fa[MAX_N];
int dep[MAX_N];
void dfs(int k, int p, int d) {
    fa[k] = p;
    dep[k] = d;
    vs.push_back(k);
    for (int i = 0; i < g[k].size(); i++) {
        if (g[k][i] != p) {
            dfs(g[k][i], k, d + 1);
            vs.push_back(k);
        }
    }
}
void build(int k) {
    if (k >= MAX) return;
    build(k << 1);
    build(k << 1 | 1);
    tree[k] = min(tree[k << 1], tree[k << 1 | 1]);
}
pair<int, int> query(int k, int s, int e, int l, int r) {
    if (e < l || r < s) return pair<int, int>(INT_MAX, 0);
    if (l <= s && e <= r) return tree[k];
    return min(query(k << 1, s, (s + e) >> 1, l, r),
               query(k << 1 | 1, ((s + e) >> 1) + 1, e, l, r));
}
void init() {
    dfs(root, root, 0);
    for (int i = 0; i < MAX * 2 + 10; i++) tree[i] = pair<int, int>(INT_MAX, 0);
    for (int i = MAX; i < MAX + vs.size(); i++)
        tree[i] = pair<int, int>(dep[vs[i - MAX]], vs[i - MAX]);
    for (int i = 0; i < vs.size(); i++) {
        if (fir[vs[i]] == 0) fir[vs[i]] = i + 1;
    }
    build(1);
}
int lca(int a, int b) {
    return query(1, 1, MAX, min(fir[a], fir[b]), max(fir[a], fir[b])).second;
}
int main() {
    scanf("%d%d%d", &n, &m, &root);
    for (int i = 1; i < n; i++) {
        int a, b;
        scanf("%d%d", &a, &b);
        g[a].push_back(b);
        g[b].push_back(a);
    }
    init();
    for (int i = 1; i <= m; i++) {
        int a, b;
        scanf("%d%d", &a, &b);
        printf("%d\n", lca(a, b));
    }
}
