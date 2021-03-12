int pa[N];

void init(int n) {
    for (int i = 0; i <= n; i++) {
        pa[i] = i;
    }
}

int find(int a) {
    return pa[a] == a ? a : pa[a] = find(pa[a]);
}

int kruskal() {
    int kcnt = n;
    init(n);
    sort(e + 1, e + 1 + m, [](edge a, edge b) { return a.l < b.l; });
    for (int i = 1; i <= m; i++) {
        int u = find(e[i].from);
        int v = find(e[i].to);
        if (u == v) continue;
        w[++kcnt] = e[i].l;
        pa[kcnt] = pa[u] = pa[v] = kcnt;
        g[u].push_back(kcnt);
        g[v].push_back(kcnt);
        g[kcnt].push_back(u);
        g[kcnt].push_back(v);
    }
    return kcnt;
}
