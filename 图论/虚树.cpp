ll fa[N], son[N], dep[N], siz[N], dfn[N], rnk[N], top[N];
ll dfscnt;
vector<ll> g[N];
ll mmin[N];

void dfs1(ll u, ll f, ll d) {
    son[u] = -1;
    siz[u] = 1;
    fa[u] = f;
    dep[u] = d;
    for (auto v:g[u]) {
        if (v == f) continue;
        dfs1(v, u, d + 1);
        siz[u] += siz[v];
        if (son[u] == -1 || siz[v] > siz[son[u]]) son[u] = v;
    }
}

void dfs2(ll u, ll t) {
    dfn[u] = ++dfscnt;
    rnk[dfscnt] = u;
    top[u] = t;
    if (son[u] == -1) return;
    dfs2(son[u], t);
    for (auto v:g[u]) {
        if (v == son[u] || v == fa[u]) continue;
        dfs2(v, v);
    }
}

ll lca(ll a, ll b) {
    while (top[a] != top[b]) {
        if (dep[top[a]] < dep[top[b]]) swap(a, b);
        a = fa[top[a]];
    }
    return dep[a] < dep[b] ? a : b;
}

struct edge {
    ll s, t, v;
};
edge e[N];

vector<int> vg[N];
int sta[N], tot;
int h[N];

void build(int *H, int num) {
    sort(H + 1, H + 1 + num, [](int a, int b) { return dfn[a] < dfn[b]; });
    sta[tot = 1] = 1, vg[1].clear();// 1 号节点入栈，清空 1 号节点对应的邻接表，设置邻接表边数为 1
    for (int i = 1, l; i <= num; ++i) {
        if (H[i] == 1) continue; //如果 1 号节点是关键节点就不要重复添加
        l = lca(H[i], sta[tot]); //计算当前节点与栈顶节点的 LCA
        if (l != sta[tot]) { //如果 LCA 和栈顶元素不同，则说明当前节点不再当前栈所存的链上
            while (dfn[l] < dfn[sta[tot - 1]]) {//当次大节点的 Dfs 序大于 LCA 的 Dfs 序
                vg[sta[tot - 1]].push_back(sta[tot]);
                vg[sta[tot]].push_back(sta[tot - 1]);
                tot--;
            } //把与当前节点所在的链不重合的链连接掉并且弹出
            if (dfn[l] > dfn[sta[tot - 1]]) { //如果 LCA 不等于次大节点（这里的大于其实和不等于没有区别）
                vg[l].clear();
                vg[l].push_back(sta[tot]);
                vg[sta[tot]].push_back(l);
                sta[tot] = l;//说明 LCA 是第一次入栈，清空其邻接表，连边后弹出栈顶元素，并将 LCA 入栈
            } else {
                vg[l].push_back(sta[tot]);
                vg[sta[tot]].push_back(l);
                tot--; //说明 LCA 就是次大节点，直接弹出栈顶元素
            }
        }
        vg[H[i]].clear();
        sta[++tot] = H[i];
        //当前节点必然是第一次入栈，清空邻接表并入栈
    }
    for (int i = 1; i < tot; ++i) {
        vg[sta[i]].push_back(sta[i + 1]);
        vg[sta[i + 1]].push_back(sta[i]);
    } //剩余的最后一条链连接一下
    return;
}
