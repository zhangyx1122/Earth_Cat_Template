#include <bits/stdc++.h>

using namespace std;
typedef long long ll;
const int N = 2e5 + 10;

int vis[N], now;

vector<int> g[N];
int fa[N], son[N], siz[N], ans[N];

void insert(int pos) {
    vis[pos] = 1;
    now = now + 1 - vis[pos - 1] - vis[pos + 1];
}

void remove(int pos) {
    vis[pos] = 0;
    now = now - 1 + vis[pos - 1] + vis[pos + 1];
}

void dfs1(ll u, ll f) {
    siz[u] = 1;
    fa[u] = f;
    son[u] = -1;
    for (auto v:g[u]) {
        if (v == f) continue;
        dfs1(v, u);
        siz[u] += siz[v];
        if (son[u] == -1 || siz[v] > siz[son[u]]) son[u] = v;
    }
}

void add(int u, int exc, int op) {
    if (op) insert(u);
    else remove(u);
    for (auto x:g[u]) {
        if (x == fa[u] || x == exc) continue;
        add(x, exc, op);
    }
}

void dfs(ll u, ll opt) {
    for (auto x:g[u]) {
        if (x == fa[u] || x == son[u]) continue;
        dfs(x, 0);
    }
    if (son[u] != -1) dfs(son[u], 1);
    add(u, son[u], 1);
    ans[u] = now;
    if (!opt) {
        add(u, 0, 0);
    }
}


int main() {
    ios::sync_with_stdio(false),
            cin.tie(nullptr),
            cout.tie(nullptr);
    int t;
    cin >> t;
    int test = 0;
    while (t--) {
        int n;
        cin >> n;
        for (int i = 1; i < n; i++) {
            int a, b;
            cin >> a >> b;
            g[a].push_back(b);
            g[b].push_back(a);
        }
        cout << "Case #" << ++test << ": ";
        dfs1(1, -1);
        dfs(1, 0);
        for (int i = 1; i <= n; i++) {
            if (i != 1) cout << ' ';
            cout << ans[i];
        }
        cout << endl;
        for (int i = 1; i <= n; i++) g[i].clear();
    }
}

