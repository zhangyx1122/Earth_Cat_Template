/*
 仙人掌:任意一条边至多只出现在一条简单回路的无向连通图称为仙人掌。
 转化为圆方树，然后根据树的算法来做一些问题，注意区分圆点和方点
 这题:求带环（环和环之间无公共边）无向图两点间的最短路径
 */

#include <iostream>
#include <cstring>
#include <algorithm>

using namespace std;

const int N = 12010, M = N * 3;

int n, m, Q, new_n;
int h1[N], h2[N], e[M], w[M], ne[M], idx;
int dfn[N], low[N], cnt;
int s[N], stot[N], fu[N], fw[N];
int fa[N][14], depth[N], d[N];
int A, B;

void add(int h[], int a, int b, int c)
{
    e[idx] = b, w[idx] = c, ne[idx] = h[a], h[a] = idx ++ ;
}

void build_circle(int x, int y, int z)
{
    int sum = z;
    for (int k = y; k != x; k = fu[k])
    {
        s[k] = sum;
        sum += fw[k];
    }
    s[x] = stot[x] = sum;
    add(h2, x, ++ new_n, 0);
    for (int k = y; k != x; k = fu[k])
    {
        stot[k] = sum;
        add(h2, new_n, k, min(s[k], sum - s[k]));
    }
}

void tarjan(int u, int from)
{
    dfn[u] = low[u] = ++ cnt;
    for (int i = h1[u]; ~i; i = ne[i])
    {
        int j = e[i];
        if (!dfn[j])
        {
            fu[j] = u, fw[j] = w[i];
            tarjan(j, i);
            low[u] = min(low[u], low[j]);
            if (dfn[u] < low[j]) add(h2, u, j, w[i]);
        }
        else if (i != (from ^ 1)) low[u] = min(low[u], dfn[j]);
    }
    for (int i = h1[u]; ~i; i = ne[i])
    {
        int j = e[i];
        if (dfn[u] < dfn[j] && fu[j] != u)
            build_circle(u, j, w[i]);
    }
}

void dfs_lca(int u, int father)
{
    depth[u] = depth[father] + 1;
    fa[u][0] = father;
    for (int k = 1; k <= 13; k ++ )
        fa[u][k] = fa[fa[u][k - 1]][k - 1];
    for (int i = h2[u]; ~i; i = ne[i])
    {
        int j = e[i];
        d[j] = d[u] + w[i];
        dfs_lca(j, u);
    }
}

int lca(int a, int b)
{
    if (depth[a] < depth[b]) swap(a, b);
    for (int k = 13; k >= 0; k -- )
        if (depth[fa[a][k]] >= depth[b])
            a = fa[a][k];
    if (a == b) return a;
    for (int k = 13; k >= 0; k -- )
        if (fa[a][k] != fa[b][k])
        {
            a = fa[a][k];
            b = fa[b][k];
        }
    A = a, B = b;
    return fa[a][0];
}

int main()
{
    scanf("%d%d%d", &n, &m, &Q);
    new_n = n;
    memset(h1, -1, sizeof h1);
    memset(h2, -1, sizeof h2);
    while (m -- )
    {
        int a, b, c;
        scanf("%d%d%d", &a, &b, &c);
        add(h1, a, b, c), add(h1, b, a, c);
    }
    tarjan(1, -1);
    dfs_lca(1, 0);

    while (Q -- )
    {
        int a, b;
        scanf("%d%d", &a, &b);
        int p = lca(a, b);
        if (p <= n) printf("%d\n", d[a] + d[b] - d[p] * 2);
        else
        {
            int da = d[a] - d[A], db = d[b] - d[B];
            int l = abs(s[A] - s[B]);
            int dm = min(l, stot[A] - l);
            printf("%d\n", da + dm + db);
        }
    }

    return 0;
}
