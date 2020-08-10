#include <bits/stdc++.h>
using namespace std;
struct SCC {
    static const int MAXN = 5000;
    static const int MAXM = 2000000;
    int dfs_clock, edge_cnt = 1, scc_cnt;
    int head[MAXN];
    int dfn[MAXN], lowlink[MAXN];
    int sccno[MAXN];
    stack<int> s;

    struct edge {
        int v, next;
    } e[MAXM];

    void add_edge(int u, int v) {
        e[edge_cnt].v = v;
        e[edge_cnt].next = head[u];
        head[u] = edge_cnt++;
    }

    void tarjan(int u) {
        int v;
        dfn[u] = lowlink[u] = ++dfs_clock;  //每次dfs，u的次序号增加1
        s.push(u);                          //将u入栈
        for (int i = head[u]; i != -1; i = e[i].next)  //访问从u出发的边
        {
            v = e[i].v;
            if (!dfn[v])  //如果v没被处理过
            {
                tarjan(v);  // dfs(v)
                lowlink[u] = min(lowlink[u], lowlink[v]);
            } else if (!sccno[v])
                lowlink[u] = min(lowlink[u], dfn[v]);
        }
        if (dfn[u] == lowlink[u]) {
            scc_cnt++;
            do {
                v = s.top();
                s.pop();
                sccno[v] = scc_cnt;
            } while (u != v);
        }
    }

    int find_scc(int n) {
        for (int i = 1; i <= n; i++)
            if (!dfn[i]) tarjan(i);
        return scc_cnt;
    }

    void init() {
        scc_cnt = dfs_clock = 0;
        edge_cnt = 1;  //不用初始化e数组，省时间
        while (!s.empty()) s.pop();
        memset(head, -1, sizeof(head));
        memset(sccno, 0, sizeof(sccno));
        memset(dfn, 0, sizeof(dfn));
        memset(lowlink, 0, sizeof(lowlink));
    }
} scc;
