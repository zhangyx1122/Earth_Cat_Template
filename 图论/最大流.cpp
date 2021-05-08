#include <bits/stdc++.h>
using namespace std;
typedef long long ll;

struct Edge {
    ll from, to, cap, flow;
    Edge(ll a, ll b, ll c, ll d) : from(a), to(b), cap(c), flow(d) {}
};

struct Dinic {
    static const ll maxn = 10000;
    static const ll inf = 0x3f3f3f3f3f3f3f3f;
    ll N, M, S, T;
    vector<Edge> edges;
    vector<ll> G[maxn];
    bool vis[maxn];
    ll d[maxn];
    ll cur[maxn];

    void AddEdge(ll from, ll to, ll cap) {
        edges.push_back(Edge(from, to, cap, 0));
        edges.push_back(Edge(to, from, 0, 0));
        M = edges.size();
        G[from].push_back(M - 2);
        G[to].push_back(M - 1);
    }

    bool BFS() {
        memset(vis, 0, sizeof(vis));
        queue<ll> Q;
        Q.push(S);
        d[S] = 0;
        vis[S] = 1;
        while (!Q.empty()) {
            ll x = Q.front();
            Q.pop();
            for (ll i = 0; i < G[x].size(); i++) {
                Edge& e = edges[G[x][i]];
                if (!vis[e.to] && e.cap > e.flow) {
                    vis[e.to] = 1;
                    d[e.to] = d[x] + 1;
                    Q.push(e.to);
                }
            }
        }
        return vis[T];
    }

    ll DFS(ll x, ll a) {
        if (x == T || a == 0) return a;
        ll flow = 0, f;
        for (ll& i = cur[x]; i < G[x].size(); i++) {
            Edge& e = edges[G[x][i]];
            if (d[x] + 1 == d[e.to] &&
                (f = DFS(e.to, min(a, e.cap - e.flow))) > 0) {
                e.flow += f;
                edges[G[x][i] ^ 1].flow -= f;
                flow += f;
                a -= f;
                if (a == 0) break;
            }
        }
        return flow;
    }

    ll Maxflow(ll S, ll T) {
        this->S = S, this->T = T;
        ll flow = 0;
        while (BFS()) {
            memset(cur, 0, sizeof(cur));
            flow += DFS(S, inf);
        }
        return flow;
    }
} MF;

//有源汇上下界最大流，跑完可行流后，s-t的最大流即为答案

//有源汇上下届最小流，不连无穷边，s-t跑最大流，再加上t-s无穷边，再跑最大流，无穷边流量为答案

//最大权闭合子图
//构造一个新的流网络，建一个源点s和汇点t，从s向原图中所有点权为正数的点建一条容量为点权的边，
//从点权为负数的点向t建一条容量为点权绝对值的边，原图中各点建的边都建成容量为正无穷的边。
//然后求从s到t的最小割，再用所有点权为正的权值之和减去最小割，就是我们要求的最大权值和了。

//最大密度子图
//01分数规划
//addedge(S,V,m),addedge(E,1),addedge(V,T,2*g-deg(v)+m)
//h(g)=n*m-maxflow(S,T)

