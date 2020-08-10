#include <bits/stdc++.h>
using namespace std;
typedef long long ll;

struct Edge {
    ll from, to, cap, flow, cost;
    Edge(ll u, ll v, ll c, ll f, ll w)
        : from(u), to(v), cap(c), flow(f), cost(w) {}
};

struct MCMF {
    static const ll maxn = 6000;
    static const ll INF = 0x3f3f3f3f3f3f3f;
    ll n, m;
    vector<Edge> edges;
    vector<ll> G[maxn];
    ll inq[maxn];
    ll d[maxn];
    ll p[maxn];
    ll a[maxn];

    void init(ll n) {
        this->n = n;
        for (ll i = 1; i <= n; i++) G[i].clear();
        edges.clear();
    }

    void add_edge(ll from, ll to, ll cap, ll cost) {
        edges.push_back(Edge(from, to, cap, 0, cost));
        edges.push_back(Edge(to, from, 0, 0, -cost));
        m = edges.size();
        G[from].push_back(m - 2);
        G[to].push_back(m - 1);
    }

    bool BellmanFord(ll s, ll t, ll& flow, ll& cost) {
        for (ll i = 1; i <= n; ++i) d[i] = INF;
        memset(inq, 0, sizeof(inq));
        d[s] = 0, inq[s] = 1, p[s] = 0, a[s] = INF;
        queue<ll> Q;
        Q.push(s);
        while (!Q.empty()) {
            ll u = Q.front();
            Q.pop();
            inq[u] = 0;
            for (ll i = 0; i < G[u].size(); ++i) {
                Edge& e = edges[G[u][i]];
                if (e.cap > e.flow && d[e.to] > d[u] + e.cost) {
                    d[e.to] = d[u] + e.cost;
                    p[e.to] = G[u][i];
                    a[e.to] = min(a[u], e.cap - e.flow);
                    if (!inq[e.to]) {
                        Q.push(e.to);
                        inq[e.to] = 1;
                    }
                }
            }
        }
        if (d[t] == INF) return false;
        flow += a[t];
        cost += (ll)d[t] * (ll)a[t];
        for (ll u = t; u != s; u = edges[p[u]].from) {
            edges[p[u]].flow += a[t];
            edges[p[u] ^ 1].flow -= a[t];
        }
        return true;
    }

    //需要保证初始网络中没有负权圈
    ll MincostMaxflow(ll s, ll t, ll& cost) {
        ll flow = 0;
        cost = 0;
        while (BellmanFord(s, t, flow, cost))
            ;
        return flow;
    }
} mcmf;  //  若固定流量k，增广时在flow+a>=k的时候只增广k-flow单位的流量，然后终止程序
