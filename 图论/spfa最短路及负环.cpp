#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
const int N = 1 << 20;
struct edge {
    ll to, len;
};

vector<edge> g[N];
ll d[N], cnt[N], vis[N];

bool spfa(ll s, ll n) {
    queue<int> que;
    for (int i = 1; i <= n; i++) {  //防止不连通，全加进去
        que.push(i);
        vis[i] = 1;
    }
    while (!que.empty()) {
        ll p = que.front();
        que.pop();
        vis[p] = 0;
        for (auto x:g[p]) {
            if (d[x.to] > d[p] + x.len) {
                d[x.to] = d[p] + x.len;
                cnt[x.to] = cnt[p] + 1;
                if (!vis[x.to]) {
                    if (cnt[x.to] > n) return 0;
                    vis[x.to] = 1;
                    que.push(x.to);
                }
            }
        }
    }
    return 1;
}