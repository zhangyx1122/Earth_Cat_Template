#include <bits/stdc++.h>

using namespace std;
typedef long long ll;
const ll N = 2e5 + 10;

ll age[N];
struct edge {
    ll to, val;
};

struct father {
    ll u, num;
    ll dist;
};

struct son {
    ll age, dist;

    bool operator<(const son &s) const {
        return age < s.age;
    }
};

vector<father> f[N];
vector<vector<son> > s[N];
vector<edge> g[N];
bool st[N];
ll siz[N];

ll getsiz(ll u, ll fa) {
    if (st[u]) return 0;
    siz[u] = 1;
    for (auto x:g[u]) {
        if (x.to == fa) continue;
        if (st[x.to]) continue;
        siz[u] += getsiz(x.to, u);
    }
    return siz[u];
}

void getwc(ll u, ll fa, ll tot, ll &wc) {
    if (st[u]) return;
    ll mmax = 0, sum = 1;
    for (auto x:g[u]) {
        if (x.to == fa) continue;
        if (st[x.to]) continue;
        getwc(x.to, u, tot, wc);
        mmax = max(mmax, siz[x.to]);
        sum += siz[x.to];
    }
    mmax = max(mmax, tot - sum);
    if (2 * mmax <= tot) wc = u;
}

void getdist(ll u, ll fa, ll now, ll rt, ll kth, vector<son> &v) {
    if (st[u]) return;
    f[u].push_back({rt, kth, now});
    v.push_back({age[u], now});
    for (auto x:g[u]) {
        if (x.to == fa || st[x.to]) continue;
        getdist(x.to, u, now + x.val, rt, kth, v);
    }
}

void calc(ll u) {
    if (st[u]) return;
    getwc(u, -1, getsiz(u, -1), u);

    st[u] = 1;

    for (auto x: g[u]) {
        if (st[x.to]) continue;
        s[u].push_back(vector<son>(0));
        auto &v = s[u].back();
        v.push_back({-0x3f3f3f3f, 0});
        v.push_back({0x3f3f3f3f, 0});
        getdist(x.to, u, x.val, u, (ll) s[u].size() - 1, v);
        sort(v.begin(), v.end(), [](son a, son b) { return a.age < b.age; });
        for (ll i = 1; i < v.size(); i++) {
            v[i].dist += v[i - 1].dist;
        }
    }
    for (auto x:g[u]) {
        calc(x.to);
    }
}

ll query(ll u, ll l, ll r) {
    ll ans = 0;
    for (auto x:f[u]) {
        if (l <= age[x.u] && age[x.u] <= r) ans += x.dist;
        for (ll i = 0; i < s[x.u].size(); i++) {
            if (i == x.num) continue;
            auto &v = s[x.u][i];
            ll btn = lower_bound(v.begin(), v.end(), (son) {l, 0}) - v.begin() - 1;
            ll top = upper_bound(v.begin(), v.end(), (son) {r, 0}) - v.begin() - 1;
            ans += v[top].dist - v[btn].dist;
            ans += (top - btn) * x.dist;
        }
    }
    for (auto v:s[u]) {
        ll btn = lower_bound(v.begin(), v.end(), (son) {l, 0}) - v.begin() - 1;
        ll top = upper_bound(v.begin(), v.end(), (son) {r, 0}) - v.begin() - 1;
        ans += v[top].dist - v[btn].dist;

    }
    return ans;
}

signed main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    cout.tie(nullptr);

    ll n, q, a;
    cin >> n >> q >> a;
    for (ll i = 1; i <= n; i++) cin >> age[i];
    for (ll i = 1; i < n; i++) {
        ll x, y, z;
        cin >> x >> y >> z;
        g[x].push_back({y, z});
        g[y].push_back({x, z});
    }

    calc(1);

    ll ans = 0;
    while (q--) {
        ll u, l, r;
        cin >> u >> l >> r;
        l = (l + ans) % a;
        r = (r + ans) % a;
        if (l > r) swap(l, r);
        ans = query(u, l, r);
        cout << ans << endl;
    }
}
