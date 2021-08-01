#include<bits/stdc++.h>

using namespace std;
typedef long long ll;
const ll maxN = 310;
const ll INF = 1e16;

struct KM {
    ll mp[maxN][maxN], link_x[maxN], link_y[maxN], N;
    bool visx[maxN], visy[maxN];
    ll que[maxN << 1], top, fail, pre[maxN];
    ll hx[maxN], hy[maxN], slk[maxN];

    inline ll check(ll i) {
        visx[i] = true;
        if (link_x[i]) {
            que[fail++] = link_x[i];
            return visy[link_x[i]] = true;
        }
        while (i) {
            link_x[i] = pre[i];
            swap(i, link_y[pre[i]]);
        }
        return 0;
    }

    void bfs(ll S) {
        for (ll i = 1; i <= N; i++) {
            slk[i] = INF;
            visx[i] = visy[i] = false;
        }
        top = 0;
        fail = 1;
        que[0] = S;
        visy[S] = true;
        while (true) {
            ll d;
            while (top < fail) {
                for (ll i = 1, j = que[top++]; i <= N; i++) {
                    if (!visx[i] && slk[i] >= (d = hx[i] + hy[j] - mp[i][j])) {
                        pre[i] = j;
                        if (d) slk[i] = d;
                        else if (!check(i)) return;
                    }
                }
            }
            d = INF;
            for (ll i = 1; i <= N; i++) {
                if (!visx[i] && d > slk[i]) d = slk[i];
            }
            for (ll i = 1; i <= N; i++) {
                if (visx[i]) hx[i] += d;
                else slk[i] -= d;
                if (visy[i]) hy[i] -= d;
            }
            for (ll i = 1; i <= N; i++) {
                if (!visx[i] && !slk[i] && !check(i)) return;
            }
        }
    }

    void init() {
        for (ll i = 1; i <= N; i++) {
            link_x[i] = link_y[i] = 0;
            visy[i] = false;
        }
        for (ll i = 1; i <= N; i++) {
            hx[i] = 0;
            for (ll j = 1; j <= N; j++) {
                if (hx[i] < mp[i][j]) hx[i] = mp[i][j];
            }
        }
    }
} km;

int main() {
    ios::sync_with_stdio(0);

    ll n;
    cin >> n;
    ll ans = 0;
    for (int i = 1; i <= n; i++) {
        ll a, b, c, d;
        cin >> a >> b >> c >> d;
        ans += a * a + b * b;
        for (int j = 1; j <= n; j++) {
            km.mp[i][j] = -(c + d * (j - 1)) * (c + d * (j - 1));
//            cout << -km.mp[i][j] << ' ';
//            cin >> km.mp[i][j];
//            km.mp[i][j] = -km.mp[i][j];
        }
//        cout << endl;
    }
    km.N = n;
    km.init();
    for (int i = 1; i <= km.N; i++) km.bfs(i);
    for (int i = 1; i <= n; i++) ans -= km.mp[i][km.link_x[i]];
    cout << ans << endl;
}
