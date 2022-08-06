#include<bits/stdc++.h>

using namespace std;

const int inf = 0x3f3f3f3f;
const int maxN = 505;

namespace KM {
    int mp[maxN][maxN], link_x[maxN], link_y[maxN], N;
    bool visx[maxN], visy[maxN];
    int que[maxN << 1], top, fail, pre[maxN];
    int hx[maxN], hy[maxN], slk[maxN];

    inline int check(int i) {
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

    void bfs(int S) {
        for (int i = 1; i <= N; i++) {
            slk[i] = inf;
            visx[i] = visy[i] = false;
        }
        top = 0;
        fail = 1;
        que[0] = S;
        visy[S] = true;
        while (true) {
            int d;
            while (top < fail) {
                for (int i = 1, j = que[top++]; i <= N; i++) {
                    if (!visx[i] && slk[i] >= (d = hx[i] + hy[j] - mp[i][j])) {
                        pre[i] = j;
                        if (d) slk[i] = d;
                        else if (!check(i)) return;
                    }
                }
            }
            d = inf;
            for (int i = 1; i <= N; i++) {
                if (!visx[i] && d > slk[i]) d = slk[i];
            }
            for (int i = 1; i <= N; i++) {
                if (visx[i]) hx[i] += d;
                else slk[i] -= d;
                if (visy[i]) hy[i] -= d;
            }
            for (int i = 1; i <= N; i++) {
                if (!visx[i] && !slk[i] && !check(i)) return;
            }
        }
    }

    void prework() {
        for (int i = 1; i <= N; i++) {
            link_x[i] = link_y[i] = 0;
            visy[i] = false;
        }
        for (int i = 1; i <= N; i++) {
            hx[i] = 0;
            for (int j = 1; j <= N; j++) {
                if (hx[i] < mp[i][j]) hx[i] = mp[i][j];
            }
        }
    }

    void init(int n) {
        N = n;
        top = fail = 0;
        for (int i = 1; i <= N; i++) {
            link_x[i] = link_y[i] = visx[i] = visy[i] = pre[i] = hx[i] = hy[i] = slk[i] = 0;
            for (int j = 1; j <= N; j++) {
                mp[i][j] = 0;
            }
        }
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(0);
    cout.tie(0);

    int n, m;
    cin >> n >> m;
    KM::init(max(n, m));
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            cin >> KM::mp[i][j];
        }
    }
    KM::prework();
    int ans = 0;
    for (int i = 1; i <= KM::N; i++) KM::bfs(i);
    for (int i = 1; i <= KM::N; i++) ans += KM::mp[i][KM::link_x[i]];

}

