#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
const ll mod = 1e9 + 7;
struct Matrix {
    static const ll MAXN = 300;
    ll a[MAXN][MAXN];

    void init() { memset(a, 0, sizeof(a)); }

    ll det(ll n) {
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) a[i][j] = (a[i][j] + mod) % mod;
        ll res = 1;
        for (int i = 0; i < n; i++) {
            if (!a[i][i]) {
                bool flag = false;
                for (int j = i + 1; j < n; j++) {
                    if (a[j][i]) {
                        flag = true;
                        for (int k = i; k < n; k++) {
                            swap(a[i][k], a[j][k]);
                        }
                        res = -res;
                        break;
                    }
                }
                if (!flag) return 0;
            }

            for (int j = i + 1; j < n; j++) {
                while (a[j][i]) {
                    ll t = a[i][i] / a[j][i];
                    for (int k = i; k < n; k++) {
                        a[i][k] = (a[i][k] - t * a[j][k]) % mod;
                        swap(a[i][k], a[j][k]);
                    }
                    res = -res;
                }
            }
            res *= a[i][i];
            res %= mod;
        }
        return (res + mod) % mod;
    }
} mat;
