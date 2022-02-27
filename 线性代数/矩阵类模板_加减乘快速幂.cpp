#include <bits/stdc++.h>

using namespace std;
typedef long long ll;
const ll N = 305;
const ll mod = 998244353;

//矩阵类模板
struct Matrix {
    ll n, m;
    ll a[N][N];

    void set(ll _a, ll _b) {
        n = _a, m = _b;
    }

    Matrix() {
        clear();
    }

    void clear() {
        n = m = 0;
        memset(a, 0, sizeof(a));
    }

    Matrix operator+(const Matrix &b) const {
        Matrix tmp;
        tmp.n = n;
        tmp.m = m;
        for (ll i = 0; i < n; ++i)
            for (ll j = 0; j < m; ++j)
                tmp.a[i][j] = (a[i][j] + b.a[i][j]) % mod;
        return tmp;
    }

    Matrix operator-(const Matrix &b) const {
        Matrix tmp;
        tmp.n = n;
        tmp.m = m;
        for (ll i = 0; i < n; ++i) {
            for (ll j = 0; j < m; ++j)
                tmp.a[i][j] = (a[i][j] - b.a[i][j] + mod) % mod;
        }

        return tmp;
    }

    Matrix operator*(const Matrix &b) const {
        Matrix tmp;
        tmp.clear();
        tmp.n = n;
        tmp.m = b.m;
        for (ll i = 0; i < n; ++i)
            for (ll j = 0; j < b.m; ++j)
                for (ll k = 0; k < m; ++k) {
                    tmp.a[i][j] += a[i][k] * b.a[k][j];
                    tmp.a[i][j] %= mod;
                }
        return tmp;
    }

    Matrix get(ll x) {//幂运算
        Matrix E;
        E.clear();
        E.set(n, m);
        for (ll i = 0; i < n; ++i)
            E.a[i][i] = 1;
        if (x == 0) return E;
        else if (x == 1) return *this;
        Matrix tmp = get(x / 2);
        tmp = tmp * tmp;
        if (x % 2) tmp = tmp * (*this);
        return tmp;
    }

    void exgcd(ll _a, ll _b, ll &x, ll &y) {
        if (!_b)return x = 1, y = 0, void();
        exgcd(_b, _a % _b, y, x);
        y -= x * (_a / _b);
    }

    ll inv(ll p) {
        ll x, y;
        exgcd(p, mod, x, y);
        return (x + mod) % mod;
    }

    Matrix inv() {
        Matrix E = *this;
        ll is[N], js[N];
        for (ll k = 0; k < E.n; k++) {
            is[k] = js[k] = -1;
            for (ll i = k; i < E.n; i++) // 1
                for (ll j = k; j < E.n; j++)
                    if (E.a[i][j]) {
                        is[k] = i, js[k] = j;
                        break;
                    }
            if (is[k] == -1) {
                E.clear();
                return E;
            }
            for (ll i = 0; i < E.n; i++) // 2
                swap(E.a[k][i], E.a[is[k]][i]);
            for (ll i = 0; i < E.n; i++)
                swap(E.a[i][k], E.a[i][js[k]]);
            if (!E.a[k][k]) {
                E.clear();
                return E;
            }
            E.a[k][k] = inv(E.a[k][k]); // 3
            for (ll j = 0; j < E.n; j++)
                if (j != k) // 4
                    (E.a[k][j] *= E.a[k][k]) %= mod;
            for (ll i = 0; i < E.n; i++)
                if (i != k) // 5
                    for (ll j = 0; j < E.n; j++)
                        if (j != k)
                            (E.a[i][j] += mod - E.a[i][k] * E.a[k][j] % mod) %= mod;
            for (ll i = 0; i < E.n; i++)
                if (i != k) // 就是这里不同
                    E.a[i][k] = (mod - E.a[i][k] * E.a[k][k] % mod) % mod;
        }
        for (ll k = E.n - 1; k >= 0; k--) { // 6
            for (ll i = 0; i < E.n; i++)
                swap(E.a[js[k]][i], E.a[k][i]);
            for (ll i = 0; i < E.n; i++)
                swap(E.a[i][is[k]], E.a[i][k]);
        }
        return E;
    }
};
//矩阵模板结束
