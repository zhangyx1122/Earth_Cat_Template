#include <bits/stdc++.h>

using namespace std;
#define de(x) cout << #x << " = " << x << endl
#define dd(x) cout << #x << " = " << x << " "
typedef long long ll;
const __int128 N = 1e6 + 10;
const __int128 M = 41;
const __int128 mod = 4179340454199820289ll;
const __int128 inv2 = 2089670227099910145ll;

__int128 isp[N], pri[N], pcnt;
vector<__int128> h[N];

void getPrime() {
    fill(isp + 2, isp + N, 1);
    for (__int128 i = 2; i < N; i++) {
        if (isp[i]) {
            pri[pcnt++] = i;
        }
        for (__int128 j = 0; j < pcnt && i * pri[j] < N; j++) {
            isp[i * pri[j]] = 0;
            if (i % pri[j] == 0) break;
        }
    }
}

__int128 qpow(__int128 x, __int128 y) {
    __int128 ans = 1;
    while (y) {
        if (y & 1) ans = ans * x % mod;
        x = x * x % mod;
        y >>= 1;
    }
    return ans;
}

__int128 g(__int128 x) {
    return x;
}

__int128 G(__int128 x) {
    return x * (x + 1) % mod * inv2 % mod;
}

__int128 f(__int128 x, __int128 c) {
    return x * qpow(c, mod - 2) % mod;
}

__int128 ans;
ll n;

void dfs(__int128 deep, __int128 hpn, __int128 pn, bool flag) {
    if (flag) {
        ans = (ans + hpn * G(n / pn)) % mod;
    }
    if (deep >= pcnt) return;
    if (pri[deep] * pri[deep] * pn > n) return;
    dfs(deep + 1, hpn, pn, false);
    for (__int128 i = 2, pi = pri[deep] * pri[deep] % mod; pn * pi <= n; i++, pi = pi * pri[deep] % mod) {
        dfs(deep + 1, hpn * h[deep][i] % mod, pn * pi, true);
    }
}

signed main() {
    ios::sync_with_stdio(0);
    getPrime();
    for (__int128 pid = 0; pid < pcnt; pid++) {
        h[pid].push_back(1);
        h[pid].push_back(0);
        __int128 invp = qpow(pri[pid], mod - 2);
        for (__int128 c = 2, pc = pri[pid] * pri[pid]; c < M && pc <= 1e12; c++, pc = pc * pri[pid]) {
            h[pid].push_back(f(pri[pid], c));//唯一f使用，传入参数类型自定义
            __int128 pci = qpow(pri[pid], c);
            for (__int128 i = 1, pi = pri[pid]; i <= c; i++, pi = pi * pri[pid] % mod) {
                pci = pci * invp % mod;
                h[pid][c] = (h[pid][c] - g(pi) * h[pid][c - i] % mod + mod) % mod;
            }
        }
    }

    ll t;
    cin >> t;
    while (t--) {
        cin >> n;
        ans = G(n);
        dfs(0, 1, 1, false);
        cout << (ll) ans << endl;
    }
}
