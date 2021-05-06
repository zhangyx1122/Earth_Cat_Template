ll fac[N];// n!
ll invfac[N]; // n!的inv
ll invn[N]; //n的inv

inline void init() {
    fac[0] = fac[1] = invfac[0] = invfac[1] = invn[0] = invn[1] = 1;
    for (int i = 2; i < N; ++i) {
        fac[i] = fac[i - 1] * i % mod;
        invn[i] = (mod - mod / i) * invn[mod % i] % mod;
        invfac[i] = invfac[i - 1] * invn[i] % mod;
    }
}

ll C(ll up, ll down) {
    if (up > down) return 0;
    if (up < 0 || down < 0) return 0;
    ll res = fac[down];
    res = res * invfac[down - up] % mod;
    res = res * invfac[up] % mod;
    return res;
}

//先init
