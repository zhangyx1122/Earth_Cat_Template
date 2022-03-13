#include<bits/stdc++.h>

using namespace std;

ll qpow(ll n, ll m, ll p) {
    ll ans = 1;
    n %= p;
    for (; m; m >>= 1, n = n * n % p)
        if (m & 1)ans = ans * n % p;
    return ans;
}

ll pri[N], tot;

ll getRoot(ll p) //求质数p的最小原根
{
    tot = 0;
    ll n = p - 1, sq = sqrt(p + 0.5);
    for (ll i = 2; i <= sq; i++)
        if (n % i == 0) {
            pri[tot++] = i;
            while (n % i == 0)n /= i;
        }
    if (n > 1)pri[tot++] = n;
    for (ll g = 2; g <= p - 1; g++) //试探每一个g是否原根
    {
        ll flag = 1;
        for (ll i = 0; i < tot; i++)
            if (qpow(g, (p - 1) / pri[i], p) == 1) {
                flag = 0;
                break;
            }
        if (flag)return g;
    }
    return -1; //没有原根
}
