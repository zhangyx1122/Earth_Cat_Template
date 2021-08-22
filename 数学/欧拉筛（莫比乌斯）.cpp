#include <bits/stdc++.h>

using namespace std;
typedef long long ll;
const int N = 1e5 + 10;

bool vis[N];
ll prime[N], mu[N];

void init_mu() {
    ll cnt = 0;
    mu[1] = 1;
    for (ll i = 2; i < N; i++) {
        if (!vis[i]) {
            prime[cnt++] = i;
            mu[i] = -1;
        }
        for (ll j = 0; j < cnt && i * prime[j] < N; j++) {
            vis[i * prime[j]] = 1;
            if (i % prime[j] == 0) {
                mu[i * prime[j]] = 0;
                break;
            } else { mu[i * prime[j]] = -mu[i]; }
        }
    }
}

int main() {
    init_mu();
}
