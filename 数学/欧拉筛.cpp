#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
const int N = 1000005;
int phi[N], prime[N], cnt;
bool st[N];

void get_eulers() {
    phi[1] = 1;
    for (int i = 2; i < N; i++) {
        if (!st[i]) {
            prime[cnt++] = i;
            phi[i] = i - 1;
        }
        for (int j = 0; prime[j] * i < N; j++) {
            st[prime[j] * i] = 1;
            if (i % prime[j] == 0) {
                phi[prime[j] * i] = phi[i] * prime[j];
                break;
            }
            phi[prime[j] * i] = phi[i] * (prime[j] - 1);
        }
    }
}

int main() {
    get_eulers();
    ll n;
    cin >> n;
    ll ans = 0;
    for (int i = 1; i <= n; i++) ans += phi[i];
    cout << ans;
}
