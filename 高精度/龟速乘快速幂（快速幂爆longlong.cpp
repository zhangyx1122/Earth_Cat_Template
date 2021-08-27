#include <bits/stdc++.h>
using namespace std;

typedef long long ll;

ll qmul(ll a, ll b, ll p) {
	ll res = 0;
	while(b) {
		if(b & 1) res = (res + a) % p;
		a = (a + a) % p;
		b >>= 1; 
	}
	return res;
}

ll qpow(ll x, ll n, ll p) {
	ll res = 1;
	while(n) {
		if(n & 1) res = qmul(res, x, p);
		x = qmul(x, x, p);
		n >>= 1;
	}
	return res % p; // 1 0 1
}

int main() {
	ll b, p, k;
	cin >> b >> p >> k;
	ll ans = qpow(b, p, k);
	printf("%lld^%lld mod %lld=%lld", b, p, k, ans);
	
	return 0;
} 
