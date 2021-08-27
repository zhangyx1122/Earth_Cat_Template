#include <bits/stdc++.h>
using namespace std;

typedef long long ll;

const int maxn = 20;

ll A[maxn], B[maxn];

ll exgcd(ll a, ll b, ll & x, ll & y) {
	if(b == 0) {
		x = 1, y = 0;
		return a;
	}
	
	ll d = exgcd(b, a % b, y, x);
	
	y -= (a / b) * x;
	
	return d; 
}

int main() {
	int n;
	cin >> n;
	ll M = 1ll;
	for(int i = 0; i < n; ++ i) {
		cin >> A[i] >> B[i];
		M = M * A[i];
	}
	
	ll ans = 0;
	
	ll x, y;
	
	for(int i = 0; i < n; ++ i) {
		ll Mi = M / A[i];
		exgcd(Mi, A[i], x, y);
		ans += B[i] * Mi * x;
	} 
	
	cout << (ans % M + M) % M;
	
} 
