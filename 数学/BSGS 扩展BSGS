# BSGS

求$a^t \equiv b ( \mod p)$ (a,p) = 1的最小的t

$t = x \times k - y, x \in[1, k], y \in[0, k -1 ]$
$t \in [1,k^2]$

$a^kx \equiv b \times a^y (\mod p)$

对 $b \times a^y$ 建立hash表，枚举x看是否有解

```cpp
#include <bits/stdc++.h>
using namespace std;

typedef long long ll;

unordered_map<int , int> mp;

int bsgs(int a, int p, int b) {
	
	if (1 % p == b % p) return 0; // 特判0是不是解
	mp.clear();
	
	int k = sqrt(p) + 1;
	
	for(int i = 0, j = b % p; i < k; ++ i, j = (ll)j * a % p) {
		mp[j] = i;
	}
	
	int ak = 1;
	for(int i = 0; i < k; ++i) {
		ak = (ll)ak * a % p;
	}
	
	for(int i = 1, j = ak % p; i <= k; ++ i, j = (ll)j * ak % p) {
		if(mp.count(j)) return (ll)i * k - mp[j];
	}
	
	return -1;
}

int main() {
	ios::sync_with_stdio(0);
	cin.tie(0); cout.tie(0);
	
	int a, p, b;
	while(cin >> a >> p >> b, a | p | b) {
		int res;
		res = bsgs(a, p, b);
		if(res == -1) {
			cout << "No Solution\n"; 
		}
		else {
			cout << res << endl;
		}
	}
	
	return 0;
}
```

# 扩展BSGS

求$a^t \equiv b(\mod p)$ 的最小的t

当$(a, p) \,!= 1$

$(a, p) = d$ $d \nmid b$ 无解

$a^t \equiv b(\mod p)$ ，$a^{t} + kp = b$  两边同时除以d， $\frac{a}{d}a^{t - 1} + k \frac{p}{d} = \frac{b}{d}$

$a^{t - 1} \equiv \frac{b}{d}(\frac{a}{d})^{-1}$

$t' = t - 1, p' = \frac{p}{d}, b' = \frac{b}{a}(\frac{a}{d})^{-1}$

```cpp
#include <bits/stdc++.h>
using namespace std;

typedef long long ll;

unordered_map<ll, ll> mp;

ll bsgs(ll a, ll p, ll b) {
	
	if(1 % p == b % p) return 0; // 特判0是不是解
	mp.clear(); 
	
	ll k = sqrt(p) + 1;
	
	for(ll i = 0, j = b % p; i < k; ++i, j = (ll)j * a % p) {
		mp[j] = i;
	}
	
	ll ak = 1;
	for(ll i = 0; i < k; ++i) {
		ak = (ll) ak * a % p;
	}
	
	for(ll i = 1, j = ak % p;i <= k; ++i, j = (ll)j * ak % p) {
		if(mp.count(j)) return (ll) i * k - mp[j];
	}
	
	return -1;
}

ll gcd(ll x, ll y) {
	return x % y == 0 ? y : gcd(y, x % y); 
}

void extgcd(ll a,ll b,ll& d,ll& x,ll& y){
    if(!b){
        d = a; x = 1; y = 0;
    }
    else{ 
        extgcd(b, a%b, d, y, x); 
        y -= x * (a / b); 
    }
}

ll inverse(ll a,ll n){
    ll d,x,y;
    extgcd(a,n,d,x,y);
    return d == 1 ? (x + n) % n : -1;

}

int main() {
	ll a, p, b;
	
	while(cin >> a >> p >> b, a | p | b) {
		ll d = gcd(a, p);
		if(d == 1) {
			ll res = bsgs(a, p, b);
			if(res == -1) {
				cout << "No Solution\n";
			}
			else {
				cout << res << endl;
			}
		}
		else {
			if(b % d != 0) {
				cout << "No Solution\n";
				continue;
			}
			else {
				p = p / d;
				b = (b / d) * inverse(a / d, p);
				ll res = bsgs(a, p, b);
				if(res == -1) {
					cout << "No Solution\n";
				}
				else {
					cout << res + 1 << endl;
				}
			} 
		}
	} 
	
	return 0;
}
```

