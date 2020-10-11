### fread快读

```c++
#include <bits/stdc++.h>
using namespace std;

char next_char() {
	static char buf[1 << 20], *first, *last;
	if(first == last) {
		last = buf + fread(buf, 1, 1 << 20, stdin);
		first = buf;
	}
	return first == last ? EOF : *first ++;
}

inline int read(){
	int x = 0, w = 0; char ch = 0;
	while(!isdigit(ch)) {w |= ch == '-'; ch = next_char(); }
	while(isdigit(ch)) {x = (x << 3) + (x << 1) + (ch ^ 48), ch = next_char(); }
	return w ? -x : x;
}

int main(){
	freopen("1.txt", "r", stdin); // 交代码的时候一定要去掉aaa 
	int T;
	cin >> T;
	while(T --){
		int x = read();
		cout << x << endl;
	}
} 
```

---

### lucas求组合数

```c++
#include <bits/stdc++.h>
using namespace std;

typedef long long ll;

ll p;

const int maxn = 1e5 + 10;

ll qpow(ll x, ll n){
	ll res = 1;
	while(n){
		if(n & 1) res = (res * x) % p;
		x = (x * x) % p;
		n >>= 1;
	}
	
	return res;
}

ll C(ll up, ll down){
	if(up > down) return 0;
	ll res = 1;

//	for(int i = up + 1; i <= down; ++ i){
//		res = (res * i) % p;
//	}
//	for(int i = 1; i <= down - up; ++ i){
//		res = (res * qpow(i, p - 2)) % p; 
//	}

	for(int i = 1, j = down; i <= up; ++ i, -- j){
		res = (res * j) % p;
		res = (res * qpow(i, p - 2)) % p;
	}
	
	return res;
}


ll lucas(ll up, ll down){
	if(up < p && down < p) return C(up, down);
	return C(up % p, down % p) * lucas(up / p, down / p) % p; 
}

int main(){
	ios::sync_with_stdio(0); cin.tie(0); cout.tie(0);
	
	int T;
	cin >> T;
	while (T --){
		ll down, up;
		cin >> down >> up >> p;
		
		cout << lucas(up, down) % p << endl;
	}
	
	return 0;
} 
```

