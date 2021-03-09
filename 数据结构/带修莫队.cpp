#include <bits/stdc++.h>
using namespace std;

const int N = 10010;

int a[N], cnt[1000010], ans[N];

int len, mq, mc;

struct Query {
	int id, l, r, t;
} q[N];

struct Modify {
	int p, c;
} c[N];

int getNum(int x) {
	return x / len;
}

// l所在块的编号，r所在块的编号，t升序

bool cmp(const Query& a, const Query& b) {
	if(getNum(a.l) == getNum(b.l) && getNum(a.r) == getNum(b.r)) {
		return a.t < b.t;
	} 
	if(getNum(a.l) == getNum(b.l)) return a.r < b.r;
	return a.l < b.l; 
}

void add(int x, int& res) {
    if (!cnt[x]) res ++ ;
    cnt[x] ++ ;
}

void del(int x, int& res) {
    cnt[x] -- ;
    if (!cnt[x]) res -- ;
}


int main() {
	ios::sync_with_stdio(0); cin.tie(0); cout.tie(0);
	
	int n, m;
	cin >> n >> m;
	char op;
	int x, y;
	for(int i = 1; i <= n; ++ i) {
		cin >> a[i];
	}
	for(int i = 1; i <= m; ++ i) {
		cin >> op >> x >> y;
        if (op == 'Q') q[++ mq] = {mq, x, y, mc};
        else c[ ++ mc] = {x, y};
	}
	
  ///
	len = cbrt((double)n * mc) + 1;
  sort(q + 1, q + mq + 1, cmp);
	
	int i = 1, j = 0, t = 0, res = 0;
	for(int k = 1; k <= mq; ++ k) {
		int id = q[k].id, l = q[k].l, r = q[k].r, tm = q[k].t;
		while(j < r) add(a[++ j], res);
		while(j > r) del(a[j --], res);
		while(i < l) del(a[i ++], res);
		while(i > l) add(a[-- i], res);
		while(t < tm) {
			++ t;
			if(c[t].p >= i && c[t].p <= j) {
				del(a[c[t].p], res);
				add(c[t].c, res);
			}
			swap(a[c[t].p], c[t].c);
		}
		while(t > tm) {
			if(c[t].p >= i && c[t].p <= j) {
				del(a[c[t].p], res);
				add(c[t].c, res);
			}
			swap(a[c[t].p], c[t].c);
			-- t;
		}
		ans[id] = res;
	}
	
	for(int i = 1; i <= mq; ++ i) {
		cout << ans[i] << endl;
	}
}
