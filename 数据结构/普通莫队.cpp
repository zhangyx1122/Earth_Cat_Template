#include <bits/stdc++.h>
using namespace std;

const int N = 1e6 + 10, M = 1e6 + 10;
int a[N];

struct node { 
	int id, l, r;
} mp[M];

int len;
int ans[M], cnt[1000010];

int getNum(int l) {
	return l / len;
}

//左指针的分块，右指针的大小
//有个奇数分块右指针从大到小，偶数分块右指针从小到大的优化 
bool cmp (const node &a, const node & b) {
	if(getNum(a.l) == getNum(b.l)) return a.r < b.r;
	return a.l < b.l;
}

void add(int x, int& res) {
	if(cnt[x] == 0) res++;
	cnt[x] ++;
}

void del(int x, int& res) {
	cnt[x] --;
	if(cnt[x] == 0) res --;
}

int main() {
	ios::sync_with_stdio(0); cin.tie(0); cout.tie(0);
	
	int n;
	cin >> n;
	for(int i = 1; i <= n; ++ i) {
		cin >> a[i];
	}
	int m;
	cin >> m;
	len = sqrt((double)n * n / m);
	for(int i = 1; i <= m; ++ i) {
		mp[i].id = i;
		cin >> mp[i].l >> mp[i].r;
	}
	sort(mp + 1, mp + m + 1, cmp);
	
	//离线处理询问 
	int res = 0, i = 0, j = 0;
	for(int k = 1; k <= m; ++ k) {
		int id = mp[k].id, l = mp[k].l, r = mp[k].r;
		while(j < r) add(a[++j], res);
		while(j > r) del(a[j--], res);
		while(i < l) del(a[i++], res);
		while(i > l) add(a[--i], res);
		ans[id] = res;
	}
	
	for(int i = 1; i <= m; ++ i) {
		cout << ans[i] << endl;
	}
	return 0;
} 
