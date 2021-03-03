# 2020-2021 Winter Petrozavodsk Camp, Belarusian SU Contest (XXI Open Cup, Grand Prix of Belarus)

## AC 7/14  目标：- （罚时

[传送门](https://codeforces.com/gym/102956)

**赛后补题目标** ：



---

### WA：

J :  读错题了...

D：有个数字在debug忘记改回去了

M：想错了

---

多交流...



---

## 需学习的知识点：

G：矩阵树定理：n个不同的点可以构成$n^{n - 2}$棵不同的树

----

# ~

## I - Binary Supersonic Utahraptors

```cpp
#include <bits/stdc++.h>
 
using namespace std;
typedef long long ll;
const ll N = 1 << 20;
const ll N2 = 5005;
const ll mod = 1e9 + 7;
const double dinf = 1e99;
const int inf = 0x3f3f3f3f;
const ll linf = 0x3f3f3f3f3f3f3f3f;
 
 
int main() {
    ios::sync_with_stdio(false), cin.tie(nullptr), cout.tie(nullptr);
    int n, m, k;
    cin >> n >> m >> k;
    int sco = 0;
    for (int i = 1; i <= n; i++) {
        int a;
        cin >> a;
        if (a == 0) sco++;
    }
    for (int i = 1; i <= m; i++) {
        int a;
        cin >> a;
        if (a == 1) sco--;
    }
    cout << abs(sco) << endl;
}
```



## J - Burnished Security Updates

二分图染色

```cpp
#include <bits/stdc++.h>
using namespace std;

const int N = 3e5 + 10;
vector<int> g[N];
bool vis[N];
int co[N];

int bfs(int s) {
	queue<int> que;
	que.push(s);
	vis[s] = 1;
	vector<int> ans;
	ans.push_back(s);
	co[s] = 0;
	
	while(!que.empty()) {
		int tmp = que.front();
		que.pop();
		int siz = g[tmp].size();
		for(int i = 0; i < siz; ++ i) {
			int u = g[tmp][i];
			if(vis[u] && co[u] == co[tmp]) return -1;
			if(vis[u]) continue;
			co[u] = 1 - co[tmp];
			vis[u] = 1;
			que.push(u);
			ans.push_back(u);
		}
	}
	
	int cnt = 0, siz = ans.size();
	for(int i = 0; i < siz; ++ i) {
		if(co[ans[i]] == 0) cnt++;
	}
	return min(cnt, siz - cnt);
}

int main() {
	ios::sync_with_stdio(0); cin.tie(0); cout.tie(0);
	memset(co, -1, sizeof co); 
	int n, m;
	cin >> n >> m;
	int x, y;
	for(int i = 0; i < m; ++ i){
		cin >> x >> y;
		g[x].push_back(y);
		g[y].push_back(x); 
	}  
	
	bool flag = 0;
	int ans = 0;
	for(int i = 1; i <= n; ++ i) {
		if(!vis[i]) {
			int res = bfs(i);
			if(res == -1) {
				flag = 1;
				break;
			}
			ans += res;
		}
	}
	if(flag) {
		cout << -1 << endl;
	}
	else cout << ans << endl; 
}
```



## M - Brilliant Sequence of Umbrellas

**题意：** 构造数列a[] 每个a[i]都在1-n的范围内并且数列的长度至少为$\lceil \frac{2}{3} \sqrt{n} \rceil$

使得$a_i > a_{i - 1} $ $ \gcd(a_i, a_{i - 1}) > \gcd(a_{i - 1}, a_{i - 2})$

**思路：** 构造公约数g[]，令g[] 与他前面两项互质并尽量小，a[i] = g[i] * g[i + 1]

其实有点猜的性质，刚开始构造了一些1， 2， 3， 5， 7， 8， 9， 11 ， 13， 14， 15 ... 看密度差不多但是不知道具体是多少，代码试了下1e12好像刚好（

```cpp
#include <bits/stdc++.h>
 
using namespace std;
typedef long long ll;
const ll N = 1 << 20;
const ll N2 = 5005;
const ll mod = 1e9 + 7;
const double dinf = 1e99;
const ll inf = 0x3f3f3f3f;
const ll linf = 0x3f3f3f3f3f3f3f3f;
 
ll gcd(ll a, ll b) {
    return b == 0 ? a : gcd(b, a % b);
}
 
vector<ll> zhi;
vector<ll> ans;
 
int main() {
    ios::sync_with_stdio(false), cin.tie(nullptr), cout.tie(nullptr);
    ll n;
    cin >> n;
    zhi.push_back(1);
    zhi.push_back(2);
    for (ll i = 3; i <= N * 2; i++) {
        if (gcd(i, zhi[(ll) zhi.size() - 1]) == 1 && gcd(i, zhi[(ll) zhi.size() - 2]) == 1) {
            zhi.push_back(i);
        }
    }
    ans.push_back(1);
    for (ll i = 1; i < zhi.size(); i++) {
        if (zhi[i] * zhi[i - 1] > n) break;
        ans.push_back(zhi[i] * zhi[i - 1]);
    }
    cout << ans.size() << endl;
    for (auto x:ans) cout << x << ' ';
}
```



## G - Biological Software Utilities

**题意：** n个节点的树有多少不同（每个节点都有编号，形状相同编号不同算不同的树）的plausible tree，这个的定义是能通过删一些边变成若干个节点数为2的树

 **思路：** 先算出有多少不同的若干个节点数为2的树的结果，然后算能组合成多少不同的plausible tree， 按顺序把二元的这个树加到原树内 $2^{\frac{n}{2} - 1} \times (1 \times 3 \times 5 \dots n - 1) $

中间缺少的数字用小数据打表得数列然后找规律得最后解（..

$a(n) = 2*n^{n-3}$ 

正解应该是==矩阵树定理==：n个不同的点可以构成$n^{n - 2}$棵不同的树

```cpp
#include<bits/stdc++.h>
using namespace std;
typedef long long ll;

const ll mod=998244353;
ll p(ll x){
    if(x==0)return 1;
    ll res=1,a=x+1;
    x--;
    while(x){
        if(x&1)res=res*a%mod;
        a=a*a%mod;
        x>>=1;
    }
    return res;
}
int main(){
    #ifdef ONLINE_JUDGE
        //std::ios::sync_with_stdio(false);
    #else
        //freopen("in.txt","r",stdin);
        //freopen("out.txt","w",stdout);
    #endif
    ll n;
    scanf("%lld",&n);
    if(n&1)printf("0");
    else{
        ll res=1;
        for(ll i=1;i<n;i+=2)res=res*i%mod;
        ll x=n-2;
        ll a=2;
        while(x>0){
            if(x&1)res=res*a%mod;
            a=a*a%mod;
            x>>=1;
        }
        res=res*p(n/2-1)%mod;
        printf("%lld",res);
    }
    return 0;
}

```



## N - Best Solution Unknown

**题意：** 有n个人排成一排，每个人有一个武力值，两个相邻的人可以进行对战，武力值高者获胜，武力值相等则都有可能获胜，胜者武力值+1，败者淘汰，问有多少人有可能成为冠军

n 1e6

**思路：** 在一个区间内，武力值最高的人一定可以打败其他所有人成为冠军

取区间$[L, R]$，设区间内武力值最高的人所在的位置为m，将区间进一步分成左区间$[L, m - 1]$ 和右区间$[m + 1, R]$ ，以左区间为例，左区间武力最大值一定可以打败左区间其他人，并获得值为区间大小-1的武力值增益，此时他的武力若大于m这个人的武力，他就可以打败整个区间的人，右区间同理

刚开始区间大小为$[1,n]$，逐步缩小区间范围处理



**2 ：** 把所有人的武力值和本身所在位置的编号封装成struct按武力值从小到大排序，从武力值最小的开始向左右扩展，记录左右边界，因为后续的武力值更大，来到已处理过的区间时可以直接全部覆盖，左右边界处理为线性复杂度，



```cpp
#include<bits/stdc++.h>
using namespace std;
typedef long long ll;

struct Stu{
    int num,sc;
}stu[1000005];
bool cmp(Stu a,Stu b){
    return a.sc<b.sc;
}
int a[1000005];
int l[1000005],r[1000005];

void judge(int x){
    bool f=0;
    int lx=l[x]-1;
    if(a[x]+r[x]-l[x]>=a[lx]){
        f=1;
        l[x]=min(lx,l[lx]);
        r[x]=max(r[x],r[lx]);
    }
    int rx=r[x]+1;
    if(a[x]+r[x]-l[x]>=a[rx]){
        f=1;
        l[x]=min(l[x],l[rx]);
        r[x]=max(rx,r[rx]);
    }
    if(f)judge(x);
}
int main(){
    #ifdef ONLINE_JUDGE
        //std::ios::sync_with_stdio(false);
    #else
        //freopen("in.txt","r",stdin);
        //freopen("out.txt","w",stdout);
    #endif
    int n;
    scanf("%d",&n);
    for(int i=1;i<=n;i++){
        scanf("%lld",&a[i]);
        l[i]=i,r[i]=i;
        stu[i].num=i;
        stu[i].sc=a[i];
    }
    sort(stu+1,stu+n+1,cmp);
    a[0]=a[n+1]=0x3f3f3f3f;
    for(int i=1;i<=n;i++){
        judge(stu[i].num);
    }
    int sum=0;
    for(int i=1;i<=n;i++){
        if(l[i]==1&&r[i]==n)sum++;
    }
    printf("%d\n",sum);
    for(int i=1;i<=n;i++){
        if(l[i]==1&&r[i]==n)printf("%d ",i);
    }
    return 0;
}

```





## D - Bank Security Unification

**题意：**  n 个数，相对顺序不变，选择其中k 个数，这几个数的权值为 $\sum _{i = 1}^{k - 1} f_i \& f_{i + 1}$， 求权值最大值

**思路：** dp, dp[i] : 表示到第i个位置并且选了a[i]  的最大可能结果。枚举前面的dp[j]进行转移 这里是$O(n^2)$

每个数转移为二进制比如是num位，同样是num位的数字，后面的这个一定比前面的优，因为他可以用前面的所有资源并且两个num可以直接加上（..

优化对于二进制的每一位，记录最新的dp的位置，（

```cpp
#include <bits/stdc++.h>

using namespace std;
typedef long long ll;
const ll N = 1 << 20;
const ll N2 = 5005;
const ll mod = 1e9 + 7;
const double dinf = 1e99;
const ll inf = 0x3f3f3f3f;
const ll linf = 0x3f3f3f3f3f3f3f3f;
ll a[N];
ll sum[N];
ll dp[42];
ll pos[42];

int main() {
    ios::sync_with_stdio(false), cin.tie(nullptr), cout.tie(nullptr);
    ll n;
    cin >> n;
    for (ll i = 1; i <= n; i++) {
        cin >> a[i];
    }

    for (ll i = 1; i <= n; i++) {
        ll top = 0;
        ll mmax = 0;
        for (ll j = 0; j < 42; j++) {
            mmax = max(mmax, dp[j] + (a[i] & a[pos[j]]));
        }
        ll k = a[i];
        while (k) {
            top++;
            k >>= 1;
        }
        dp[top] = mmax;
        pos[top] = i;
    }
    ll ans = 0;
    for (ll i = 0; i < 42; i++) ans = max(ans, dp[i]);
    cout << ans << endl;
}

```





## C - Brave Seekers of Unicorns

**题意：** 数字小于等于n且满足没有连续三个数字异或和为0的递增序列的个数

**思路：** dp，

对于dp[i] ：计算$y < y \bigoplus i < i $   的和，

```cpp
#include<bits/stdc++.h>
using namespace std;
typedef long long ll;

const ll mod=998244353;
ll a[1000005],b[1000005];
int main(){
    #ifdef ONLINE_JUDGE
        //std::ios::sync_with_stdio(false);
    #else
        //freopen("in.txt","r",stdin);
        //freopen("out.txt","w",stdout);
    #endif
    int n;
    scanf("%d",&n);
    b[0]=0;
    for(int i=1;i<=n;i++){
        ll x=1;
        a[i]=b[i-1]+1;
        while(x*2<i){
            if(x&i){
                a[i]=(a[i]-b[x*2-1]+b[x-1]+mod)%mod;
            }
            x*=2;
        }
        b[i]=(b[i-1]+a[i])%mod;
    }
    printf("%lld",b[n]);
    return 0;
}

```

