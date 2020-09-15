# ccpc网络赛选拔赛Round3（2016）

## AC 5/11   目标：6/11

**赛后补题目标 ：H - Special Tetrahedron（hdu5839）** (已)

### WA：3

D ：有种情况忘记乘2了

B :   $a_i$范围$10^{18}$, 没有开 long long

B :   高斯消元写炸了

### TLE: 5

A ：求幂次的地方用了快速幂，改线性递推

B ：// 高斯消元写炸了

B ：// WA 报 T

C ：//

C ：没有考虑菊花图



## 需学习的知识点：

**高斯消元** （B）

**四点共面判断**（？ （H）

---

# ~

## A - A water problem hdu5832

**题意** ：大数取模

```c++
#include<bits/stdc++.h>
using namespace std;
typedef long long ll;

int main(){
    ios::sync_with_stdio(false),cin.tie(nullptr),cout.tie(nullptr);
    int t=0;
    ll n;
        string s;
    while(cin>>s){
        ++t;
        ll m1=0,m2=0;
        reverse(s.begin(),s.end());
        int n1=1,n2=1;
        for(int i=0;i<s.length();i++){
            m1=(m1+n1*(s[i]-'0'))%73;
            m2=(m2+n2*(s[i]-'0'))%137;
            n1=(n1*10)%73;
            n2=(n2*10)%137;
        }
        cout<<"Case #"<<t<<": ";
        m1=(m1+1)%73;
        m2=(m2+1)%137;
        if(m1==1&&m2==1){
            cout<<"YES"<<endl;
        }else {
            cout<<"NO"<<endl;
        }
    }
}
```

---

 ## K - Lweb and String hdu 5842

**题意** ： 假的LIS，每个字母代表的数字可以自定义，于是直接判断一共出现过多少个字母就可

```c++
#include<bits/stdc++.h>
using namespace std;
typedef long long ll;

int main(){
    ios::sync_with_stdio(false),cin.tie(nullptr),cout.tie(nullptr);
    int t;
    cin>>t;
    for(int i=1;i<=t;i++){
        set<char> s;
        string str;
        cin>>str;
        for(auto x:str){
            s.insert(x);
        }
        cout<<"Case #"<<i<<": "<<s.size()<<endl;
    }
}
```

---

## D - Danganronpa hdu 5835

**题意** ：有n种礼物，第i个有$a_i$个。现在要把他们分给一排同学，每个人发普通礼物和神秘礼物各一个，要求相邻两人的普通礼物不能为同一种，而神秘礼物随意。问最多能发给多少人。

**思路** ：贪心，找数量最多的礼物，如果该礼物数小于等于总和的一半，将该礼物作为普通礼物时插空，总存在使$ 总礼物个数 / 2 $的学生的方法，否则先分配普通礼物：我们把其他礼物按间隔为1摆放，然后让绝对大的礼物去插空。然后多余的去填补神秘礼物，如果神秘礼物也能填补完，那么总数就为$其他礼物数*2+1 $；如果填补不完，就从普通礼物后面截取然后往神秘礼物填，直到神秘和普通数量相同，总数为$(其他礼物+总数最大的礼物个数)/2$

```c++
#include<bits/stdc++.h>
using namespace std;
typedef long long ll;

int a[1000005];
int main(){
    ios::sync_with_stdio(false),cin.tie(nullptr),cout.tie(nullptr);
    int t;
    cin>>t;
    for(int j=1;j<=t;j++){
        int n;
        cin>>n;
        int sum=0;
        int mmax=0;
        for(int i=1;i<=n;i++){
            cin>>a[i];
            sum+=a[i];
            mmax=max(mmax,a[i]);
        }
        int now=sum/2;
        int ans=0;
        if(mmax>sum-mmax){
            ans=min(sum/2,(sum-mmax)*2);
        }else{
            ans=sum/2;
        }
        cout<<"Case #"<<j<<": "<<ans<<endl;
    }
}
```

---

## B - Zhu and 772002 hdu5833

**题意** ： 给出一堆数$a_i$ ，保证每个数的质因数的值小于2000， 求使选取若干个数（可以为1）的乘积为完全平方数的个数

**思路** ：高斯消元求自由变元的个数

筛出2000内的质数，对每个数$a_i$，分解放入矩阵$g[i, j]$ ，（$g[i, j]$表示第i个数含有第j个素数的个数），因为要构成完全平方数，mod2后为01矩阵，异或求解，，最后减去全为0的情况



题目没告诉T的范围，然后报错tle的时候还以为卡常了，就疯狂优化~~一些没啥用的地方~~吃罚时..

中间还自己造数据测了一下本地时间，但是这个时间也有点迷...

后来才发现是有个判定写错了，（但是写错并不会增加循环次数等 ）不明白为什么报T不是WA

下次还是要先看看代码是否存在逻辑性错误（？

```c++
#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
const ll mod=1000000007;

bool p[2005];
int tot;
ll prime[2005];
bool g[305][305];
bool vis[400];

inline ll qpow(ll a,ll b){
    ll ans=1;
    while(b){
        if(b&1) ans=(ans*a)%mod;
        a=(a*a)%mod;
        b>>=1;
    }
    return ans;
}

int main(){
    ios::sync_with_stdio(0),cin.tie(0),cout.tie(0);
    p[1]=1;
    for(int i=2;i<2001;i++){
        if(p[i]==0){
            prime[++tot]=i;
            for(int j=i+i;j<2001;j+=i){
                p[j]=1;
            }
        }
    }
    int t;
    cin>>t;
    int s;
    for(int o=1;o<=t;o++){
        int n;
        cin>>n;
        for(int i=1;i<=n;i++){
            memset(g[i],0,sizeof g[i]);
        }
        memset(vis,0,sizeof vis);
        for(int i=1;i<=n;i++){
            ll a;
            cin>>a;
            for(int j=1;j<=tot;j++){
                int cnt=0;
                while(a%prime[j]==0){
                    cnt++;
                    a/=prime[j];
                }
                g[i][j]=cnt%2;
            }
        }
        s=clock();
        ll zhi=0;
        for(int j=1;j<=tot;j++){
            int p=0;
            for(int i=1;i<=n;i++) {
                if(g[i][j]&&!vis[i]){
                    vis[i]=1;
                    p=i;
                    break;
                }
            }
            if(!p) continue;
            zhi++;
            for(int i=1;i<=n;i++){
                if(g[i][j]&&i!=p){
                    for(int k=j;k<=tot;k++){
                        g[i][k]^=g[p][k];
                    }
                }
            }
        }
        cout<<"Case #"<<o<<":\n";
        ll ans=(qpow(2,n-zhi)+mod-1)%mod;
        cout<<ans<<endl;
        //cout<<clock()-s<<endl;
    }
}
```

---

## C - Magic boy Bi Luo with his excited tree hdu 5834

**题意** ： 给出一张无向图，每个点每条边都有一个权值，求从每个点出发可以获得的最多的点的权值（边可以经过多次，每经过一次就要支付一次费用）

**思路** : 树形DP，以1为根

对于每个点维护往下走并且回来的最大值，往下走不回来的最大值，往上走并且回来的最大值，往上走不回来的最大值，

```c++
#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
struct edge {
    int to, c;
    edge() {}
    edge(int a, int b) : to(a), c(b) {}
};

vector<edge> g[200000];
int v[200000];
int fa[200000];
int len[200000];
int w[200000][4];
int to[200000];

void dfs1(int u, int f, int l) {
    fa[u] = f;
    len[u] = l;
    w[u][0] = w[u][1] = v[u];
    for (auto x : g[u]) {
        if (x.to == f) continue;
        dfs1(x.to, u, x.c);
        if (w[x.to][0] - x.c * 2 > 0) w[u][0] += w[x.to][0] - x.c * 2;
    }
    to[u] = 0;
    for (auto x : g[u]) {
        if (x.to == f) continue;
        if (w[x.to][0] - x.c * 2 > 0) {
            int kk = w[u][0] - (w[x.to][0] - x.c * 2) + w[x.to][1] - x.c;
            if (kk > w[u][1]) {
                w[u][1] = kk;
                to[u] = x.to;
            }
        } else {
            int kk = w[u][0] + w[x.to][1] - x.c;
            if (w[u][1] < kk) {
                w[u][1] = kk;
                to[u] = x.to;
            }
        }
    }
}

void dfs2(int u, int f) {
    if (f == -1)
        w[u][2] = w[u][3] = 0;
    else {
        int now = w[f][0];
        if (w[f][2] > 0) now += w[f][2];
        if (w[u][0] - len[u] * 2 > 0) now -= w[u][0] - len[u] * 2;

        w[u][2] = max(now - len[u] * 2, 0);
        int nnow = now;
        if (to[f] == u) {
            for (auto x : g[f]) {
                if (x.to == u) {
                    continue;
                }
                if (x.to == fa[f]) {
                    if (w[f][2] > 0) {
                        nnow = max(nnow, now - w[f][2] + w[f][3]);
                    } else {
                        nnow = max(nnow, now + w[f][3]);
                    }
                } else {
                    if (w[x.to][0] - x.c * 2 > 0) {
                        nnow = max(nnow, now - (w[x.to][0] - x.c * 2) +
                                             w[x.to][1] - x.c);
                    } else {
                        nnow = max(nnow, now + w[x.to][1] - x.c);
                    }
                }
            }
        } else {
            nnow = max(w[f][0] + w[f][3], w[f][1] + w[f][2]);
            if (w[u][0] - len[u] * 2 > 0) nnow -= w[u][0] - len[u] * 2;
        }
        w[u][3] = max(nnow - len[u], 0);
    }
    for (auto x : g[u]) {
        if (x.to == f) continue;
        dfs2(x.to, u);
    }
}

int main() {
    ios::sync_with_stdio(0), cin.tie(0), cout.tie(0);
    
    int t;
    cin >> t;
    for (int y = 1; y <= t; y++) {
        int n;
        cin >> n;
        for (int i = 1; i <= n; i++) g[i].clear();
        for (int i = 1; i <= n; i++) cin >> v[i];
        for (int i = 1; i < n; i++) {
            int a, b, c;
            cin >> a >> b >> c;
            g[a].push_back(edge(b, c));
            g[b].push_back(edge(a, c));
        }
        dfs1(1, -1, 0);
        dfs2(1, -1);
        cout << "Case #" << y << ":" << endl;
        for (int i = 1; i <= n; i++) {
            int ans = max(w[i][0] + w[i][3], w[i][1] + w[i][2]);
            cout << ans << endl;
        }
    }
}

```

---

## H - Special Tetrahedron hdu5839

**题意** ： 给出三维平面的n个点，问有多少个满足以下两个条件的四面体：至少四条边相等； 若正好四条边相等，则剩余两条边不能相邻。

**思路** ：每次暴力枚举两个点A, B，再暴力枚举剩余点 C，看C是否在A B的中垂面上，（C到AB距离是否相等），再暴力枚举在中垂面上的点，判断是否共面，是否为正四面体，是否为四面体

正四面体被计算了6次，四面体被计算了2次，分开记录最后除一下就好

思路没问题，可能四点共面写炸了（，中间判断中垂面上的点的时候也不用排序，直接暴力复杂度也不会爆

**四点共面**

[判断四点共面](https://blog.csdn.net/a892573486/article/details/79143669)

 已知四个点坐标判断是否共面
  可以用行列式来判断
  用四个点求出三个向量分别为$(x_1,y_1,z_1),(x_2,y_2,z_2),(x_3,y_3,z_3)$
$$
判断行列式\begin{vmatrix}
x_1 & x_2 & x_3 \\
y_1 & y_2 & y_3 \\
z_1 & z_2 & z_3 \\
\end{vmatrix}是否为零， 若为零则四点共面
$$

```c++
#include <bits/stdc++.h>
using namespace std;

const int maxn = 210;

struct node{
	int x, y, z;
}po[maxn];

int dis(node a, node b){
	return (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) 
				+ (a.z - b.z) * (a.z - b.z); 
}

bool gongmian(node a, node b, node c, node d){
	int x1, x2, x3, y1, y2, y3, z1, z2, z3;
	x1 = b.x - a.x; y1 = b.y - a.y; z1 = b.z - a.z;
	x2 = c.x - a.x; y2 = c.y - a.y; z2 = c.z - a.z;
	x3 = d.x - a.x; y3 = d.y - a.y; z3 = d.z - a.z;
	if((x1*y2*z3)+(x2*y3*z1)+(x3*y1*z2)-(x3*y2*z1)-(y3*z2*x1)-(z3*x2*y1)==0) return 1;
	return 0;
} 

bool zhengsimianti(node a, node b, node c, node d){
	if(dis(a, b) != dis(a, c)) return 0;
	if(dis(a, b) != dis(a, d)) return 0;
	if(dis(a, b) != dis(b, c)) return 0;
	if(dis(a, b) != dis(b, d)) return 0;
	if(dis(a, b) != dis(c, d)) return 0;
	return 1;
}

bool simianti(node a, node b, node c, node d){
	return (dis(a, c) == dis(a, d));
}

int main(){
	ios::sync_with_stdio(0); cin.tie(0); cout.tie(0);
	int T;
	cin >> T;
	int cas = 0;
	while(T --){
		++ cas;
		int n;
		cin >> n;
		for(int i = 1; i <= n; ++ i){
			cin >> po[i].x >> po[i].y >> po[i].z;
		}
		
		vector<node> vec;
		int ans1 = 0, ans2 = 0;
		
		for(int a = 1; a <= n; ++ a){
			for(int b = a + 1; b <= n; ++ b){
				vec.clear();
				for(int c = 1; c <= n; ++ c){
					if(a == c || b == c) continue;
					if(dis(po[a], po[c]) == dis(po[b], po[c])) 
						vec.push_back(po[c]); 
				}
				
				int siz = vec.size();
				for(int i = 0; i < siz; ++ i){
					for(int j = i + 1; j < siz; ++ j){
						if(gongmian(po[a], po[b], vec[i], vec[j])) continue;
						if(zhengsimianti(po[a], po[b], vec[i], vec[j])){
							ans2 ++;
							continue;
						} 
						if(simianti(po[a], po[b], vec[i], vec[j])) ans1 ++;
					} 
				}
				
			}
		}
		cout << "Case #" << cas << ": "; 
		cout << ans1 / 2 + ans2 / 6 << endl;
	}
}
```

---