# 2020-2021 Winter Petrozavodsk Camp, Day 9 Contest (XXI Open Cup, Grand Prix of Suwon)

## AC 2/12  目标：3/12

[传送门](https://codeforces.com/gym/102979)

**赛后补题目标** ：

F 

---

### WA:

J ： 想错了

F：没开longlong，中间有个判定条件爆了

---



---

## 需学习的知识点：



----

# ~



## I. Integer Array Shuffle

**题意：** 每次将A数组的最左边或者最右边的元素加入到B，问至少多少次之后数组变成有序的数列



```cpp
#include<bits/stdc++.h>
using namespace std;
typedef long long ll;

int a[300005];
int main(){
    #ifdef ONLINE_JUDGE
        //std::ios::sync_with_stdio(false);
    #else
        //freopen("in.txt","r",stdin);
        //freopen("out.txt","w",stdout);
    #endif
    int n;
    scanf("%d",&n);
    for(int i=0;i<n;i++){
        scanf("%d",&a[i]);
    }
    int sum=1;
    for(int i=1;i<n;i++){
        if(sum%2==1&&a[i]<a[i-1])sum++;
        if(sum%2==0&&a[i]>a[i-1])sum++;
    }
    int res=0,x=1;
    while(x<sum){
        res++,x*=2;
    }
    printf("%d",res);
    return 0;
}
```



## J. Junkyeom's Contest

**题意：** 在给出的数列里面找到七个数，满足$p_1 \geq p_2 \geq p_3 \geq p_4 \geq p_5 \geq p_6 \geq p_7$ 并且 $p_1 < p_2 + p_ 3 < p_4 + p_5 + p_ 6 + p_7$ ，问和最大的方案，

**思路：** 刚开始想错了。这七个数，最多分为$p_1$，$p_2$ ，$p_3- p_7$，三个部分，枚举$p_2$，然后根据给出的条件2维护一个后四个数的优先队列，然后两边二分去找

```cpp
#include<bits/stdc++.h>
using namespace std;
typedef long long ll;

ll a[500005],b[500005],c[500005],d[500005];
bool cmp(ll a,ll b){
    return a>b;
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
    for(int i=0;i<n;i++)scanf("%lld",&a[i]);
    sort(a,a+n,cmp);
    int cnt=0;
    ll res=0;
    for(int i=n-6;i>=1;i--){
        d[cnt]=i+1;
        b[cnt]=a[i+2]+a[i+3]+a[i+4]+a[i+5]-a[i+1];
        c[cnt]=a[i+1]+a[i+2]+a[i+3]+a[i+4]+a[i+5];
        while(cnt&&b[cnt-1]<b[cnt]){
            d[cnt-1]=d[cnt];
            b[cnt-1]=b[cnt];
            c[cnt-1]=c[cnt];
            cnt--;
        }
        cnt++;
        ll sum=0;
        int l=0,r=cnt-1;
        if(b[l]-a[i]<=0)continue;
        while(l<r){
            int mid=(l+r+1)/2;
            if(b[mid]-a[i]>0)l=mid;
            else r=mid-1;
        }
        ll p3=d[l];
        sum=c[l]+a[i];
        if(a[p3]+a[i]-a[i-1]<=0)continue;
        l=0,r=i-1;
        while(l<r){
            int mid=(l+r)/2;
            if(a[p3]+a[i]-a[mid]>0)r=mid;
            else l=mid+1;
        }
        sum+=a[l];
        res=max(res,sum);
    }
    if(res==0)res=-1;
    printf("%lld",res);
    return 0;
}
```



## F. Find the XOR

**题意：** 给一个n个点m条带权边的连通图（有环有自环有重边），定义$d(u,v)$为u到v的一条路径中每条边xor的和的最大值，给q次询问，每次询问一个$l, r$，求**所有**$l \leq i < j \leq r$ 的$d(i,j)$ 的最大值

**思路：**// 

1.从u到v的每一条路径上的边异或和 都可以用 u到v的任意一条路径 和 图中若干环 的异或和 表示

2.u到v的任意一条路径的异或和可以用从根到u的路径（P(root, u) ）异或上从根到v的路径得到，表示为$P(u, v) = P(root, u) \bigoplus P(root, v)$

3.L到R区间内的所有$d(i,j)$的异或和可以表示为 $P(root, i) \bigoplus P(root, j)$ 的异或和，根据异或的性质，为0或从P(root, L) 异或到 P(root, R) 的区间异或和

4.设mx（x）为x与线性基异或的最大值
mn（x）为x与线性基异或的最小值
$mx(x1)\bigoplus mx(x2)=mn(x1\bigoplus x2)$
$mx(x1)\bigoplus mx(x2)\bigoplus mx(x3)=mx(x1\bigoplus x2\bigoplus x3)$

5.线性基可以将很大的一个数的集合等价成  最高有效二进制位  个数的集合，并且两个集合中取若干个数能生成的异或和的集合是相同的，线性基从高位向低位枚举可以做到取一个数和这个集合的异或的最大或最小

6.环可以由dfs序得到

```cpp
#include <bits/stdc++.h>

using namespace std;
typedef long long ll;

const int N = 1 << 17;

const int maxbit = 30;  

struct L_B {
    int lba[maxbit];
    bool empty;

    L_B() {
        memset(lba, 0, sizeof(lba));
        empty = 1;
    }

    void insert(int val) { 
        empty = 0;
        for (int i = maxbit - 1; i >= 0; --i)
            if (val & (1 << i)) { //
                if (!lba[i]) {
                    lba[i] = val;
                    break;
                }
                val ^= lba[i];
            }
    }

    int query1(int ans) {
        for (int i = maxbit - 1; i >= 0; i--) {
            if (ans < (ans ^ lba[i])) ans ^= lba[i];
        }
        return ans;
    }

    int query2(int ans) {
        for (int i = maxbit - 1; i >= 0; i--) {
            if (ans > (ans ^ lba[i])) ans ^= lba[i];
        }
        return ans;
    }

} base;

struct edge {
    int to, v;
    int num;
};
int  dfn[N], dfs_clock;
vector<edge> g[N];
int val[N];

void tarjan(int u, int fa) {
     dfn[u] = ++dfs_clock;
    for (int i = 0; i < g[u].size(); i++) {
        int v = g[u][i].to;
        if (!dfn[v]) {
            val[v] = g[u][i].v ^ val[u];
            tarjan(v, g[u][i].num);
        } else if (g[u][i].num != fa) {
            base.insert(val[v] ^ val[u] ^ g[u][i].v);
        }
    }
}

int xr[N];

int main() {
    ios::sync_with_stdio(false), cin.tie(nullptr), cout.tie(nullptr);
    
    int n, m, q;
    cin >> n >> m >> q;
    for (int i = 1; i <= m; i++) {
        int a, b, c;
        cin >> a >> b >> c;
        g[a].push_back({b, c, i});
        g[b].push_back({a, c, i});
    }
    tarjan(1, 0);
    for (int i = 1; i <= n; i++) xr[i] = xr[i - 1] ^ val[i];

    while (q--) {
        ll a, b;
        cin >> a >> b;
        ll x = xr[b] ^ xr[a - 1];
        if ((b - a + 1) % 2 == 1) x = 0;
        if ((b - a + 1) * (b - a) / 2 % 2 == 1) {
            cout << base.query1(x) << endl;
        } else {
            cout << base.query2(x) << endl;
        }
    }
}
```