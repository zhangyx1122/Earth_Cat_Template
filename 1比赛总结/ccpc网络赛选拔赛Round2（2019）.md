# ccpc网络赛选拔赛Round2（2019）

## AC 4/11   目标：5/11（

**赛后补题目标 ：**

**B - array （hdu6703）**(已)

**E - huntian oy（hdu6706)**  (-)

### WA：3

F ：操作遍历反了

H :   //看不出然后换人重写了

H :   没开long long

### TLE: 1

H :  //看不出然后换人重写了

---

## 需学习的知识点：

**主席树** （B）

**莫比乌斯反演** （E）

**杜教筛** （E）

---

# ~

##A - ^ & ^ hdu 6702

**题意** ：给出A，B，求使$(A xor C) $  &  $(BxorC)$ 最小的C，如果当表达式为0时C为0，则输出1

```c++
#include<bits/stdc++.h>
using namespace std;
int main(){
    ios::sync_with_stdio(0),cin.tie(0),cout.tie(0);
    int t;
    cin>>t;
    while(t--){
        long long a,b;
        cin>>a>>b;
        if((a&b)==0) cout<<1<<endl;
        else  cout<<(a&b)<<endl;
    }
}
```

----

## G - Windows Of CCPC hdu 6708

**题意** :构造`ccpc`矩阵

递推

```c++
#include<bits/stdc++.h>
using namespace std;
char g[1100][1100];

void dfs(int k){
    if(k==1) return;
    dfs(k-1);
    for(int i=1;i<=(1<<(k-1));i++){
        for(int j=1;j<=(1<<(k-1));j++){
            g[(1<<(k-1))+i][(1<<(k-1))+j]=g[i][j];
            g[i][(1<<(k-1))+j]=g[i][j];
            g[(1<<(k-1))+i][j]=(g[i][j]=='C'?'P':'C');
        }

    }
}

int main(){
    ios::sync_with_stdio(0),cin.tie(0),cout.tie(0);
    int t;
    cin>>t;
    while(t--){
        int n;
        cin>>n;
        g[1][1]='C';
        g[1][2]='C';
        g[2][1]='P';
        g[2][2]='C';
        dfs(n);
        for(int i=1;i<=(1<<n);i++){
            for(int j=1;j<=(1<<n);j++){
                cout<<g[i][j];
            }
            cout<<endl;
        }
    }
}
```

---

## F - Shuffle Card hdu 6707

**题意** ：有n张有序排列的牌， 每次选一个数字的牌放到最前面，问m次操作后的牌的顺序

```c++
#include<bits/stdc++.h>
using namespace std;
vector<int> v;
int num[200000];
int a[200000];
int main(){
    ios::sync_with_stdio(0),cin.tie(0),cout.tie(0);
    int n,m;
    cin>>n>>m;
    for(int i=1;i<=n;i++) cin>>num[i];
    set<int> s;
    for(int i=1;i<=m;i++){
        cin>>a[i];
    }
    for(int i=m;i>0;i--){
        if(s.count(a[i])) continue;
        v.push_back((a[i]));
        s.insert((a[i]));
    }
    for(auto x:v) cout<<x<<' ';
    for(int i=1;i<=n;i++){
        if(s.count(num[i])==0){
            cout<<num[i]<<' ';
        }
    }
}
```

---

## H - Fishing Master hdu 6709

**题意** ：有n条鱼，每条鱼的钓鱼时间为k，煮鱼时间为$a_i$， 问将n条鱼都钓上来煮熟所需要的最短时间

贪心+优先队列

　假设你当前煮的鱼需要花费 t 时间，钓鱼需要花费 k 时间；你可以在这 t 时间内钓 $t/k$条鱼上来，在钓鱼的时间，锅处于煮鱼状态；但是剩下的 $t%k$ 时间不足以再钓一条上来；此时，你就有两个决策可以选择：决策1：去钓下一条鱼；决策2：等待$ t%k$ 时间往锅中放入下一条鱼；当然，选择 决策2 的前题是你得有鱼可煮；如果你手中有鱼的话，肯定要选择 决策2，因为等待的这 t%k 时间是必须的；

于是优先选择时间长的鱼

```c++
#include<bits/stdc++.h>
using namespace std;
typedef long long ll;

int main(){
    ios::sync_with_stdio(0),cin.tie(0),cout.tie(0);
    long long t;
    cin>>t;
    while(t--){
         long long  n;
         long long  k;
        cin>>n>>k;
        priority_queue< long long > que;
         long long ans=0;
        ans=n*k;
        for(int i=1;i<=n;i++) {
            long long a;
            cin>>a;
            que.push(a);
            ans+=a;
        }
        for(int i=1;i<n;i++){
             long long  p=que.top();
            que.pop();
            if(p>k){
                ans-=k;
                que.push(p-k);
            }else {
            ans-=min(k,p);
            }
        }
        cout<<ans<<endl;
    }
}
```

---

## B - array hdu6703

**题意** ：有n个数，$a_1$ 到 $a_n$ ，保证不重复且$\in [1,n]$，有2种操作，第一种$a_{pos} \to  a_{pos} + 10^7$, 第二种求与前r个数都不同且至少为k的最小的数

第一种操作$（1，t_1）$ 其中$pos = t_1 \bigoplus LASTANS$

第二种操作 $(2,t_2,t_3)$  其中$r = t_2 \bigoplus LASTANS, k = t_3 \bigoplus LASTANS$



因为 $1 \leq n \leq 10^5$,  $ans \in [1, n+ 1]$， 操作1就相当于删数

[sol1 主席树+set   sol2 权值线段树](https://blog.csdn.net/qq_40482358/article/details/100053523)

**主席树+set**

对于每个删除的元素，将其插入到set中，查询的时候，查询[r+1,n+1]中第一个>=k的元素，再和set中第一个>=k的元素相比较，选择最小的
```c++
#include <bits/stdc++.h>
using namespace std;
typedef long long ll;

const int N = 1 << 17;
int tree[N << 5], lc[N << 5], rc[N << 5];
int root[N], tot;
int shu[N];

int build(int l, int r) {
    int num = ++tot;
    tree[num] = l;
    if (l == r) return num;
    int mid = (l + r) >> 1;
    lc[num] = build(l, mid);
    rc[num] = build(mid + 1, r);
    return num;
}

int update(int k, int l, int r, int x) {
    int num = ++tot;
    lc[num] = lc[k], rc[num] = rc[k];
    if (l < r) {
        int mid = (l + r) >> 1;
        if (x <= mid)
            lc[num] = update(lc[k], l, mid, x);
        else
            rc[num] = update(rc[k], mid + 1, r, x);
        tree[num] = min(tree[lc[num]], tree[rc[num]]);
    } else {
        tree[num] = INT_MAX;
    }
    return num;
}

int query(int k, int l, int r, int ql, int qr) {
    if (ql <= l && r <= qr) return tree[k];
    if (qr < l || r < ql) return INT_MAX;
    int mid = (l + r) >> 1;
    return min(query(lc[k], l, mid, ql, qr), query(rc[k], mid + 1, r, ql, qr));
}

void init() { tot = 0; }

int main() {
    int t;
    scanf("%d", &t);
    while (t--) {
        init();
        int lastans = 0;
        int n, m;
        scanf("%d%d", &n, &m);
        root[0] = build(1, N);
        set<int> s;
        for (int i = 1; i <= n; i++) {
            scanf("%d", &shu[i]);
            root[i] = update(root[i - 1], 1, N, shu[i]);
        }
        for (int i = 1; i <= m; i++) {
            int a;

            scanf("%d", &a);
            if (a == 1) {
                scanf("%d", &a);
                a ^= lastans;
                s.insert(shu[a]);
            } else {
                int b;
                scanf("%d%d", &a, &b);
                a ^= lastans;
                b ^= lastans;
                int ans = query(root[a], 1, N, b, N);
                if (ans > n) ans = n + 1;
                set<int>::iterator it = s.lower_bound(b);
                if (it != s.end()) ans = min(ans, *it);
                printf("%d\n", ans);
                lastans = ans;
            }
        }
    }
}
```


---

## E - huntian oy hdu 6706

**题意** ： 求 $f(n,a,b) = \sum_{i = 1}^{n}\sum_{j= 1}^{i}{gcd(i^a-j^a,i^b-j^b)[gcd(i,j)= 1]} % (10^9+7)$



**莫比乌斯反演** ：

[OI wiki 莫比乌斯反演](https://oi-wiki.org/math/mobius/)

**dirichlet卷积**

定义两个数论函数$f,g$的**dirichlet卷积**为$(f*g)(n)=\sum_{d|n}f(d)g(\frac{n}{d})$

**莫比乌斯函数**
$$
\mu(x)=\begin{cases}1, &  n =1    \\0, &     \exists d >1 : d^2 | n \\(-1)^{w(n)}, & otherwise (w(n)表示n的本质不同的因数个数)\end{cases}
$$
**性质** ：
$$
\sum_{d|n} \mu(d) = \begin{cases}
1, & n = 1\\
0, & n \neq 1 \\
\end{cases}
$$
线性筛莫比乌斯函数，（线性筛基本可以求所有的积性函数   ~~但不会~~

**莫比乌斯反演**

公式：

设$f(n), g(n)$为两个数论函数，

若有$f(n) = \sum_{d|n}g(d)$，则有$f(n) = \sum_{d|n}{\mu(d)f(\frac{n}{d})}$

若有$f(n) = \sum_{n|d}{g(d)}$，则有$g(n)= \sum_{n|d}{\mu(\frac{d}{n})f(d)}$

**杜教筛**

[OI wiki 杜教筛](https://oi-wiki.org/math/du/)

**积性函数** ：对于所有**互质**的a和b，总有$f(ab)=f(a)*f(b)$

常见的有：

+ **约数个数**   $d(x) = \sum_{i|n}1$         

+ **约数的和**   $\sigma(x) = \sum_{i|n} i$

+ **欧拉函数**   $\psi(x) = \sum_{i = 1} ^{x} 1 [gcd(x,i) = 1]$

+ $$
  莫比乌斯函数  \mu(x)=
  \begin{cases}
  1, &  n =1    \\
  0, &     \exists d >1 : d^2 | n \\
  (-1)^{w(n)}, & otherwise (w(n)表示n的本质不同的因数个数)
  \end{cases}
  $$

<img src="sdfgaw e2342.png" style="zoom: 70%;" />

**积性函数的性质** ：若$f(x),g(x)$均为积性函数，则$h(x)$也为积性函数

+ $h(x) = f(x^p)$
+ $h(x)=f^p(x)$
+ $h(x) = f(x)g(x)$
+ $h(x)=\sum_{d|x}f(d)g(\frac{x}{d})$

杜教筛被用来处理**数论函数的前缀和**问题。对于求解一个前缀和，杜教筛可以在低于线性时间的复杂度内求解

$O(n^{\frac{2}{3}})$ (?)

[P3768 简单的数学题](https://www.luogu.com.cn/problem/P3768)

应用： 莫比乌斯函数前缀和，欧拉函数前缀和





[推导  题解](https://blog.csdn.net/tirion_chenrui/article/details/100112163)



$x^n - y^n = (x - y)(x^{n - 1} + x^{n - 2}y + \dots + y^{n - 1})$

$gcd(i^a-j^a,i^b-j^b)[gcd(i,j)= 1] = i - j$$\sum_{i = 1}^{n}\sum_{j = 1}^{i}{j[gcd(i,j) = 1]} = \sum_{i = 1}^{n}{(\frac{i\psi(i)}{2}+[i = 1])}$

于是问题变成如何求$\sum_{i = 1}^{n}i\psi(i), n \leq 1e9$





