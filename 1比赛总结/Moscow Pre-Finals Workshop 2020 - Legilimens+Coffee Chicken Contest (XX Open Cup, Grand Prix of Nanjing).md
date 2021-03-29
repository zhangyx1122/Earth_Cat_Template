# Moscow Pre-Finals Workshop 2020 - Legilimens+Coffee Chicken Contest (XX Open Cup, Grand Prix of Nanjing)

## AC 3/13  目标：4/13

[Moscow Pre-Finals Workshop 2020 - Legilimens+Coffee Chicken Contest (XX Open Cup, Grand Prix of Nanjing)](https://codeforces.ml/gym/102994)

---

**赛后补题目标** ：（没有题解，咕咕（不是



---

### WA：

L : 在分类讨论情况的时候前面那种把后面一起带过了，交换顺序后通过



---

感觉实力还是有点差距（... 要加油鸭

---

## 需学习的知识点：



----

# ~

## J

**题意：** 斐波那契数列第n行有多少个奇数

```cpp
#include<iostream>
using namespace std;
typedef long long ll;

ll a[1005][1005];

 ll popcount(ll x){
     ll ans=0;
     while(x){
         x-=x&-x;
         ans++;
     }
     return ans;
 }

 ll qpow(ll a,ll b){
     ll ans=1;
     while(b){
         if(b&1) ans=ans*a;
         b>>=1;
         a=a*a;
     }
     return ans;
 }

int main(){
    ios::sync_with_stdio(0),cin.tie(0),cout.tie(0);
    ll t;
    cin>>t;
    while(t--){
        ll n;
        cin>>n;
        n--;
        cout<<qpow(2,popcount(n))<<endl;
    }
}
```





## L

**题意：** 给两个矩形的位置，问把二维平面分成几部分

```cpp
#include<iostream>

using namespace std;
typedef long long ll;

int main() {
    ios::sync_with_stdio(0), cin.tie(0), cout.tie(0);
    int t;
    cin >> t;
    while (t--) {
        int a[4], b[4];
        for (int i = 0; i < 4; i++) cin >> a[i];
        for (int i = 0; i < 4; i++) cin >> b[i];
        int cnt = 0;
        for (int i = 0; i < 4; i++) if (a[i] == b[i]) cnt++;
        if (a[0] == b[0] && a[1] == b[1] && a[2] == b[2] && a[3] == b[3]) {
            cout << 2 << endl;
        } else if (a[2] <= b[0] || a[3] <= b[1] || b[2] <= a[0] || b[3] <= a[1]) {
            cout << 3 << endl;
        }  else if (a[0] <= b[0] && b[2] <= a[2] && b[1] <= a[1] && a[3] <= b[3]
                   || b[0] <= a[0] && a[2] <= b[2] && a[1] <= b[1] && b[3] <= a[3]) {
            cout << 6 - cnt << endl;
        } else if (a[2] <= b[2] && a[3] <= b[3] && b[0] <= a[0] && b[1] <= a[1]
                   || b[2] <= a[2] && b[3] <= a[3] && a[0] <= b[0] && a[1] <= b[1]) {
            cout << 3 << endl;
        }else {
            cout << 4 << endl;
        }
    }
}
```



## A

**题意：** A有n对（x, y）数，可以选择x或y异或到答案上，B有m对数，同样可以选择x或y异或到答案上，A想使答案尽可能地大，B想使答案尽可能的小，问答案会是多少

把xy异或起来得到的数放入线性基，（这时候会有一个基准的数，直接异或到答案里），得到A的线性基和B的线性基，对于每一位AB线性基的情况进行讨论



```cpp
#include<iostream>
#include<cstring>

using namespace std;
typedef long long ll;
const ll maxbit = 62;
const ll N = 1 << 20;

struct L_B {
    ll lba[maxbit];

    L_B() { memset(lba, 0, sizeof lba); }

    void insert(ll val) {
        for (ll i = maxbit - 1; i >= 0; --i) {
            if (val & (1ll << i)) {
                if (!lba[i]) {
                    lba[i] = val;
                    break;
                }
                val ^= lba[i];
            }
        }
    }

    ll query(ll val) {
        for (ll i = maxbit - 1; i >= 0; --i) {
            if (val & (1ll << i)) {
                continue;
            }
            val ^= lba[i];
        }
        return val;
    }
};

int main() {
    ios::sync_with_stdio(0), cin.tie(0), cout.tie(0);
    ll t;
    cin >> t;
    while (t--) {
        ll n, m;
        cin >> n >> m;
        L_B la, lb;
        ll init = 0;
        for (ll i = 1; i <= n; i++) {
            ll a, b;
            cin >> a >> b;
            init ^= a;
            la.insert(a ^ b);
        }
        for (ll i = 1; i <= m; i++) {
            ll a, b;
            cin >> a >> b;
            init ^= a;
            lb.insert(a ^ b);
        }
        for (ll i = maxbit - 1; i >= 0; --i) {
            if (init & (1ll << i)) {
                if (la.lba[i] == 0 && lb.lba[i] == 0) {
                    continue;
                } else if (la.lba[i] && !lb.lba[i]) {
                    continue;
                } else if (!la.lba[i] && lb.lba[i]) {
                    init ^= lb.lba[i];
                } else {
                    init ^= la.lba[i];
                    ll tmp = la.lba[i] ^lb.lba[i];
                    la.lba[i] = 0;
                    la.insert(tmp);
                }
            } else {
                if (la.lba[i] == 0 && lb.lba[i] == 0) {
                    continue;
                } else if (la.lba[i] && !lb.lba[i]) {
                    init ^= la.lba[i];
                    la.lba[i] = 0;
                } else if (!la.lba[i] && lb.lba[i]) {
                    continue;
                } else {
                    ll tmp = la.lba[i] ^lb.lba[i];
                    la.lba[i] = 0;
                    la.insert(tmp);
                }
            }
        }
        cout << la.query(init) << endl;
    }
}
```



## I

**题意：** 有n个粉丝m个队伍，每个粉丝至少粉一个队伍但是没有人是所有队伍的粉丝，设一个队伍 i 的粉丝集合为$T_i$ ， $\forall i, j , \exist k,  Ti \bigcap Tj = Tk $ ,  $\forall i, j , \exist k,  Ti \bigcup Tj = Tk $

问方案数，

注意到队伍数非常小，

分析发现根据传递性，必定会有一个队伍所有人都是他的粉丝，没有两个队伍的粉丝集合相同，并且由于没有一个粉丝是所有队伍的集合，必定会有一个队伍没有粉丝

中间算的时候要用到容斥原理，想到后面就有点乱（

（隔壁队根据给的样例解系数方程真的给我看傻了)
