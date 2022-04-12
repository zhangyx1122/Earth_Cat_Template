#include<cstdio>

using namespace std;
typedef long long ll;

ll n;
ll a[100010], b[100010];

ll mul(ll A, ll B, ll mod) //快速乘取余 模板
{
    ll ans = 0;
    while (B > 0) {
        if (B & 1) ans = (ans + A % mod) % mod;
        A = (A + A) % mod;
        B >>= 1;
    }
    return ans;
}

ll exgcd(ll A, ll B, ll &x, ll &y) //扩展欧几里得 模板
{
    if (!B) {
        x = 1, y = 0;
        return A;
    }
    ll d = exgcd(B, A % B, x, y);
    ll tmp = x;
    x = y, y = tmp - A / B * y;
    return d;
}

ll lcm(ll A, ll B) //求最小公倍数
{
    ll xxx, yyy;
    ll g = exgcd(A, B, xxx, yyy);
    return (A / g * B);
}

ll excrt() //重点:求解同余方程组
{
    ll x, y;
    ll M = b[1], ans = a[1]; //赋初值
    //M为前k-1个数的最小公倍数，ans为前k-1个方程的通解
    for (int i = 2; i <= n; i++) {
        ll A = M, B = b[i];
        ll C = (a[i] - ans % B + B) % B; //代表同余方程 ax≡c(mod b) 中a,b,c

        ll g = exgcd(A, B, x, y);
        //求得A,B的最大公约数，与同余方程ax≡gcd(a,b)(mod b)的解，

        if (C % g) return -1; //无解的情况

        x = mul(x, C / g, B); //求得x的值,x即t
        ans += x * M;  //获得前k个方程的通解
        M = lcm(M, B); //更改M的值
        ans = (ans % M + M) % M;
    }
    return ans;
}

int main() {
    scanf("%lld", &n);
    for (int i = 1; i <= n; i++)
        scanf("%lld%lld", &b[i], &a[i]);
    ll ans = excrt();
    printf("%lld", ans);
}
