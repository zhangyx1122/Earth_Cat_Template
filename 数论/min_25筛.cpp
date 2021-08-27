/*
https://loj.ac/p/6053
筛积性函数f的前缀和
f(1)=1
f(p^e)=f xor e
n<=1e10，LOJ 347ms本地1100ms
*/
#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
const ll mod=1e9+7,inv3=333333336;
const int N=1e5+5;//开到sqrt(n)即可

ll prime[N],sp0[N],sp1[N],sp2[N],g0[N<<1],g1[N<<1],g2[N<<1];
ll pnum,min25n,sqrn,w[N<<1],ind1[N],ind2[N];
bool notp[N];

void pre() { //预处理，线性筛
    notp[1]=1;
    for(int i=1; i<N; i++) {
        if(!notp[i]) {
            prime[++pnum]=i;
            sp0[pnum]=(sp0[pnum-1]+1)%mod;//p^0前缀和（p指质数），可以按需增删，下标意义为第pnum个质数的前缀和，而g的实际下标意义为w之前的前缀和，两者有所区别
            sp1[pnum]=(sp1[pnum-1]+i)%mod;//p^1前缀和
            sp2[pnum]=(sp2[pnum-1]+1ll*i*i)%mod;//p^2前缀和
        }
        for(int j=1; j<=pnum&&prime[j]*i<N; j++) {
            notp[i*prime[j]]=1;
            if(i%prime[j]==0)break;
        }
    }
}

void min25(ll n) {
    ll tot=0;
    min25n=n;
    sqrn=sqrt(n);
    for(ll i=1; i<=n; i=n/(n/i)+1) {
        w[++tot]=n/i;//实际下标
        ll x=w[tot]%mod;
        g0[tot]=x-1;//x^0前缀和
        g1[tot]=x*(x+1)/2%mod-1;//x^1前缀和
        g2[tot]=x*(x+1)/2%mod*(2*x+1)%mod*inv3%mod-1;//x^2前缀和
        if(n/i<=sqrn)ind1[n/i]=tot;//离散下标
        else ind2[n/(n/i)]=tot;//离散下标
    }
    for(int i=1; i<=pnum; i++) {//扩展埃氏筛，筛质数部分前缀和
        for(int j=1; j<=tot&&prime[i]*prime[i]<=w[j]; j++) {
            int id=w[j]/prime[i]<=sqrn?ind1[w[j]/prime[i]]:ind2[n/(w[j]/prime[i])];
            g0[j]-=(g0[id]-sp0[i-1]+mod)%mod;
            g1[j]-=prime[i]*(g1[id]-sp1[i-1]+mod)%mod;
            g2[j]-=prime[i]*prime[i]%mod*(g2[id]-sp2[i-1]+mod)%mod;
            g0[j]%=mod,g1[j]%=mod,g2[j]%=mod;
            if(g0[j]<0)g0[j]+=mod;
            if(g1[j]<0)g1[j]+=mod;
            if(g2[j]<0)g2[j]+=mod;
        }
    }
}

//该前缀和不计算f(1)，需要自行加上
ll S(ll x,int y) {//x以内最小质因子大于第y个因子的前缀和
    if(prime[y]>=x)return 0;
    int id=x<=sqrn?ind1[x]:ind2[min25n/x];
    ll ans=(((g1[id]-g0[id])-(sp1[y]-sp0[y]))%mod+mod)%mod;//x以内大于第y个因子的质数部分前缀和
    if(x>=2&&y<1)ans=(ans+2)%mod;//特判包含f(2)的情况
    for(int i=y+1; i<=pnum&&prime[i]*prime[i]<=x; i++) {//筛合数部分前缀和
        ll pe=prime[i];
        for(int e=1; pe<=x; e++,pe=pe*prime[i]) {
            ll fpe=prime[i]^e;//f(p^e)
            ans=(ans+fpe%mod*(S(x/pe,i)+(e!=1)))%mod;
        }
    }
    return ans%mod;
}

int main() {
    pre();//预处理一次即可
    ll n;
    scanf("%lld",&n);
    min25(n);//每个不同的n都要调用一次该函数，再调用S(n,0)
    printf("%lld\n",S(n,0)+1);//加上f(1)
    return 0;
}
