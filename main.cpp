#include <bits/stdc++.h>
using namespace std;
typedef long long ll;

const ll mod=1e9+7;
const int maxn=2e5+10;
ll invn[maxn];
int init(){
    int len=(int)(2e5+5);
    invn[0]=invn[1]=1;
    for(int i=2;i<=len;++i){
       invn[i]=(mod-mod/i)*invn[mod%i]%mod;
    }
    return 0;
}
ll qpow(ll a,ll b){
    ll res=1;
    while(b){
        if(b&1){
            res=(res*a)%mod;
            b--;
        }
        else{
            a=(a*a)%mod;
            b>>=1;
        }
    }
    return res;
}
ll sum[200005],a[200005];
int main(){
#ifdef ONLINE_JUDGE
#else
    freopen("1.in.txt","r",stdin);
    freopen("out.txt","w",stdout);
#endif
    init();
    ll nn,n;
    scanf("%lld",&nn);
    while(nn--){
        scanf("%lld",&n);
        sum[1]=0;
        for(ll i=1;i<=n;i++){
            scanf("%lld",&a[i]);
            sum[1]+=a[i];
        }
        for(ll i=2;i<=(n+1)/2;i++){
            sum[i]=sum[i-1]-a[i-1]-a[n-i+2];
        }
        sum[1]%=mod;
        for(ll i=2;i<=(n+1)/2;i++)sum[i]=(sum[i]+sum[i-1])%mod;
        ll res=0;
        for(ll i=1;i<=n;i++){
            ll ii=min(i,n-i+1);
            res=(res+sum[ii]*invn[i])%mod;
        }
        ll num=n*(n+1)/2;
        num%=mod;
        num=qpow(num,mod-2);
        printf("%lld\n",(res*num)%mod);
    }
    return 0;
}
