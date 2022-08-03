///2022杭电多校5-1002
#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
typedef pair<int,int> pii;

const int inf=0x3f3f3f3f,N=1e7+9;
const ll mod=4179340454199820289;

const int PMAX=N,PN=N;//PN开到n以内P的最大数量可以省空间
int prime[PN],pcnt;//prime[0]=1,prime[1]=2
bool notp[PMAX];//motp[1]=1,notp[2]=0,notp[4]=1
void Prime(){
    pcnt=0;
    for(int i=2;i<PMAX;i++){
		if(!notp[i])prime[++pcnt]=i;
		for(int j=1;j<=pcnt&&i*prime[j]<PMAX;j++){
			notp[i*prime[j]]=1;
			if(i%prime[j]==0)break;
		}
	}
	notp[1]=1;
}

ll qpow(ll a,ll b){
    ll ans=1;
    while(b){
        if(b&1)ans=(__int128)ans*a%mod;
        a=(__int128)a*a%mod;
        b>>=1;
    }
    return ans;
}

struct PNS{//修改G，f，g即可使用，不同的n只需要初始化ans和n。
    ll ans,n;
    ll G(ll x){
        return (__int128)x*(x+1)/2%mod;
    }
    //f和g均只在h函数中使用，如果可以直接公式求h函数则不用定义这两个函数
    ll f(int pid,ll c){
        return (__int128)qpow(prime[pid],c)*qpow(c,mod-2)%mod;
    }
    //g函数不用记忆化，实际调用次数很少
    ll g(int pid,ll c){
        return qpow(prime[pid],c);
    }

    vector<ll>vh[PN];//这里要确保c>2时调用h(pid,c)前调用过h(pid,c-1)
    ll h(int pid,ll c){
        if(c==0)return 1;
        if(c==1)return 0;
        if(c-2>=(ll)vh[pid].size()){//n=1e12、1e13、1e14时会进80070、230567、670121次，跑不满根号n，所以递推的耗时也是不高的。
            //vh[pid].push_back((ll)((-(__int128)qpow(prime[pid],c)*qpow(c*(c-1),mod-2)%mod+mod)%mod));
            //递推h函数，需要配合f函数和g函数一起使用
            ll ans=f(pid,c);
            for(ll i=1;i<=c;i++){
                ans=((ans-(__int128)g(pid,i)*h(pid,c-i))%mod+mod)%mod;
            }
            vh[pid].push_back(ans);
        }
        return vh[pid][c-2];
    }
    void dfs(ll prod,ll hprod,int pid){
        ans=(ans+(__int128)hprod*G(n/prod))%mod;
        for(int i=pid;i<=pcnt;i++){
            if(prod>n/prime[i]/prime[i])break;
            for(ll c=2,x=prod*prime[i]*prime[i];x<=n;c++,x*=prime[i]){
                dfs(x,(__int128)hprod*h(i,c)%mod,i+1);
                if(x>n/prime[i])break;
            }
        }
    }
}pns;

int main(){
    #ifdef ONLINE_JUDGE
        //std::ios::sync_with_stdio(false);
    #else
        freopen("1002.in","r",stdin);
        //freopen("out.txt","w",stdout);
    #endif
    Prime();
    int t;
    scanf("%d",&t);
    while(t--){
        pns.ans=0;
        scanf("%lld",&pns.n);
        pns.dfs(1,1,1);
        //printf("%lld\n",pns.ans);
        printf("%lld\n",(ll)((__int128)pns.ans*qpow(pns.n,mod-2)%mod));
    }
    return 0;
}
