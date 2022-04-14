#include <bits/stdc++.h>

using namespace std;
typedef long long ll;

const int N = 4e6 + 10;
const ll mod = 998244353, G = 3, Gi = 332748118;

int limit = 1, L, r[N];
ll a[N], b[N];

ll qpow(ll _a, ll _b) {
    ll ans = 1;
    while (_b) {
        if (_b & 1) ans = (ans * _a) % mod;
        _b >>= 1;
        _a = (_a * _a) % mod;
    }
    return ans;
}

void ntt(ll *A, int type) {
    auto swap = [](ll &_a, ll &_b) {
        _a ^= _b, _b ^= _a, _a ^= _b;
    };
    for (int i = 0; i < limit; i++)
        if (i < r[i]) swap(A[i], A[r[i]]);
    for (int mid = 1; mid < limit; mid <<= 1) {
        ll Wn = qpow(type == 1 ? G : Gi, (mod - 1) / (mid << 1));
        for (int j = 0; j < limit; j += (mid << 1)) {
            ll w = 1;
            for (int k = 0; k < mid; k++, w = (w * Wn) % mod) {
                int x = A[j + k], y = w * A[j + k + mid] % mod;
                A[j + k] = (x + y) % mod,
                        A[j + k + mid] = (x - y + mod) % mod;
            }
        }
    }
}

void NTT(int n, int m) {
    limit = 1;
    L = 0;
    while (limit <= n + m) limit <<= 1, L++;
    for (int i = 0; i < limit; i++) r[i] = (r[i >> 1] >> 1) | ((i & 1) << (L - 1));
    ntt(a, 1), ntt(b, 1);
    for (int i = 0; i < limit; i++) a[i] = (a[i] * b[i]) % mod;
    ntt(a, -1);
    ll inv = qpow(limit, mod - 2);
    for (int i = 0; i <= n + m; i++) a[i] = a[i] * inv % mod;
}

int main() {
    int n, m;
    cin >> n >> m;
    for (int i = 0; i <= n; i++) {
        cin >> a[i];
        a[i] = (a[i] + mod) % mod;
    }
    for (int i = 0; i <= m; i++) {
        cin >> b[i];
        b[i] = (b[i] + mod) % mod;
    }
    NTT(n, m);
    for (int i = 0; i <= n + m; i++) cout << a[i] << ' ';
}

/*
#include<cstdio>
#include<cctype>
#include<cstring>
#include<cmath>
namespace fast_IO
{
	const int IN_LEN=10000000,OUT_LEN=10000000;
	char ibuf[IN_LEN],obuf[OUT_LEN],*ih=ibuf+IN_LEN,*oh=obuf,*lastin=ibuf+IN_LEN,*lastout=obuf+OUT_LEN-1;
	inline char getchar_(){return (ih==lastin)&&(lastin=(ih=ibuf)+fread(ibuf,1,IN_LEN,stdin),ih==lastin)?EOF:*ih++;}
	inline void putchar_(const char x){if(oh==lastout)fwrite(obuf,1,oh-obuf,stdout),oh=obuf;*oh++=x;}
	inline void flush(){fwrite(obuf,1,oh-obuf,stdout);}
}
using namespace fast_IO;
#define getchar() getchar_()
#define putchar(x) putchar_((x))
typedef long long LL;
#define rg register
template <typename T> inline T max(const T a,const T b){return a>b?a:b;}
template <typename T> inline T min(const T a,const T b){return a<b?a:b;}
template <typename T> inline T mind(T&a,const T b){a=a<b?a:b;}
template <typename T> inline T maxd(T&a,const T b){a=a>b?a:b;}
template <typename T> inline T abs(const T a){return a>0?a:-a;}
template <typename T> inline void swap(T&a,T&b){T c=a;a=b;b=c;}
template <typename T> inline void swap(T*a,T*b){T c=a;a=b;b=c;}
template <typename T> inline T gcd(const T a,const T b){if(!b)return a;return gcd(b,a%b);}
template <typename T> inline T square(const T x){return x*x;};
template <typename T> inline void read(T&x)
{
    char cu=getchar();x=0;bool fla=0;
    while(!isdigit(cu)){if(cu=='-')fla=1;cu=getchar();}
    while(isdigit(cu))x=x*10+cu-'0',cu=getchar();
    if(fla)x=-x;  
}
template <typename T> void printe(const T x)
{
    if(x>=10)printe(x/10);
    putchar(x%10+'0');
}
template <typename T> inline void print(const T x)
{
    if(x<0)putchar('-'),printe(-x);
    else printe(x);
}
const int maxn=262145;
int n,m;
struct Ntt
{
	LL mod,a[maxn],b[maxn];;
	inline LL pow(LL x,LL y)
	{
		rg LL res=1;
		for(;y;y>>=1,x=x*x%mod)if(y&1)res=res*x%mod;
		return res;
	}
	int lenth,Reverse[maxn];
	inline void init(const int x)
	{
		rg int tim=0;lenth=1;
		while(lenth<=x)lenth<<=1,tim++;
		for(rg int i=0;i<lenth;i++)Reverse[i]=(Reverse[i>>1]>>1)|((i&1)<<(tim-1));
	}
	inline void NTT(LL*A,const int fla)
	{
		for(rg int i=0;i<lenth;i++)if(i<Reverse[i])swap(A[i],A[Reverse[i]]);
		for(rg int i=1;i<lenth;i<<=1)
		{
			LL w=pow(3,(mod-1)/i/2);
			if(fla==-1)w=pow(w,mod-2);
			for(rg int j=0;j<lenth;j+=(i<<1))
			{
				LL K=1;
				for(rg int k=0;k<i;k++,K=K*w%mod)
				{
					const LL x=A[j+k],y=A[j+k+i]*K%mod;
					A[j+k]=(x+y)%mod;
					A[j+k+i]=(mod+x-y)%mod;
				}
			}
		}
		if(fla==-1)
		{
			const int inv=pow(lenth,mod-2);
			for(rg int i=0;i<lenth;i++)A[i]=A[i]*inv%mod;
		}	
	}
}Q[3];
LL EXgcd(const LL a,const LL b,LL &x,LL &y)  
{  
    if(!b)
    {
		x=1,y=0;
        return a;  
    }
    const LL res=EXgcd(b,a%b,y,x);
    y-=a/b*x;
    return res;
}
inline LL msc(LL a,LL b,LL mod)
{
    LL v=(a*b-(LL)((long double)a/mod*b+1e-8)*mod);
    return v<0?v+mod:v;
}
int N,a[3],p[3];
LL CRT()
{  
    LL P=1,sum=0;  
    for(rg int i=1;i<=N;i++)P*=p[i];
    for(rg int i=1;i<=N;i++)  
	{
    	const LL m=P/p[i];
    	LL x,y;
    	EXgcd(p[i],m,x,y);
    	sum=(sum+msc(msc(y,m,P),a[i],P))%P;
	}
	return sum;
}
int P;
int main()
{
	read(n),read(m),read(P);
	Q[0].mod=469762049,Q[0].init(n+m);
	Q[1].mod=998244353,Q[1].init(n+m);
	Q[2].mod=1004535809,Q[2].init(n+m);
	for(rg int i=0;i<=n;i++)read(Q[0].a[i]),Q[2].a[i]=Q[1].a[i]=Q[0].a[i];
	for(rg int i=0;i<=m;i++)read(Q[0].b[i]),Q[2].b[i]=Q[1].b[i]=Q[0].b[i];
	Q[0].NTT(Q[0].a,1),Q[0].NTT(Q[0].b,1);
	Q[1].NTT(Q[1].a,1),Q[1].NTT(Q[1].b,1);
	Q[2].NTT(Q[2].a,1),Q[2].NTT(Q[2].b,1);
	for(rg int i=0;i<Q[0].lenth;i++)
		Q[0].a[i]=(LL)Q[0].a[i]*Q[0].b[i]%Q[0].mod,
		Q[1].a[i]=(LL)Q[1].a[i]*Q[1].b[i]%Q[1].mod,
		Q[2].a[i]=(LL)Q[2].a[i]*Q[2].b[i]%Q[2].mod;
	Q[0].NTT(Q[0].a,-1);
	Q[1].NTT(Q[1].a,-1);
	Q[2].NTT(Q[2].a,-1);
	N=2,p[1]=Q[0].mod,p[2]=Q[1].mod;
	const int INV=Q[2].pow(Q[0].mod,Q[2].mod-2)*Q[2].pow(Q[1].mod,Q[2].mod-2)%Q[2].mod;
	for(rg int i=0;i<=n+m;i++)
	{
		a[1]=Q[0].a[i],a[2]=Q[1].a[i];
		const LL ans1=CRT();
		const LL ans2=((Q[2].a[i]-ans1)%Q[2].mod+Q[2].mod)%Q[2].mod*INV%Q[2].mod;
		print((ans2*Q[0].mod%P*Q[1].mod%P+ans1)%P),putchar(' ');
	}
	return flush(),0;
}

*/
