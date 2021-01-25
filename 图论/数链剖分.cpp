//洛谷p3384
#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
const ll N=1<<20;
const ll inf=0x3f3f3f3f3f3f3f3f;

int fa[N],son[N],dep[N],siz[N],dfn[N],rnk[N],top[N];
int dfscnt;
int n,m,r;
ll p;
vector<int> g[N];

ll tree[N<<2],lazy[N<<2];

void pushdown(int k,int l,int r){
    if(k>=N) return;
    lazy[k<<1]=(lazy[k<<1]+lazy[k])%p;
    lazy[k<<1|1]=(lazy[k<<1|1]+lazy[k])%p;
    int mid=(l+r)>>1;
    tree[k<<1]=(tree[k<<1]+lazy[k]*(mid-l+1))%p;
    tree[k<<1|1]=(tree[k<<1|1]+lazy[k]*(r-mid))%p;
    lazy[k]=0;
}

void pointadd(int k,int l,int r,int pos,int x){
    if(l==pos&&r==pos){
        tree[k]=(tree[k]+x)%p;
        lazy[k]=(lazy[k]+x)%p;
        return;
    }
    if(r<pos||pos<l) return;
    pushdown(k, l, r);
    int mid=(l+r)>>1;
    pointadd(k<<1,l,mid,pos,x);
    pointadd(k<<1|1,mid+1,r,pos,x);
    tree[k]=(tree[k<<1]+tree[k<<1|1])%p;
}

void invadd(int k,int l,int r,int ql,int qr,int x){
    if(ql<=l&&r<=qr){
        tree[k]=(tree[k]+(r-l+1)*x)%p;
        lazy[k]=(lazy[k]+x)%p;
        return;
    }
    if(r<ql||qr<l) return;
    pushdown(k, l, r);
    int mid=(l+r)>>1;
    invadd(k<<1,l,mid,ql,qr,x);
    invadd(k<<1|1,mid+1,r,ql,qr,x);
    tree[k]=(tree[k<<1]+tree[k<<1|1])%p;
}

ll getsum(int k,int l,int r,int ql,int qr){
    if(ql<=l&&r<=qr) return tree[k];
    if(r<ql||qr<l) return 0;
    pushdown(k, l, r);
    int mid=(l+r)>>1;
    return (getsum(k<<1,l,mid,ql,qr)+getsum(k<<1|1,mid+1,r,ql,qr))%p;
}

void init(int n){
    for(int i=0;i<=n;i++) g[i].clear();
    for(int i=0;i<(N<<2);i++) tree[i]=lazy[i]=0;
    dfscnt=0;
}

void dfs1(int u,int f,int d){
    son[u]=-1;
    siz[u]=1;
    fa[u]=f;
    dep[u]=d;
    for(auto v:g[u]){
        if(v==f) continue;
        if(son[u]==-1||siz[v]>siz[son[u]]) son[u]=v;
        dfs1(v, u,d+1);
        siz[u]+=siz[v];
    }
}

void dfs2(int u,int t){
    dfn[u]=++dfscnt;
    rnk[dfscnt]=u;
    top[u]=t;
    if(son[u]==-1) return;
    dfs2(son[u],t);
    for(auto v:g[u]){
        if(v==son[u]||v==fa[u]) continue;
        dfs2(v,v);
    }
}

int lca(int a,int b){
    while(top[a]!=top[b]){
        if(dep[top[a]]<dep[top[b]]) swap(a,b);
        a=fa[top[a]];
    }
    return dep[a]<dep[b]?a:b;
}

void treeadd(int u,int x){
    invadd(1,1,N,dfn[u],dfn[u]+siz[u]-1,x);
}

ll treesum(int u){
    return getsum(1,1,N,dfn[u],dfn[u]+siz[u]-1);
}

void pathadd(int a,int b,int x){
    while(top[a]!=top[b]){
        if(dep[top[a]]<dep[top[b]]) swap(a,b);
        invadd(1,1,N,dfn[top[a]],dfn[a],x);
        a=fa[top[a]];
    }
    if(dep[a]>dep[b]) swap(a,b);
    invadd(1,1,N,dfn[a],dfn[b],x);
}

ll pathsum(int a,int b){
    ll ans=0;
    while(top[a]!=top[b]){
        if(dep[top[a]]<dep[top[b]]) swap(a,b);
        ans=(ans+getsum(1,1,N,dfn[top[a]],dfn[a]))%p;
        a=fa[top[a]];
    }
    if(dep[a]>dep[b]) swap(a,b);
    ans=(ans+getsum(1,1,N,dfn[a],dfn[b]))%p;
    return ans;
}
int chu[N];
int main(){
    ios::sync_with_stdio(0);
    cin>>n>>m>>r>>p;
    for(int i=1;i<=n;i++) cin>>chu[i];
    init(n);
    for(int i=1;i<n;i++){
        int a,b;
        cin>>a>>b;
        g[a].push_back(b);
        g[b].push_back(a);
    }
    dfs1(r,-1,1);
    dfs2(r,r);
    for(int i=1;i<=n;i++){
        pointadd(1,1,N,dfn[i],chu[i]);
    }
    for(int i=1;i<=m;i++){
        int cho;
        cin>>cho;
        if(cho==1){
            int a,b,x;
            cin>>a>>b>>x;
            pathadd(a, b, x);
        }else if(cho==2){
            int a,b;
            cin>>a>>b;
            cout<<pathsum(a,b)<<endl;
        }else if(cho==3){
            int u,x;
            cin>>u>>x;
            treeadd(u,x);
        }else {
            int u;
            cin>>u;
            cout<<treesum(u)<<endl;
        }
    }
}
