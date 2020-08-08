//大量使用了memset，但常数貌似很小？HDU6808跑了998ms（限制5000ms），然而这个代码不是HDU6808的
#include<bits/stdc++.h>
using namespace std;

const int maxn=505;// 最大点数
const int inf=0x3f3f3f3f;// 距离初始值
struct HK_Hungary{//这个板子从1开始，0点不能用,nx为左边点数，ny为右边点数
    int nx,ny;//左右顶点数量
    vector<int>bmap[maxn];
    int cx[maxn];//cx[i]表示左集合i顶点所匹配的右集合的顶点序号
    int cy[maxn]; //cy[i]表示右集合i顶点所匹配的左集合的顶点序号
    int dx[maxn];
    int dy[maxn];
    int dis;
    bool bmask[maxn];
    void init(int a,int b){
        nx=a,ny=b;
        for(int i=0;i<=nx;i++){
            bmap[i].clear();
        }
    }
    void add_edge(int u,int v){
        bmap[u].push_back(v);
    }
    bool searchpath(){//寻找 增广路径
        queue<int>Q;
        dis=inf;
        memset(dx,-1,sizeof(dx));
        memset(dy,-1,sizeof(dy));
        for(int i=1;i<=nx;i++){//cx[i]表示左集合i顶点所匹配的右集合的顶点序号
            if(cx[i]==-1){//将未遍历的节点 入队 并初始化次节点距离为0
                Q.push(i);
                dx[i]=0;
            }
        }//广度搜索增广路径
        while(!Q.empty()){
            int u=Q.front();
            Q.pop();
            if(dx[u]>dis) break;//取右侧节点
            for(int i=0;i<bmap[u].size();i++){
                int v=bmap[u][i];//右侧节点的增广路径的距离
                if(dy[v]==-1){
                    dy[v]=dx[u]+1;//v对应的距离 为u对应距离加1
                    if(cy[v]==-1)dis=dy[v];
                    else{
                        dx[cy[v]]=dy[v]+1;
                        Q.push(cy[v]);
                    }
                }
            }
        }
        return dis!=inf;
    }
    int findpath(int u){//寻找路径 深度搜索
        for(int i=0;i<bmap[u].size();i++){
            int v=bmap[u][i];//如果该点没有被遍历过 并且距离为上一节点+1
            if(!bmask[v]&&dy[v]==dx[u]+1){//对该点染色
                bmask[v]=1;
                if(cy[v]!=-1&&dy[v]==dis)continue;
                if(cy[v]==-1||findpath(cy[v])){
                    cy[v]=u;cx[u]=v;
                    return 1;
                }
            }
        }
        return 0;
    }
    int MaxMatch(){//得到最大匹配的数目
        int res=0;
        memset(cx,-1,sizeof(cx));
        memset(cy,-1,sizeof(cy));
        while(searchpath()){
            memset(bmask,0,sizeof(bmask));
            for(int i=1;i<=nx;i++){
                if(cx[i]==-1){
                    res+=findpath(i);
                }
            }
        }
        return res;
    }
}HK;

int main(){
    int nn,n,m;
    cin>>nn;
    while(nn--){
        scanf("%d%d",&n,&m);
        HK.init(n,m);//左端点和右端点数量
        for(int i=1;i<=n;i++){
            int snum;
            cin>>snum;
            int v;
            for(int j=1;j<=snum;j++){
                cin>>v;
                HK.add_edge(i,v);//连边
            }
        }
        cout<<HK.MaxMatch()<<endl;//求最大匹配
    }
    return 0;
}
