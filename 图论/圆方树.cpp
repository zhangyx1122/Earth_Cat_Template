///
#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
typedef pair<int,int> pii;

const int inf=0x3f3f3f3f,N=2e5+9;

int n,m;
vector<int>e[N],g[N*2];

int cnt,dfn[N],low[N],dfc;
int stk[N],tp;

void tarjan(int u) {
    low[u] = dfn[u] = ++dfc;                // low 初始化为当前节点 dfn
    stk[++tp] = u;                          // 加入栈中
    for (int v : e[u]) {                    // 遍历 u 的相邻节点
        if (!dfn[v]) {                        // 如果未访问过
            tarjan(v);                          // 递归
            low[u] = std::min(low[u], low[v]);  // 未访问的和 low 取 min
            if (low[v] == dfn[u]) {  // 标志着找到一个以 u 为根的点双连通分量
                ++cnt;                 // 增加方点个数
                // 将点双中除了 u 的点退栈，并在圆方树中连边
                for (int x = 0; x != v; --tp) {
                    x = stk[tp];
                    g[cnt].push_back(x);
                    g[x].push_back(cnt);
                }
                // 注意 u 自身也要连边（但不退栈）
                g[cnt].push_back(u);
                g[u].push_back(cnt);
            }
        }
        else low[u] = std::min(low[u], dfn[v]); // 已访问的和 dfn 取 min
    }
}

int main() {
    #ifdef ONLINE_JUDGE
        //std::ios::sync_with_stdio(false);
    #else
        //freopen("in.txt","r",stdin);
        //freopen("out.txt","w",stdout);
    #endif
    scanf("%d%d", &n, &m);
    for (int i = 0; i < m; ++i) {
        int u, v;
        scanf("%d%d", &u, &v);
        e[u].push_back(v);  // 加双向边
        e[v].push_back(u);
    }
    cnt = n;  // 点双 / 方点标号从 N 开始
    for (int u = 1; u <= n; ++u)// 处理非连通图
        if (!dfn[u]){
            tarjan(u),--tp;
        }
    return 0;
}
