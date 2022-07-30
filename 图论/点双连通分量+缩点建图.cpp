#include <bits/stdc++.h>

#define ll long long
using namespace std;
const int N = 10010;
const int M = 10010 * 4;
int head[N];
int ver[M];
int Next[M];
int tot, n, m;

void add(int x, int y) {
    ver[++tot] = y;
    Next[tot] = head[x];
    head[x] = tot;
}

int root;
vector<int> dcc[N];
int stackk[N];
int dfn[N], low[N];
int num = 0;//时间戳
int top;//stackk
int cnt = 0;//联通块数目
bool cut[N];//割点判断
void tarjan(int x) {
    dfn[x] = low[x] = ++num;
    stackk[++top] = x;
    if (x == root && head[x] == 0) {
        dcc[++cnt].push_back(x);//cnt联通块标号
        return;
    }
    int flag = 0;
    for (int i = head[x]; i; i = Next[i]) {
        int y = ver[i];
        if (!dfn[y]) {
            tarjan(y);
            low[x] = min(low[x], low[y]);
            if (low[y] >= dfn[x]) {
                flag++;
                if (x != root || flag > 1)cut[x] = true;
                cnt++;
                int z;
                do//弹出的元素与x一起构成一个联通块(或者说割点的子树中的节点+割点?)
                {
                    z = stackk[top--];
                    dcc[cnt].push_back(z);
                } while (z != y);
                dcc[cnt].push_back(x);
            }
        } else low[x] = min(low[x], dfn[y]);
    }
}

int tot2 = 1;
int new_id[N];

int hc[N];
int vc[M];
int nc[M];

void add_c(int x, int y) {
    vc[++tot2] = y;
    nc[tot2] = hc[x];
    hc[x] = tot2;
}

int main() {
    while (cin >> n >> m) {
        tot = 1;//方便用^运算访问各边的终点
        for (int i = 1; i <= m; ++i) {
            int x, y;
            cin >> x >> y;
            if (x == y)continue;
            add(x, y);
            add(y, x);
        }
        for (int i = 1; i <= n; ++i) {
            if (!dfn[i])root = i, tarjan(i);
        }
        /*for(int i=1;i<=n;++i)
            if(cut[i])printf("%d ",i);*/
        //上面求割点同时求V-DCC
        //下面输出每个联通块中的点
        for (int i = 1; i <= cnt; ++i) {
            for (int j = 0; j < dcc[i].size(); ++j)cout << i << " " << dcc[i][j] << endl;
        }

        //缩点
        tot2 = 1;
        int num2 = cnt;
        for (int i = 1; i <= n; ++i) {
            if (cut[i])new_id[i] = ++num2;//缩点后割点的新编号,相当于每个割点单独作为一个联通块
        }
        for (int i = 1; i <= cnt; ++i) {
            for (int j = 0; j < dcc[i].size(); ++j) {
                int x = dcc[i][j];
                if (cut[x])//一个联通块中有且只有一个割点，通过割点们把这些联通块连接起来;
                {
                    add_c(i, new_id[x]);
                    add_c(new_id[x], i);
                } else new_id[x] = i;//其余点均只属于一个联通块
            }
        }

        //输出缩点后的图中各点之间的邻接关系，再次注意^符号的使用，i从2开始，每次加2，<tot2而非<=；
        for (int i = 2; i < tot2; i += 2)
            cout << vc[i ^ 1] << "   " << vc[i] << endl;


    }
    return 0;
}

/*
 * tot2为边数，从2开始（？
 * num2为缩点之后的点数
 * 点双缩点，点可能越缩越多，注意N大小
*/
