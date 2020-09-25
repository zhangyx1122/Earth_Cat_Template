#include <bits/stdc++.h>
using namespace std;

struct Link_Cut_Tree {
#define ls ch[p][0]
#define rs ch[p][1]
    static const int N = 200005;
    int ch[N][2], f[N], sum[N], val[N], tag[N], laz[N], siz[N];

    inline void pushup(int p) {
        // maintain other variables
        siz[p] = siz[ls] + siz[rs];
    }

    inline void pushdown(int p) {}

    int get(int x) { return ch[f[x]][1] == x; }

    bool isroot(int x) { return ch[f[x]][0] != x && ch[f[x]][1] != x; }

    inline void rotate(int x) {
        int y = f[x], z = f[y], k = get(x);
        if (!isroot(y)) ch[z][ch[z][1] == y] = x;
        // 上面这句一定要写在前面，普通的Splay是不用的，因为 isRoot  (后面会讲)
        ch[y][k] = ch[x][!k], f[ch[x][!k]] = y;
        ch[x][!k] = y, f[y] = x, f[x] = z;
        pushup(x), pushup(y);
    }

    // 从上到下一层一层 pushDown 即可
    void update(int p) {
        if (!isroot(p)) update(f[p]);
        pushdown(p);
    }

    inline void splay(int x) {
        update(x);  // 马上就能看到啦。 在
                    // Splay之前要把旋转会经过的路径上的点都PushDown
        for (int fa; fa = f[x], !isroot(x); rotate(x)) {
            if (!isroot(fa)) rotate(get(fa) == get(x) ? fa : x);
        }
    }

    // 回顾一下代码
    inline int access(int x) {
        int p;
        for (p = 0; x; p = x, x = f[x]) {
            splay(x), ch[x][1] = p, pushup(x);
        }
        return p;
    }

    inline void makeroot(int p) {
        p = access(p);
        swap(ch[p][0], ch[p][1]);
        tag[p] ^= 1;
    }

    inline void link(int x, int p) {
        makeroot(x);
        splay(x);
        f[x] = p;
    }

    inline void cut(int x, int p) {
        makeroot(x), access(p), splay(p), ls = f[x] = 0;
    }

    inline int find(int p) {
        access(p), splay(p);
        while (ls) pushdown(p), p = ls;
        splay(p);
        return p;
    }
};
