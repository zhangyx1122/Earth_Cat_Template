ll ch[N][2], f[N], sum[N], val[N], tag[N], siz[N], siz2[N];

inline void pushup(ll p) {
    sum[p] = sum[ch[p][0]] ^ sum[ch[p][1]] ^ val[p];
    siz[p] = siz[ch[p][0]] + siz[ch[p][1]] + 1 + siz2[p];
}

inline void pushdown(ll p) {
    if (tag[p]) {
        if (ch[p][0]) swap(ch[ch[p][0]][0], ch[ch[p][0]][1]), tag[ch[p][0]] ^= 1;
        if (ch[p][1]) swap(ch[ch[p][1]][0], ch[ch[p][1]][1]), tag[ch[p][1]] ^= 1;
        tag[p] = 0;
    }
}

ll getch(ll x) { return ch[f[x]][1] == x; }

bool isroot(ll x) { return ch[f[x]][0] != x && ch[f[x]][1] != x; }

inline void rotate(ll x) {
    ll y = f[x], z = f[y], k = getch(x);
    if (!isroot(y)) ch[z][ch[z][1] == y] = x;
    // 上面这句一定要写在前面，普通的Splay是不用的，因为 isRoot  (后面会讲)
    ch[y][k] = ch[x][!k], f[ch[x][!k]] = y;
    ch[x][!k] = y, f[y] = x, f[x] = z;
    pushup(y), pushup(x);
}

// 从上到下一层一层 pushDown 即可
void update(ll p) {
    if (!isroot(p)) update(f[p]);
    pushdown(p);
}

inline void splay(ll x) {
    update(x);  // 马上就能看到啦。 在
    // Splay之前要把旋转会经过的路径上的点都PushDown
    for (ll fa; fa = f[x], !isroot(x); rotate(x)) {
        if (!isroot(fa)) rotate(getch(fa) == getch(x) ? fa : x);
    }
}

// 回顾一下代码
inline void access(ll x) {
    for (ll p = 0; x; p = x, x = f[x]) {
        splay(x), siz2[x] += siz[ch[x][1]] - siz[p], ch[x][1] = p, pushup(x);
    }
}

inline void makeroot(ll p) {
    access(p);
    splay(p);
    swap(ch[p][0], ch[p][1]);
    tag[p] ^= 1;
}

inline void split(ll a, ll b) {
    makeroot(a);
    access(b);
    splay(b);
}


inline ll find(ll p) {
    access(p), splay(p);
    while (ch[p][0]) pushdown(p), p = ch[p][0];
    splay(p);
    return p;
}

inline void link(ll x, ll y) {
    makeroot(y);
    makeroot(x);
    if (find(y) != x) {
        f[x] = y;
        siz2[y] += siz[x];
    }
}

inline void cut(ll x, ll y) {
    makeroot(x);
    if (find(y) == x && f[y] == x) {
        ch[x][1] = f[y] = 0;
        pushup(x);
    }
}

void init(int n) {
    for (int i = 1; i <= n; i++) siz[i] = 1;
}
