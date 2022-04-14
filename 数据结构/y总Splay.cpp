#include <bits/stdc++.h>

using namespace std;
const int N = 1e6 + 10;
struct node {
    int p, v, s[2];
    int siz, tag;
    void init(int _v, int _p) {
        v = _v, p = _p;
        siz = 1;
    }
};
node tr[N];
int root, idx;

void pushup(int x) { tr[x].siz = tr[tr[x].s[0]].siz + tr[tr[x].s[1]].siz + 1; }

void pushdown(int x) {
    if (tr[x].tag) {
        swap(tr[x].s[0], tr[x].s[1]);
        tr[tr[x].s[0]].tag ^= 1;
        tr[tr[x].s[1]].tag ^= 1;
        tr[x].tag = 0;
    }
}
void rotate(int x) {
	pushdown(x); 
    int y = tr[x].p, z = tr[y].p;
    int k = tr[y].s[1] == x;
    tr[y].s[k] = tr[x].s[k ^ 1], tr[tr[y].s[k]].p = y;
    tr[x].s[k ^ 1] = y, tr[y].p = x;
    tr[z].s[tr[z].s[1] == y] = x, tr[x].p = z;
    pushup(y), pushup(x);
}

void splay(int x, int k) {
    while (tr[x].p != k) {
        int y = tr[x].p, z = tr[y].p;
        if (z != k) {
            if ((tr[z].s[1] == y) ^ (tr[y].s[1] == x)) {
                rotate(x);
            } else {
                rotate(y);
            }
        }
        rotate(x);
    }
    if (!k) root = x;
}

void insert(int v) {
    int u = root, p = 0;
    while (u) p = u, u = tr[u].s[v > tr[u].v];
    u = ++idx;
    if (p) tr[p].s[v > tr[p].v] = u;
    tr[u].init(v, p);
    splay(u, 0);
}

int getk(int k) {
    int u = root;
    while (1) {
        pushdown(u);
        if (k <= tr[tr[u].s[0]].siz) {
            u = tr[u].s[0];
        } else if (k == tr[tr[u].s[0]].siz + 1) {
            splay(u, 0);
            return u;
        } else {
            k -= tr[tr[u].s[0]].siz + 1, u = tr[u].s[1];
        }
    }
}

int n, m;
void output(int u) {
    if (u == 0) return;
    pushdown(u);
    output(tr[u].s[0]);
    if (1 <= tr[u].v && tr[u].v <= n) cout << tr[u].v << ' ';
    output(tr[u].s[1]);
}

int main() {
    ios::sync_with_stdio(0), cin.tie(0), cout.tie(0);
    cin >> n >> m;
    for (int i = 0; i <= n + 1; i++) insert(i);
    while (m--) {
        int a, b;
        cin >> a >> b;
        int id1 = getk(a), id2 = getk(b + 2);
        splay(id1, 0), splay(id2, id1);
        tr[tr[id2].s[0]].tag ^= 1;
    }
    output(root);
}
