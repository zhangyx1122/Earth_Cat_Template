#include <bits/stdc++.h>

using namespace std;
typedef long long ll;
const int N = 1 << 20;

struct node {
    int mmax, semax, cnt;
    ll sum;
};

node tree[N << 1];
int init[N << 1];

node merge_range(node a, node b) {
    node ans;
    ans.sum = a.sum + b.sum;
    if (a.mmax == b.mmax) {
        ans.mmax = a.mmax;
        ans.cnt = a.cnt + b.cnt;
        ans.semax = max(a.semax, b.semax);
    } else {
        if (a.mmax < b.mmax) swap(a, b);
        ans.mmax = a.mmax;
        ans.cnt = a.cnt;
        ans.semax = max(a.semax, b.mmax);
    }
    return ans;
}

void build(int k, int l, int r) {
    if (l == r) {
        tree[k] = {init[l], -1, 1, init[l]};
        return;
    }
    int mid = (l + r) >> 1;
    build(k << 1, l, mid);
    build(k << 1 | 1, mid + 1, r);
    tree[k] = merge_range(tree[k << 1], tree[k << 1 | 1]);
}


void pushdown(int k, int l, int r) {
    if (l == r) return;
    if (tree[k].mmax < tree[k << 1].mmax) {
        tree[k << 1].sum -= 1LL * (tree[k << 1].mmax - tree[k].mmax) * tree[k << 1].cnt;
        tree[k << 1].mmax = tree[k].mmax;
    }
    if (tree[k].mmax < tree[k << 1 | 1].mmax) {
        tree[k << 1 | 1].sum -= 1LL * (tree[k << 1 | 1].mmax - tree[k].mmax) * tree[k << 1 | 1].cnt;
        tree[k << 1 | 1].mmax = tree[k].mmax;
    }
}


node query(int k, int l, int r, int ql, int qr) {
    if (qr < l || r < ql) return {0, -1, 1, 0};
    if (ql <= l && r <= qr) {
        return tree[k];
    }
    pushdown(k, l, r);
    int mid = (l + r) >> 1;
    node lq = query(k << 1, l, mid, ql, qr);
    node rq = query(k << 1 | 1, mid + 1, r, ql, qr);
    return merge_range(lq, rq);
}

void modify(int k, int l, int r, int ql, int qr, int x) {
    if (qr < l || r < ql) return;
    if (ql <= l && r <= qr && tree[k].semax < x) {
        if (x < tree[k].mmax) {
            tree[k].sum -= 1LL * (tree[k].mmax - x) * tree[k].cnt;
            tree[k].mmax = x;
        }
        return;
    }
    pushdown(k, l, r);
    int mid = (l + r) >> 1;
    modify(k << 1, l, mid, ql, qr, x);
    modify(k << 1 | 1, mid + 1, r, ql, qr, x);
    tree[k] = merge_range(tree[k << 1], tree[k << 1 | 1]);
}


signed main() {
//    freopen("data.txt", "r", stdin);
//    freopen("test1.txt", "w", stdout);
    int t;
    scanf("%d", &t);
    while (t--) {
        int n, q;
        scanf("%d%d", &n, &q);
        for (int i = 1; i <= n; i++) scanf("%d", &init[i]);
        build(1, 1, n);
        while (q--) {
            int x, y, op, val;
            scanf("%d%d%d", &op, &x, &y);
            if (op == 0) {
                scanf("%d", &val);
                modify(1, 1, n, x, y, val);
            } else if (op == 1) {
                node ans = query(1, 1, n, x, y);
                printf("%d\n", ans.mmax);
            } else {
                node ans = query(1, 1, n, x, y);
                printf("%lld\n", ans.sum);
            }
        }
    }
}
