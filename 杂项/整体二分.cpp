ll bit[N];

void add_bit(ll k, ll a) {
    while (k < N) {
        bit[k] = bit[k] + a;
        k += k & -k;
    }
}

ll query_bit(ll k) {
    ll ans = 0;
    while (k) {
        ans = ans + bit[k];
        k -= k & -k;
    }
    return ans;
}

struct node {
    ll x, y, k, id, type;
};
node q[N], q1[N], q2[N];
ll ans[N], now[N], tot, totx;

void solve(ll l, ll r, ll ql, ll qr) {
    if (ql > qr) return;
    if (l == r) {
        for (ll i = ql; i <= qr; i++) {
            if (q[i].type == 2) {
                ans[q[i].id] = l;
            }
        }
        return;
    }
    ll mid = (l + r) >> 1;
    ll cq1 = 0, cq2 = 0;
    for (ll i = ql; i <= qr; i++) {
        if (q[i].type == 1) {
            if (q[i].y <= mid) {
                add_bit(q[i].x, q[i].k);
                q1[++cq1] = q[i];
            } else {
                q2[++cq2] = q[i];
            }
        } else {
            ll sum = query_bit(q[i].y) - query_bit(q[i].x - 1);
            if (sum >= q[i].k) {
                q1[++cq1] = q[i];
            } else {
                q2[++cq2] = q[i];
                q2[cq2].k -= sum;
            }
        }
    }
    for (ll i = 1; i <= cq1; i++) if (q1[i].type == 1) add_bit(q1[i].x, -q1[i].k);
    for (ll i = 1; i <= cq1; i++) q[ql + i - 1] = q1[i];
    for (ll i = 1; i <= cq2; i++) q[ql + cq1 + i - 1] = q2[i];
    solve(l, mid, ql, ql + cq1 - 1);
    solve(mid + 1, r, ql + cq1, qr);

}

void init() {
    totx = 0;
    tot = 0;
    memset(bit, 0, sizeof bit);
}
