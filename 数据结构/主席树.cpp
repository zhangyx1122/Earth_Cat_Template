#include <bits/stdc++.h>

using namespace std;
typedef long long ll;
const ll N = 1 << 20;

ll ch[N << 5][2], rt[N], tot;
ll val[N << 5];

ll update(ll a, ll b) {
    return a + b;
}

ll build(ll l, ll r) {  // 建树
    ll p = ++tot;
    if (l == r) {
        //初始化
        val[p] = 0;
        return p;
    }
    ll mid = (l + r) >> 1;
    ch[p][0] = build(l, mid);
    ch[p][1] = build(mid + 1, r);
    val[p] = update(val[ch[p][0]], val[ch[p][1]]);
    return p;  // 返回该子树的根节点
}

ll modify(ll pre, ll l, ll r, ll pos, ll v) {  // 插入操作
    ll now = ++tot;
    ch[now][0] = ch[pre][0], ch[now][1] = ch[pre][1];
    if (l == r) {
        val[now] = val[pre] + v;
        return now;
    }
    ll mid = (l + r) >> 1;
    if (pos <= mid)
        ch[now][0] = modify(ch[now][0], l, mid, pos, v);
    else
        ch[now][1] = modify(ch[now][1], mid + 1, r, pos, v);
    val[now] = update(val[ch[now][0]], val[ch[now][1]]);
    return now;
}

ll kth(ll pre, ll now, ll l, ll r, ll k) {  // 查询操作
    ll mid = (l + r) >> 1;
    ll x = val[ch[now][0]] - val[ch[pre][0]];  // 通过区间减法得到左儿子的信息
    if (l == r) return l;
    if (k <= x)  // 说明在左儿子中
        return kth(ch[pre][0], ch[now][0], l, mid, k);
    else  // 说明在右儿子中
        return kth(ch[pre][1], ch[now][1], mid + 1, r, k - x);
}

ll query(ll pre, ll now, ll l, ll r, ll ql, ll qr) {  // 查询操作
    if (ql <= l && r <= qr) {
        return val[now] - val[pre];
    }
    if (qr < l || r < ql) {
        return 0;
    }
    ll mid = (l + r) >> 1;
    ll lv = query(ch[pre][0], ch[now][0], l, mid, ql, qr);
    ll rv = query(ch[pre][1], ch[now][1], mid + 1, r, ql, qr);
    return update(lv, rv);
}
//修改查询记得用rt[]!!!
