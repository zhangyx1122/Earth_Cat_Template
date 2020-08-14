#include <bits/stdc++.h>
using namespace std;

struct Splay {
    static const int N = 100005;
    int rt, tot, fa[N], ch[N][2], val[N], cnt[N], sz[N];
    // rt=根编号，tot=总节点，fa=父节点编号，ch=左/右儿子编号，val=节点的值，cnt=权值出现次数，sz=子树大小
    void maintain(int x) {  //更新x节点字数大小
        sz[x] = sz[ch[x][0]] + sz[ch[x][1]] + cnt[x];
    }

    bool get(int x) {
        return x == ch[fa[x]][1];
    }  //返回节点是父亲的0/1-左/右儿子

    void clear(int x) {  //销毁节点x
        ch[x][0] = ch[x][1] = fa[x] = val[x] = sz[x] = cnt[x] = 0;
    }

    void rotate(int x) {  //旋转
        int y = fa[x], z = fa[y], chk = get(x);
        ch[y][chk] = ch[x][chk ^ 1];
        fa[ch[x][chk ^ 1]] = y;
        ch[x][chk ^ 1] = y;
        fa[y] = x;
        fa[x] = z;
        if (z) ch[z][y == ch[z][1]] = x;
        maintain(x);
        maintain(y);
    }

    void splay(int x) {  //将x节点移动到根
        for (int f = fa[x]; f = fa[x], f; rotate(x))
            if (fa[f]) rotate(get(x) == get(f) ? f : x);
        rt = x;
    }

    void ins(int k) {  //插入
        if (!rt) {
            val[++tot] = k;
            cnt[tot]++;
            rt = tot;
            maintain(rt);
            return;
        }
        int cnr = rt, f = 0;
        while (1) {
            if (val[cnr] == k) {
                cnt[cnr]++;
                maintain(cnr);
                maintain(f);
                splay(cnr);
                break;
            }
            f = cnr;
            cnr = ch[cnr][val[cnr] < k];
            if (!cnr) {
                val[++tot] = k;
                cnt[tot]++;
                fa[tot] = f;
                ch[f][val[f] < k] = tot;
                maintain(tot);
                maintain(f);
                splay(tot);
                break;
            }
        }
    }

    int rk(int k) {  // k权值的排名
        int res = 0, cnr = rt;
        while (1) {
            if (k < val[cnr]) {
                cnr = ch[cnr][0];
            } else {
                res += sz[ch[cnr][0]];
                if (k == val[cnr]) {
                    splay(cnr);
                    return res + 1;
                }
                res += cnt[cnr];
                cnr = ch[cnr][1];
            }
        }
    }

    int kth(int k) {  //第k名的权值
        int cnr = rt;
        while (1) {
            if (ch[cnr][0] && k <= sz[ch[cnr][0]]) {
                cnr = ch[cnr][0];
            } else {
                k -= cnt[cnr] + sz[ch[cnr][0]];
                if (k <= 0) {
                    splay(cnr);
                    return val[cnr];
                }
                cnr = ch[cnr][1];
            }
        }
    }

    int pre() {  //前驱节点编号
        int cnr = ch[rt][0];
        while (ch[cnr][1]) cnr = ch[cnr][1];
        splay(cnr);
        return cnr;
    }  // 若需要得到前驱 tree.ins(x), printf("%d\n", tree.val[tree.pre()]),
       // tree.del(x);

    int nxt() {  //后驱节点编号
        int cnr = ch[rt][1];
        while (ch[cnr][0]) cnr = ch[cnr][0];
        splay(cnr);
        return cnr;
    }  // 若需要得到后驱 tree.ins(x), printf("%d\n", tree.val[tree.pre()]),
       // tree.del(x);

    void del(int k) {  //删除k值
        rk(k);
        if (cnt[rt] > 1) {
            cnt[rt]--;
            maintain(rt);
            return;
        }
        if (!ch[rt][0] && !ch[rt][1]) {
            clear(rt);
            rt = 0;
            return;
        }
        if (!ch[rt][0]) {
            int cnr = rt;
            rt = ch[rt][1];
            fa[rt] = 0;
            clear(cnr);
            return;
        }
        if (!ch[rt][1]) {
            int cnr = rt;
            rt = ch[rt][0];
            fa[rt] = 0;
            clear(cnr);
            return;
        }
        int cnr = rt;
        int x = pre();
        splay(x);
        fa[ch[cnr][1]] = x;
        ch[x][1] = ch[cnr][1];
        clear(cnr);
        maintain(rt);
    }
} tree;
