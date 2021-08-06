# prufer

 Prufer 序列 (Prufer code)，这是一种将带标号的树用一个唯一的整数序列表示的方法。

Prufer 序列可以将一个带标号n个结点的树用$[1, n]$中的$n - 2$个整数表示。你也可以把它理解为完全图的生成树与数列之间的双射。

显然你不会想不开拿这玩意儿去维护树结构。这玩意儿常用组合计数问题上。

## 线性建立prufer

Prufer 是这样建立的：每次选择一个编号最小的叶结点并删掉它，然后在序列中记录下它连接到的那个结点。重复n - 2次后就只剩下两个结点，算法结束。

线性构造的本质就是维护一个指针指向我们将要删除的结点。首先发现，叶结点数是非严格单调递减的。要么删一个，要么删一个得一个。

于是我们考虑这样一个过程：维护一个指针p 。初始时 p指向编号最小的叶结点。同时我们维护每个结点的度数，方便我们知道在删除结点的时侯是否产生新的叶结点。操作如下：

1. 删除 指向的结点，并检查是否产生新的叶结点。
2. 如果产生新的叶结点，假设编号为x ，我们比较 p, x的大小关系。如果 x>p，那么不做其他操作；否则就立刻删除 x，然后检查删除x 后是否产生新的叶结点，重复 2步骤，直到未产生新节点或者新节点的编号>p 。
3. 让指针 p  自增直到遇到一个未被删除叶结点为止；

循环上述操作n - 2 次，就完成了序列的构造。

```cpp
// 从原文摘的代码，同样以 0 为起点
vector<vector<int>> adj;
vector<int> parent;

void dfs(int v) {
  for (int u : adj[v]) {
    if (u != parent[v]) parent[u] = v, dfs(u);
  }
}

vector<int> pruefer_code() {
  int n = adj.size();
  parent.resize(n), parent[n - 1] = -1;
  dfs(n - 1);

  int ptr = -1;
  vector<int> degree(n);
  for (int i = 0; i < n; i++) {
    degree[i] = adj[i].size();
    if (degree[i] == 1 && ptr == -1) ptr = i;
  }

  vector<int> code(n - 2);
  int leaf = ptr;
  for (int i = 0; i < n - 2; i++) {
    int next = parent[leaf];
    code[i] = next;
    if (--degree[next] == 1 && next < ptr) {
      leaf = next;
    } else {
      ptr++;
      while (degree[ptr] != 1) ptr++;
      leaf = ptr;
    }
  }
  return code;
}
```

## 性质

1. 在构造完 Prufer 序列后原树中会剩下两个结点，其中一个一定是编号最大的点 。
2. 每个结点在序列中出现的次数是其度数减1 。（没有出现的就是叶结点）

## 线性prufer转化成树

同线性构造 Prufer 序列的方法。在删度数的时侯会产生新的叶结点，于是判断这个叶结点与指针p的大小关系，如果更小就优先考虑它

```cpp
// 原文摘代码
vector<pair<int, int>> pruefer_decode(vector<int> const& code) {
  int n = code.size() + 2;
  vector<int> degree(n, 1);
  for (int i : code) degree[i]++;

  int ptr = 0;
  while (degree[ptr] != 1) ptr++;
  int leaf = ptr;

  vector<pair<int, int>> edges;
  for (int v : code) {
    edges.emplace_back(leaf, v);
    if (--degree[v] == 1 && v < ptr) {
      leaf = v;
    } else {
      ptr++;
      while (degree[ptr] != 1) ptr++;
      leaf = ptr;
    }
  }
  edges.emplace_back(leaf, n - 1);
  return edges;
}
```

## cayley公式

完全图$K_n$ 有$n^{n - 2}$ 棵生成树。

用 Prufer 序列证:任意一个长度为n - 2的值域 [1, n] 的整数序列都可以通过 Prufer 序列双射对应一个生成树，于是方案数就是$n^{n - 2}$  。

## 图连通方案数

一个n个点m条边的带标号无向图有k个连通块。我们希望添加k - 1条边使得整个图连通。求方案数。

设$s_i$表示每个连通块的数量。我们对k个连通块构造 Prufer 序列，然后你发现这并不是普通的 Prufer 序列。因为每个连通块的连接方法很多。不能直接淦就设啊。于是设$d_i$为第 i个连通块的度数。由于度数之和是边数的两倍，于是$\sum_{i = 1}^{k} d_i = 2k - 2$ 。则对于给定的d 序列构造 Prufer 序列的方案数是
$$
\tbinom{k - 2}{d_1 - 1, d_2 - 1, \dots, d_k - 1} = \frac{(k - 2)!}{(d_1 - 1)!(d_2 - 1)! \dots (d_k - 1)!}
$$
对于第i个连通块，它的连接方式有$s_i^{d_i}$种，因此对于给定d序列使图连通的方案数是
$$
\tbinom{k - 2}{d_1 - 1, d_2 - 1, \dots, d_k - 1} \prod_{i = 1}^{k}s_i^{d_i}
$$
现在我们要枚举d序列，式子变成
$$
\sum_{d_i\geq 1. \sum_{i = 1}^{k} d_i = 2k - 2} \tbinom{k - 2}{d_1 - 1, d2 - 1, \dots ,d_k - 1} \prod_{i = 1}^{k}s_i^{d_i}
$$
根据多元二项式定理
$$
(x_1+\dots+x_m)^{p}= \sum_{c_i \geq 0, \sum_{i = 1}^{m} c_i = p} \tbinom{p}{C_1, C_2, \dots , C_m} \prod_{i= 1}^{m}x_i^{C_i}
$$

对原式换元，设$e_i = d_i - 1$ ，显然有$\sum_{i = 1} ^{k} e_i = k - 2$ 
$$
\Rightarrow \sum_{e_i\ge 0, \sum_{i= 1}^{k} e_i = k - 2} \tbinom{k -2}{e_1, e_2, \dots, e_k} \prod_{i = 1}^{k}s_i ^{e _i+ 1} \\
化简 \Rightarrow (s_1 + s_2 + \dots + s_k)^{k - 2} \prod_{i = 1}^{k}s_i \\
\Rightarrow n^{k - 2} \prod_{i = 1}^{k} s_i
$$
