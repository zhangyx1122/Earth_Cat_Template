#include <iostream>
#include <cstring>
#include <algorithm>

using namespace std;

struct Dinic {
	static constexpr int N = 10010, M = 100010, INF = 1e8;
	static constexpr double eps = 1e-8;	
//	int n, m, S, T;
	int S, T;
	int h[N], e[M], ne[M], idx;
	double f[M];
	int q[N], d[N], cur[N]; // d 表示从源点开始走到该点的路径上所有边的容量的最小值 
	
	void AddEdge(int a, int b, double c)
	{
	    e[idx] = b, f[idx] = c, ne[idx] = h[a], h[a] = idx ++ ;
	    e[idx] = a, f[idx] = 0, ne[idx] = h[b], h[b] = idx ++ ;
	}
	
	bool bfs()
	{
	    int hh = 0, tt = 0;
	    memset(d, -1, sizeof d);
	    q[0] = S, d[S] = 0, cur[S] = h[S];
	    while (hh <= tt)
	    {
	        int t = q[hh ++ ];
	        for (int i = h[t]; ~i; i = ne[i])
	        {
	            int ver = e[i];
	            if (d[ver] == -1 && f[i] > 0)
	            {
	                d[ver] = d[t] + 1;
	                cur[ver] = h[ver];
	                if (ver == T) return true;
	                q[ ++ tt] = ver;
	            }
	        }
	    }
	    return false;
	}
	
	double find(int u, double limit)
	{
	    if (u == T) return limit;
	    double flow = 0;
	    for (int i = cur[u]; ~i && flow < limit; i = ne[i])
	    {
	        cur[u] = i;
	        int ver = e[i];
	        if (d[ver] == d[u] + 1 && f[i] > 0)
	        {
	            double t = find(ver, min(f[i], limit - flow));
	            if (t < eps) d[ver] = -1;
	            f[i] -= t, f[i ^ 1] += t, flow += t;
	        }
	    }
	    return flow;
	}
	
	double Maxflow(int S, int T)
	{
		this->S = S, this->T = T;
	    double r = 0, flow;
	    while (bfs()) while (flow = find(S, INF)) r += flow;
	    return r;
	}		
	void init() //////// 
	{
		memset(h, -1, sizeof h);
		idx = 0; 
	}
} MF;

// 先init 

