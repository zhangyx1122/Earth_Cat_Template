#include <bits/stdc++.h>

using namespace std;
typedef long long ll;
const ll N = 1 << 20;
const ll mod = 1e9 + 7;
const double dinf = 1e99;
const int inf = 0x3f3f3f3f;
const ll linf = 0x3f3f3f3f3f3f3f3f;

const double eps = 1e-9;
const double PI = acos(-1.0);

struct Line;

struct Point {
    double x, y;

    Point() { x = y = 0; }

    Point(const Line &a);

    Point(const double &a, const double &b) : x(a), y(b) {}

    Point operator+(const Point &a) const {
        return {x + a.x, y + a.y};
    }

    Point operator-(const Point &a) const {
        return {x - a.x, y - a.y};
    }

    Point operator*(const double &a) const {
        return {x * a, y * a};
    }

    Point operator/(const double &d) const {
        return {x / d, y / d};
    }

    bool operator==(const Point &a) const {
        return abs(x - a.x) + abs(y - a.y) < eps;
    }

    void standardize() {
        *this = *this / sqrt(x * x + y * y);
    }
};

Point normal(const Point &a) { return Point(-a.y, a.x); }

double dist(const Point &a, const Point &b) {
    return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}

double dist2(const Point &a, const Point &b) {
    return (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y);
}

struct Line {
    Point s, t;

    Line() {}

    Line(const Point &a, const Point &b) : s(a), t(b) {}

};

struct circle {
    Point o;
    double r;

    circle() {}

    circle(Point P, double R = 0) { o = P, r = R; }
};

double length(const Point &p) {
    return sqrt(p.x * p.x + p.y * p.y);
}

double length(const Line &l) {
    Point p(l);
    return length(p);
}

Point::Point(const Line &a) { *this = a.t - a.s; }

istream &operator>>(istream &in, Point &a) {
    in >> a.x >> a.y;
    return in;
}

double dot(const Point &a, const Point &b) {
    return a.x * b.x + a.y * b.y;
}

double det(const Point &a, const Point &b) {
    return a.x * b.y - a.y * b.x;
}

int sgn(const double &x) { return fabs(x) < eps ? 0 : (x > 0 ? 1 : -1); }

double sqr(const double &x) { return x * x; }

Point rotate(const Point &a, const double &ang) {
    double x = cos(ang) * a.x - sin(ang) * a.y;
    double y = sin(ang) * a.x + cos(ang) * a.y;
    return {x, y};
}

//点在线段上 <=0 包含端点
bool sp_on(const Line &seg, const Point &p) {
    Point a = seg.s, b = seg.t;
    return !sgn(det(p - a, b - a)) && sgn(dot(p - a, p - b)) <= 0;
}

bool lp_on(const Line &line, const Point &p) {
    Point a = line.s, b = line.t;
    return !sgn(det(p - a, b - a));
}

//等于不包含共线
int andrew(Point *point, Point *convex, int n) {
    sort(point, point + n, [](Point a, Point b) {
        if (a.x != b.x) return a.x < b.x;
        return a.y < b.y;
    });
    int top = 0;
    for (int i = 0; i < n; i++) {
        while ((top > 1) && det(convex[top - 1] - convex[top - 2], point[i] - convex[top - 1]) <= 0)
            top--;
        convex[top++] = point[i];
    }
    int tmp = top;
    for (int i = n - 2; i >= 0; i--) {
        while ((top > tmp) && det(convex[top - 1] - convex[top - 2], point[i] - convex[top - 1]) <= 0)
            top--;
        convex[top++] = point[i];
    }
    if (n > 1) top--;
    return top;
}

double slope(const Point &a, const Point &b) {
    return (a.y - b.y) / (a.x - b.x);
}

double slope(const Line &a) {
    return slope(a.s, a.t);
}

Point ll_intersection(const Line &a, const Line &b) {
    double s1 = det(Point(a), b.s - a.s), s2 = det(Point(a), b.t - a.s);
    return (b.s * s2 - b.t * s1) / (s2 - s1);
}

int ss_cross(const Line &a, const Line &b, Point &p) {
    int d1 = sgn(det(a.t - a.s, b.s - a.s));
    int d2 = sgn(det(a.t - a.s, b.t - a.s));
    int d3 = sgn(det(b.t - b.s, a.s - b.s));
    int d4 = sgn(det(b.t - b.s, a.t - b.s));
    if ((d1 ^ d2) == -2 && (d3 ^ d4) == -2) {
        p = ll_intersection(a, b);
        return 1;
    }
    if (!d1 && sp_on(a, b.s)) {
        p = b.s;
        return 2;
    }
    if (!d2 && sp_on(a, b.t)) {
        p = b.t;
        return 2;
    }
    if (!d3 && sp_on(b, a.s)) {
        p = a.s;
        return 2;
    }
    if (!d4 && sp_on(b, a.t)) {
        p = a.t;
        return 2;
    }
    return 0;
}

Point project(const Line &l, const Point &p) {
    Point base(l);
    double r = dot(base, p - l.s) / sqr(length(base));
    return l.s + (base * r);
}

double sp_dist(const Line &l, const Point &p) {
    if (l.s == l.t) return dist(l.s, p);
    Point x = p - l.s, y = p - l.t, z = l.t - l.s;
    if (sgn(dot(x, z)) < 0)return length(x);//P距离A更近
    if (sgn(dot(y, z)) > 0)return length(y);//P距离B更近
    return abs(det(x, z) / length(z));//面积除以底边长
}

double lp_dist(const Line &l, const Point &p) {
    Point x = p - l.s, y = p - l.t, z = l.t - l.s;
    return abs(det(x, z) / length(z));//面积除以底边长
}

int lc_cross(const Line &l, const Point &a, const double &r, pair<Point, Point> &ans) {
    int num = 0;
    Point pr = project(l, a);
    double dis = dist(pr, a);
    double tmp = r * r - dis * dis;
    if (sgn(tmp) == 1) num = 2;
    else if (sgn(tmp) == 0) num = 1;
    else return 0;
    double base = sqrt(r * r - dis * dis);
    Point e(l);
    e.standardize();
    e = e * base;
    ans = make_pair(pr + e, pr - e);
    return num;
}

int cc_cross(const Point &c1, const double &r1, const Point &c2, const double &r2, pair<Point, Point> &ans) {
    double x1 = c1.x, x2 = c2.x, y1 = c1.y, y2 = c2.y;
    double d = length(c1 - c2);
    if (sgn(fabs(r1 - r2) - d) > 0) return -1;  //内含
    if (sgn(r1 + r2 - d) < 0) return 0; //相离
    double a = r1 * (x1 - x2) * 2, b = r1 * (y1 - y2) * 2, c = r2 * r2 - r1 * r1 - d * d;
    double p = a * a + b * b, q = -a * c * 2, r = c * c - b * b;

    double cosa, sina, cosb, sinb;
    //One Intersection
    if (sgn(d - (r1 + r2)) == 0 || sgn(d - fabs(r1 - r2)) == 0) {
        cosa = -q / p / 2;
        sina = sqrt(1 - sqr(cosa));
        Point p0(x1 + r1 * cosa, y1 + r1 * sina);
        if (sgn(dist(p0, c2) - r2)) p0.y = y1 - r1 * sina;
        ans = pair<Point, Point>(p0, p0);
        return 1;
    }
    //Two Intersections
    double delta = sqrt(q * q - p * r * 4);
    cosa = (delta - q) / p / 2;
    cosb = (-delta - q) / p / 2;
    sina = sqrt(1 - sqr(cosa));
    sinb = sqrt(1 - sqr(cosb));
    Point p1(x1 + r1 * cosa, y1 + r1 * sina);
    Point p2(x1 + r1 * cosb, y1 + r1 * sinb);
    if (sgn(dist(p1, c2) - r2)) p1.y = y1 - r1 * sina;
    if (sgn(dist(p2, c2) - r2)) p2.y = y1 - r1 * sinb;
    if (p1 == p2) p1.y = y1 - r1 * sina;
    ans = pair<Point, Point>(p1, p2);
    return 2;
}

Point lp_sym(const Line &l, const Point &p) {
    return p + (project(l, p) - p) * 2;
}

double alpha(const Point &t1, const Point &t2) {
    double theta;
    theta = atan2((double) t2.y, (double) t2.x) - atan2((double) t1.y, (double) t1.x);
    if (sgn(theta) < 0)
        theta += 2.0 * PI;
    return theta;
}

int pip(const Point *P, const int &n, const Point &a) {//【射线法】判断点A是否在任意多边形Poly以内
    int cnt = 0;
    int tmp;
    for (int i = 1; i <= n; ++i) {
        int j = i < n ? i + 1 : 1;
        if (sp_on(Line(P[i], P[j]), a))return 2;//点在多边形上
        if (a.y >= min(P[i].y, P[j].y) && a.y < max(P[i].y, P[j].y))//纵坐标在该线段两端点之间
            tmp = P[i].x + (a.y - P[i].y) / (P[j].y - P[i].y) * (P[j].x - P[i].x), cnt += sgn(tmp - a.x) > 0;//交点在A右方
    }
    return cnt & 1;//穿过奇数次则在多边形以内
}

bool pip_convex_jud(const Point &a, const Point &L, const Point &R) {//判断AL是否在AR右边
    return sgn(det(L - a, R - a)) > 0;//必须严格以内
}

bool pip_convex(const Point *P, const int &n, const Point &a) {//【二分法】判断点A是否在凸多边形Poly以内
    //点按逆时针给出
    if (pip_convex_jud(P[0], a, P[1]) || pip_convex_jud(P[0], P[n - 1], a)) return 0;//在P[0_1]或P[0_n-1]外
    if (sp_on(Line(P[0], P[1]), a) || sp_on(Line(P[0], P[n - 1]), a)) return 2;//在P[0_1]或P[0_n-1]上
    int l = 1, r = n - 2;
    while (l < r) {//二分找到一个位置pos使得P[0]_A在P[0_pos],P[0_(pos+1)]之间
        int mid = (l + r + 1) >> 1;
        if (pip_convex_jud(P[0], P[mid], a))l = mid;
        else r = mid - 1;
    }
    if (pip_convex_jud(P[l], a, P[l + 1]))return 0;//在P[pos_(pos+1)]外
    if (sp_on(Line(P[l], P[l + 1]), a))return 2;//在P[pos_(pos+1)]上
    return 1;
}
// 多边形是否包含线段
// 因此我们可以先求出所有和线段相交的多边形的顶点，然后按照X-Y坐标排序(X坐标小的排在前面，对于X坐标相同的点，Y坐标小的排在前面，
// 这种排序准则也是为了保证水平和垂直情况的判断正确)，这样相邻的两个点就是在线段上相邻的两交点，如果任意相邻两点的中点也在多边形内，
// 则该线段一定在多边形内。

int pp_judge(Point *A, int n, Point *B, int m) {//【判断多边形A与多边形B是否相离】
    for (int i1 = 1; i1 <= n; ++i1) {
        int j1 = i1 < n ? i1 + 1 : 1;
        for (int i2 = 1; i2 <= m; ++i2) {
            int j2 = i2 < m ? i2 + 1 : 1;
            Point tmp;
            if (ss_cross(Line(A[i1], A[j1]), Line(B[i2], B[j2]), tmp)) return 0;//两线段相交
            if (pip(B, m, A[i1]) || pip(A, n, B[i2]))return 0;//点包含在内
        }
    }
    return 1;
}

double area(Point *P, int n) {//【任意多边形P的面积】
    double S = 0;
    for (int i = 1; i <= n; i++) S += det(P[i], P[i < n ? i + 1 : 1]);
    return S / 2.0;
}

Line Q[N];

int judge(Line L, Point a) { return sgn(det(a - L.s, L.t - L.s)) > 0; }//判断点a是否在直线L的右边
int halfcut(Line *L, int n, Point *P) {//【半平面交】
    sort(L, L + n, [](const Line &a, const Line &b) {
        double d = atan2((a.t - a.s).y, (a.t - a.s).x) - atan2((b.t - b.s).y, (b.t - b.s).x);
        return sgn(d) ? sgn(d) < 0 : judge(a, b.s);
    });

    int m = n;
    n = 0;
    for (int i = 0; i < m; ++i)
        if (i == 0 || sgn(atan2(Point(L[i]).y, Point(L[i]).x) - atan2(Point(L[i - 1]).y, Point(L[i - 1]).x)))
            L[n++] = L[i];
    int h = 1, t = 0;
    for (int i = 0; i < n; ++i) {
        while (h < t && judge(L[i], ll_intersection(Q[t], Q[t - 1]))) --t;//当队尾两个直线交点不是在直线L[i]上或者左边时就出队
        while (h < t && judge(L[i], ll_intersection(Q[h], Q[h + 1]))) ++h;//当队头两个直线交点不是在直线L[i]上或者左边时就出队
        Q[++t] = L[i];

    }
    while (h < t && judge(Q[h], ll_intersection(Q[t], Q[t - 1]))) --t;
    while (h < t && judge(Q[t], ll_intersection(Q[h], Q[h + 1]))) ++h;
    n = 0;
    for (int i = h; i <= t; ++i) {
        P[n++] = ll_intersection(Q[i], Q[i < t ? i + 1 : h]);
    }
    return n;
}

Point V1[N], V2[N];

int mincowski(Point *P1, int n, Point *P2, int m, Point *V) {//【闵可夫斯基和】求两个凸包{P1},{P2}的向量集合{V}={P1+P2}构成的凸包
    for (int i = 0; i < n; ++i) V1[i] = P1[(i + 1) % n] - P1[i];
    for (int i = 0; i < m; ++i) V2[i] = P2[(i + 1) % m] - P2[i];
    int t = 0, i = 0, j = 0;
    V[t++] = P1[0] + P2[0];
    while (i < n && j < m) V[t] = V[t - 1] + (sgn(det(V1[i], V2[j])) > 0 ? V1[i++] : V2[j++]), t++;
    while (i < n) V[t] = V[t - 1] + V1[i++], t++;
    while (j < m) V[t] = V[t - 1] + V2[j++], t++;
    return t;
}

circle getcircle(const Point &A, const Point &B, const Point &C) {//【三点确定一圆】向量垂心法
    Point P1 = (A + B) * 0.5, P2 = (A + C) * 0.5;
    Line R1 = Line(P1, P1 + normal(B - A));
    Line R2 = Line(P2, P2 + normal(C - A));
    circle O;
    O.o = ll_intersection(R1, R2);
    O.r = length(A - O.o);
    return O;
}

struct ConvexHull {

    int op;

    struct cmp {
        bool operator()(const Point &a, const Point &b) const {
            return sgn(a.x - b.x) < 0 || sgn(a.x - b.x) == 0 && sgn(a.y - b.y) < 0;
        }
    };

    set<Point, cmp> s;

    ConvexHull(int o) {
        op = o;
        s.clear();
    }

    inline int PIP(Point P) {
        set<Point>::iterator it = s.lower_bound(Point(P.x, -dinf));//找到第一个横坐标大于P的点
        if (it == s.end())return 0;
        if (sgn(it->x - P.x) == 0) return sgn((P.y - it->y) * op) <= 0;//比较纵坐标大小
        if (it == s.begin())return 0;
        set<Point>::iterator j = it, k = it;
        --j;
        return sgn(det(P - *j, *k - *j) * op) >= 0;//看叉姬1
    }

    inline int judge(set<Point>::iterator it) {
        set<Point>::iterator j = it, k = it;
        if (j == s.begin())return 0;
        --j;
        if (++k == s.end())return 0;
        return sgn(det(*it - *j, *k - *j) * op) >= 0;//看叉姬
    }

    inline void insert(Point P) {
        if (PIP(P))return;//如果点P已经在凸壳上或凸包里就不插入了
        set<Point>::iterator tmp = s.lower_bound(Point(P.x, -inf));
        if (tmp != s.end() && sgn(tmp->x - P.x) == 0)s.erase(tmp);//特判横坐标相等的点要去掉
        s.insert(P);
        set<Point>::iterator it = s.find(P), p = it;
        if (p != s.begin()) {
            --p;
            while (judge(p)) {
                set<Point>::iterator temp = p--;
                s.erase(temp);
            }
        }
        if ((p = ++it) != s.end()) {
            while (judge(p)) {
                set<Point>::iterator temp = p++;
                s.erase(temp);
            }
        }
    }
} up(1), down(-1);

int PIC(circle C, Point a) { return sgn(length(a - C.o) - C.r) <= 0; }//判断点A是否在圆C内
void Random(Point *P, int n) { for (int i = 0; i < n; ++i)swap(P[i], P[(rand() + 1) % n]); }//随机一个排列
circle min_circle(Point *P, int n) {//【求点集P的最小覆盖圆】 O(n)
//  random_shuffle(P,P+n);
    Random(P, n);
    circle C = circle(P[0], 0);
    for (int i = 1; i < n; ++i)
        if (!PIC(C, P[i])) {
            C = circle(P[i], 0);
            for (int j = 0; j < i; ++j)
                if (!PIC(C, P[j])) {
                    C.o = (P[i] + P[j]) * 0.5, C.r = length(P[j] - C.o);
                    for (int k = 0; k < j; ++k) if (!PIC(C, P[k])) C = getcircle(P[i], P[j], P[k]);
                }
        }
    return C;
}
