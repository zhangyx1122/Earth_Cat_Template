const double Pi = acos(-1.0);

struct Complex {
    double x, y;

    Complex(double xx = 0, double yy = 0) { x = xx, y = yy; }
} a[N], b[N];

Complex operator+(Complex _a, Complex _b) { return Complex(_a.x + _b.x, _a.y + _b.y); }

Complex operator-(Complex _a, Complex _b) { return Complex(_a.x - _b.x, _a.y - _b.y); }

Complex operator*(Complex _a, Complex _b) {
    return Complex(_a.x * _b.x - _a.y * _b.y, _a.x * _b.y + _a.y * _b.x);
} //不懂的看复数的运算那部分

int L, r[N];
int limit = 1;

void fft(Complex *A, int type) {
    for (int i = 0; i < limit; i++)
        if (i < r[i]) swap(A[i], A[r[i]]); //求出要迭代的序列
    for (int mid = 1; mid < limit; mid <<= 1) { //待合并区间的长度的一半
        Complex Wn(cos(Pi / mid), type * sin(Pi / mid)); //单位根
        for (int R = mid << 1, j = 0; j < limit; j += R) { //R是区间的长度，j表示前已经到哪个位置了
            Complex w(1, 0); //幂
            for (int k = 0; k < mid; k++, w = w * Wn) { //枚举左半部分
                Complex x = A[j + k], y = w * A[j + mid + k]; //蝴蝶效应
                A[j + k] = x + y;
                A[j + mid + k] = x - y;
            }
        }
    }
}

void FFT(int n, int m) {
    limit = 1;
    L = 0;
    while (limit <= n + m) limit <<= 1, L++;
    for (int i = 0; i < limit; i++) r[i] = (r[i >> 1] >> 1) | ((i & 1) << (L - 1));
    // 在原序列中 i 与 i/2 的关系是 ： i可以看做是i/2的二进制上的每一位左移一位得来
    // 那么在反转后的数组中就需要右移一位，同时特殊处理一下奇数
    fft(a, 1), fft(b, 1);
    for (int i = 0; i <= limit; i++) a[i] = a[i] * b[i];
    fft(a, -1);
    for (int i = 0; i <= n + m; i++) a[i].x /= limit;
}
