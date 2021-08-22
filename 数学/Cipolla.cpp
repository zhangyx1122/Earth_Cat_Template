#include <bits/stdc++.h>

using namespace std;
typedef long long ll;

ll mod;
ll I_mul_I; // 虚数单位的平方

struct Complex {
    ll real, imag;

    Complex(ll real = 0, ll imag = 0) : real(real), imag(imag) {}
};

inline bool operator==(Complex x, Complex y) {
    return x.real == y.real and x.imag == y.imag;
}

inline Complex operator*(Complex x, Complex y) {
    return Complex((x.real * y.real + I_mul_I * x.imag % mod * y.imag) % mod,
                   (x.imag * y.real + x.real * y.imag) % mod);
}

Complex power(Complex x, ll k) {
    Complex res = 1;
    while (k) {
        if (k & 1) res = res * x;
        x = x * x;
        k >>= 1;
    }
    return res;
}

bool check_if_residue(ll x) {
    return power(x, (mod - 1) >> 1) == 1;
}

void solve(ll n, ll &x0, ll &x1) {

    ll a = rand() % mod;
    while (!a or check_if_residue((a * a + mod - n) % mod))
        a = rand() % mod;
    I_mul_I = (a * a + mod - n) % mod;
    x0 = ll(power(Complex(a, 1), (mod + 1) >> 1).real);
    x1 = mod - x0;
}

signed main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    cout.tie(nullptr);

    ll t;
    cin >> t;
    while (t--) {
        ll n;
        cin >> n >> mod;
        if (n == 0) {
            cout << 0 << endl;
            continue;
        }
        if (!check_if_residue(n)) {
            cout << "Hola!" << endl;
            continue;
        }
        ll x0, x1;
        solve(n, x0, x1);
        if (x0 > x1) swap(x0, x1);
        cout << x0 << ' ' << x1 << endl;
    }
}
