#include <bits/stdc++.h>
using namespace std;
string div(string a, int b)  //高精度a除以单精度b
{
    string r, ans;
    int d = 0;
    if (a == "0") return a;  //特判
    for (int i = 0; i < a.size(); i++) {
        r += (d * 10 + a[i] - '0') / b + '0';  //求出商
        d = (d * 10 + (a[i] - '0')) % b;       //求出余数
    }
    int p = 0;
    for (int i = 0; i < r.size(); i++)
        if (r[i] != '0') {
            p = i;
            break;
        }
    return r.substr(p);
}

//o(n)
