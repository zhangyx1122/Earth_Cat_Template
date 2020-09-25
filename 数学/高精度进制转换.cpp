#include <bits/stdc++.h>
using namespace std;
//将字符串表示的10进制大整数转换为m进制的大整数
//并返回m进制大整数的字符串
bool judge(string s)  //判断串是否为全零串
{
    for (int i = 0; i < s.size(); i++)
        if (s[i] != '0') return 1;
    return 0;
}
string solve(
    string s, int n,
    int m)  // n进制转m进制只限0-9进制，若涉及带字母的进制，稍作修改即可
{
    string r, ans;
    int d = 0;
    if (!judge(s)) return "0";  //特判
    while (judge(s))            //被除数不为0则继续
    {
        for (int i = 0; i < s.size(); i++) {
            r += (d * n + s[i] - '0') / m + '0';  //求出商
            d = (d * n + (s[i] - '0')) % m;       //求出余数
        }
        s = r;           //把商赋给下一次的被除数
        r = "";          //把商清空
        ans += d + '0';  //加上进制转换后数字
        d = 0;           //清空余数
    }
    reverse(ans.begin(), ans.end());  //倒置下
    return ans;
}

//o(n^2)
