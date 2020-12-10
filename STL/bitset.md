[C++ bitset 用法 ](https://www.cnblogs.com/magisk/p/8809922.html)

## bitset

> C++的 bitset 在 bitset 头文件中，它是一种类似数组的结构，它的每一个元素只能是０或１，每个元素仅用１bit空间。
> **bitset数组与vector数组区别**
> bitset声明数组:bitset<100> number[10]
> vector声明数组:vector number[10];
> **bitset<每个bitset元素的长度(没有占满前面全部自动补0)> 元素**
> **bitset内置转化函数：可将bitset转化为string,unsigned long,unsigned long long。**

####构造

```cpp
	bitset<4> bitset1;　　//无参构造，长度为４，默认每一位为０

    bitset<8> bitset2(12);　　//长度为８，二进制保存，前面用０补充

    string s = "100101";
    bitset<10> bitset3(s);　　//长度为10，前面用０补充
    
    char s2[] = "10101";
    bitset<13> bitset4(s2);　　//长度为13，前面用０补充

    cout << bitset1 << endl;　　//0000
    cout << bitset2 << endl;　　//00001100
    cout << bitset3 << endl;　　//0000100101
    cout << bitset4 << endl;　　//0000000010101
```

#### 函数

```cpp
	bitset<8> foo ("10011011");

    cout << foo.count() << endl;　　//5　　（count函数用来求bitset中1的位数，foo中共有５个１
    cout << foo.size() << endl;　　 //8　　（size函数用来求bitset的大小，一共有８位

    cout << foo.test(0) << endl;　　//true　　（test函数用来查下标处的元素是０还是１，并返回false或true，此处foo[0]为１，返回true
    cout << foo.test(2) << endl;　　//false　　（同理，foo[2]为０，返回false

    cout << foo.any() << endl;　　//true　　（any函数检查bitset中是否有１
    cout << foo.none() << endl;　　//false　　（none函数检查bitset中是否没有１
    cout << foo.all() << endl;　　//false　　（all函数检查bitset中是全部为１
```



[2019-2020 ICPC Asia Taipei-Hsinchu Regional Contest（H](https://blog.csdn.net/chitudexixi/article/details/109453360)

###H

```cpp
#include <bits/stdc++.h>
#define ll long long
using namespace std;
int t,n,m;
char str[1010];
bitset<500> number[30];
int main() {
	ios::sync_with_stdio(false); cin.tie(0); cout.tie(0);
    //freopen("test.in","r",stdin);
    //freopen("test.out","w",stdout);
	scanf("%d",&t);
	while(t--)
	{
		scanf("%d %d",&n,&m);
		for(int i=0;i<m;i++)
		{
			scanf("%s",str);
			number[i]=bitset<500>(str);
		}
		int len=1<<m,ans=m+1;
		for(int i=1;i<len;i++)
		{
			int t=i,s=0;
			bitset<500> num(0);
			for(int j=0;j<m&&t>0;j++)
			{
				if(t&1) 
				{
					num=num|number[j];
					s++;
				}
				t>>=1;
			}
			if(num.count()==n) ans=min(ans,s);
		}
		if(ans==m+1) printf("-1\n");
		else printf("%d\n",ans);
	}
	return 0;
}

```

