#include <bits/stdc++.h>
using namespace std;

char next_char() {
	static char buf[1 << 20], *first, *last;
	if(first == last) {
		last = buf + fread(buf, 1, 1 << 20, stdin);
		first = buf;
	}
	return first == last ? EOF : *first ++;
}

inline int read(){
	int x = 0, w = 0; char ch = 0;
	while(!isdigit(ch)) {w |= ch == '-'; ch = next_char(); }
	while(isdigit(ch)) {x = (x << 3) + (x << 1) + (ch ^ 48), ch = next_char(); }
	return w ? -x : x;
}

int main(){
	freopen("1.txt", "r", stdin); // 交代码的时候一定要去掉aaa 
	int T;
	cin >> T;
	while(T --){
		int x = read();
		cout << x << endl;
	}
} 
