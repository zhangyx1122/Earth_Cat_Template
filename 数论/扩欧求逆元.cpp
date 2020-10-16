#include <bits/stdc++.h>
using namespace std;

typedef  long long ll;

void extgcd(ll a,ll b,ll& d,ll& x,ll& y){
    if(!b){ d=a; x=1; y=0;}
    else{ extgcd(b,a%b,d,y,x); y-=x*(a/b); }
}

ll inverse(ll a,ll n){
    ll d,x,y;
    extgcd(a,n,d,x,y);
    return d==1?(x+n)%n:-1;
}

int main(){
	int x, y;
	//cin >> x >> y;
	while(1){
		cin >> x >> y;
		cout << inverse(x, y) << endl;
	} 
	//cout << inverse(x, y) << endl;
} 
