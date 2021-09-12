#include<bits/stdc++.h>
#define fi first
#define sf scanf
#define se second
#define pf printf
#define pb push_back
#define mp make_pair
#define sz(x) ((int)(x).size())
#define all(x) (x).begin(),(x).end()
#define mem(x,y) memset((x),(y),sizeof(x))
#define fup(i,x,y) for(int i=(x);i<=(y);++i)
#define fdn(i,x,y) for(int i=(x);i>=(y);--i)
typedef long long ll;
typedef long double ld;
typedef unsigned long long ull;
typedef std::pair<int,int> pii;
using namespace std;
 
const ld pi=acos(-1);
 
ld pow2(ld x){return x*x;}
 
ld pow3(ld x){return x*x*x;}
 
ld dis(ld x1,ld y1,ld z1,ld x2,ld y2,ld z2)
{
    return pow2(x1-x2)+pow2(y1-y2)+pow2(z1-z2);
}
 
ld cos(ld a,ld b,ld c){return (b*b+c*c-a*a)/(2*b*c);}
 
ld cap(ld r,ld h){return pi*(r*3-h)*h*h/3;} // 球缺体积公式，h为球缺的高
 
//2球体积交
ld sphere_intersect(ld x1,ld y1,ld z1,ld r1,ld x2,ld y2,ld z2,ld r2)
{
    ld d=dis(x1,y1,z1,x2,y2,z2);
    //相离
    if(d>=pow2(r1+r2))return 0;
    //包含
    if(d<=pow2(r1-r2))return pow3(min(r1,r2))*4*pi/3;
    //相交
    ld h1=r1-r1*cos(r2,r1,sqrt(d)),h2=r2-r2*cos(r1,r2,sqrt(d));
    return cap(r1,h1)+cap(r2,h2);
}
 
//2球体积并
ld sphere_union(ld x1,ld y1,ld z1,ld r1,ld x2,ld y2,ld z2,ld r2)
{
    ld d=dis(x1,y1,z1,x2,y2,z2);
    //相离
    if(d>=pow2(r1+r2))return (pow3(r1)+pow3(r2))*4*pi/3;
    //包含
    if(d<=pow2(r1-r2))return pow3(max(r1,r2))*4*pi/3;
    //相交
    ld h1=r1+r1*cos(r2,r1,sqrt(d)),h2=r2+r2*cos(r1,r2,sqrt(d));
    return cap(r1,h1)+cap(r2,h2);
}
 
int main()
{
    double x1,y1,z1,r1,x2,y2,z2,r2;
    sf("%lf%lf%lf%lf%lf%lf%lf%lf",&x1,&y1,&z1,&r1,&x2,&y2,&z2,&r2);
    pf("%.12Lf\n",sphere_union(x1,y1,z1,r1,x2,y2,z2,r2));
    return 0;
}
