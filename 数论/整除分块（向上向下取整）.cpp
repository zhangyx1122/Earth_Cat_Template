int x;
scanf("%d",&x);
int ans1=0,ans2=0;
//向下取整
for(int l=1,r;l<=x;l=r+1){
    int m=x/l;
    r=x/m;
    ans1+=(r-l+1)*m;
}
//向上取整
int R=1e5;
for(int l=1,r;l<=R;l=r+1){
    int m=(x+l-1)/l;
    r=m!=1?(x-1)/(m-1):R;
    ans2+=(r-l+1)*m;
}
