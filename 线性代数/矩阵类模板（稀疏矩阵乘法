struct Matrix{
    int n,m;
    int a[maxn][maxn];////
    void clear(){
        n=m=0;
        memset(a,0,sizeof(a));
    }
    Matrix operator * (const Matrix &b) const{
        Matrix tmp;
        tmp.clear();
        tmp.n=n;tmp.m=b.m;
        for (int k=0;k<m;++k){
            for (int i=0;i<n;++i){
            	if(a[i][k]==0) continue;
            	for(int j=0;j<b.m;++j){
            		if(b.a[k][j]==0) continue;
            		tmp.a[i][j]+=a[i][k]*b.a[k][j];
                    tmp.a[i][j]%=mod;
				}       
			}         
        }
        return tmp;
    }
};
//稀疏矩阵乘法 
