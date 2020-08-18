//

const int maxbit = 62;		//maxbit不能太大

struct L_B{
	ll lba[maxbit];
	L_B(){
        memset(lba, 0, sizeof(lba));
    }
    
	void Insert(ll val){		//插入
        for(int i = maxbit - 1; i >= 0; -- i) // 从高位向低位扫  
            if(val & (1ll << i)){ // 
                if(!lba[i]){
                    lba[i] = val;
                    break;
                }
                val ^= lba[i];
            }
    }
};
//对原集合的每个数val转为2进制，从高位向低位扫，对于当前位为1的，若lba[i]不存在就令lba[i]=x，否则令val=val`xor`lba[i]
// --------------线性基模板
