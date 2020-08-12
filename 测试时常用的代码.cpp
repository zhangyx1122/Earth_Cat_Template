//文件输入输出、时间函数与随机数

#ifdef ONLINE_JUDGE
#else
    freopen("in.txt","r",stdin);
    //freopen("out.txt","w",stdout);
#endif
//重新定向stdin/stdout到in.txt/out.txt，"r"和"w"为只读和只写


#include<ctime>
    clock_t ST,ED;
    ST=clock();
    //这里填测试的程序
    ED=clock();
    cout<<ED-ST<<"ms"<<endl;


#include<ctime>
#include<cstdlib>
    srand(time(0));//初始化
    rand();//返回[0,RAND_MAX]之间的随机整数(int)，RAND_MAX是cstdlib中的宏定义，一般为0x7fff(32767)