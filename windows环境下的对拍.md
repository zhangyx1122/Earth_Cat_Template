

```
@echo off
:loop
	dataa.exe > data.txt
	biaocheng.exe < data.txt > ac.txt
	A.exe < data.txt > test.txt
	fc ac.txt test.txt
	if not errorlevel 1 goto loop
pause
goto loop
```

**其中要改的部分（标红辽）**：

@echo off
:loop
	dataa.exe > data.txt
	$\color{red}{biaocheng.exe}$ < data.txt > ac.txt
	$\color{red}{A.exe}$ < data.txt > test.txt
	fc ac.txt test.txt
	if not errorlevel 1 goto loop
pause
goto loop



文件以`.bat`作为后缀

---

将三个程序（数据生成文件（dataa），标程或暴力代码（biaocheng）, 要看的代码（A））放在同一目录下，

记得加 `freopen`

随机数记得加`srand((int)time(0));`

---

随机数生成code

```c++
#include <iostream>
#include <cstdlib>
#include <ctime>
using namespace std;

int main(){
	freopen("data.txt", "w", stdout);
	
    srand((int)time(0));
    int T = rand() % 100000;
    cout << T << endl;
	 
    for (int i = 0; i < T; i++){
    	cout << rand() % 100;
    }
}
```



`rand()` 似乎只有三万多，需要更大的数的话要乘一下