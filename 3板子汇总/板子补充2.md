## 杂项

### fread快读

```c++
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
```

---

## STL

### set

**begin()** ,返回set容器的第一个元素

**end()**  ,返回set容器的最后一个元素

**clear()**   ,删除set容器中的所有的元素

**empty()**,判断set容器是否为空

**max_size()** ,返回set容器可能包含的元素最大个数

**size()**  ,返回当前set容器中的元素个数

**rbegin** ,返回的值和end()相同

**rend()**,返回的值和rbegin()相同

**count()** 用来查找set中某个某个键值出现的次数。

**equal_range()** ，返回一对定位器，分别表示第一个大于或等于给定关键值的元素和 第一个大于给定关键值的元素，这个返回值是一个pair类型，如果这一对定位器中哪个返回失败，就会等于end()的值。

**erase(iterator)** ,删除定位器iterator指向的值

**erase(first,second)**,删除定位器first和second之间的值

**erase(key_value)**,删除键值key_value的值

**find()** ，返回给定值值得定位器，如果没找到则返回end()。

**insert(key_value);** 将key_value插入到set中 ，返回值是pair<set<int>::iterator,bool>，bool标志着插入是否成功，而iterator代表插入的位置，若key_value已经在set中，则iterator表示的key_value在set中的位置。

**inset(first,second);**将定位器first到second之间的元素插入到set中，返回值是void.

**lower_bound(key_value)** ，返回第一个大于等于key_value的定位器

**upper_bound(key_value)，**返回最后一个大于等于key_value的定位器

---

### map

**插入操作**

使用[ ]进行单个插入

```c++
map<int, string> ID_Name;

// 如果已经存在键值2015，则会作赋值修改操作，如果没有则插入
ID_Name[2015] = "Tom";1234
```

使用insert进行单个和多个插入 (insert共有4个重载函数：

```c++
// 插入单个键值对，并返回插入位置和成功标志，插入位置已经存在值时，插入失败
pair<iterator,bool> insert (const value_type& val);

//在指定位置插入，在不同位置插入效率是不一样的，因为涉及到重排
iterator insert (const_iterator position, const value_type& val);

// 插入多个
void insert (InputIterator first, InputIterator last);

//c++11开始支持，使用列表插入多个   
void insert (initializer_list<value_type> il);
```

**取值**

Map中元素取值主要有at和[ ]两种操作，at会作下标检查，而[]不会。

```c++
map<int, string> ID_Name;

//ID_Name中没有关键字2016，使用[]取值会导致插入
//因此，下面语句不会报错，但打印结果为空
cout<<ID_Name[2016].c_str()<<endl;

//使用at会进行关键字检查，因此下面语句会报错
ID_Name.at(2016) = "Bob";
```

**容量查询**

```c++
// 查询map是否为空
bool empty();

// 查询map中键值对的数量
size_t size();

// 查询map所能包含的最大键值对数量，和系统和应用库有关。
// 此外，这并不意味着用户一定可以存这么多，很可能还没达到就已经开辟内存失败了
size_t max_size();

// 查询关键字为key的元素的个数，在map里结果非0即1
size_t count( const Key& key ) const; //
```

**迭代器**

共有八个获取迭代器的函数：**begin, end, rbegin,rend** 以及对应的 **cbegin, cend, crbegin,crend**。

二者的区别在于，后者一定返回 const_iterator，而前者则根据map的类型返回iterator 或者 const_iterator。const情况下，不允许对值进行修改。如下面代码所示：

```c++
map<int,int>::iterator it;
map<int,int> mmap;
const map<int,int> const_mmap;

it = mmap.begin(); //iterator
mmap.cbegin(); //const_iterator

const_mmap.begin(); //const_iterator
const_mmap.cbegin(); //const_iterator123456789
```

返回的迭代器可以进行加减操作，此外，如果map为空，则 begin = end。

 **删除**

```c++
// 删除迭代器指向位置的键值对，并返回一个指向下一元素的迭代器
iterator erase( iterator pos )

// 删除一定范围内的元素，并返回一个指向下一元素的迭代器
iterator erase( const_iterator first, const_iterator last );

// 根据Key来进行删除， 返回删除的元素数量，在map里结果非0即1
size_t erase( const key_type& key );

// 清空map，清空后的size为0
void clear();
```

**交换**

```c++
// 就是两个map的内容互换
void swap( map& other );
```

**顺序比较**

```c++
// 比较两个关键字在map中位置的先后
key_compare key_comp() const;
```

**查找**

```c++
// 关键字查询，找到则返回指向该关键字的迭代器，否则返回指向end的迭代器
// 根据map的类型，返回的迭代器为 iterator 或者 const_iterator
iterator find (const key_type& k);
const_iterator find (const key_type& k) const;
```

**操作符**

operator: == != < <= > >=
**注意** 对于==运算符, 只有键值对以及顺序完全相等才算成立。

---

### unordered_map

**查找元素是否存在**

若有unordered_map <int, int> mp;查找x是否在map中
    方法1:  若存在  mp.find(x)!=mp.end()
    方法2:  若存在  mp.count(x)!=0123

**插入数据**

```c++
mp.insert(Map::value_type(1,"Raoul"));1
```

**遍历map**

```c++
 unordered_map<key,T>::iterator it;
    (*it).first;   //the key value
    (*it).second   //the mapped value
    for(unordered_map<key,T>::iterator iter=mp.begin();iter!=mp.end();iter++)
          cout<<"key value is"<<iter->first<<" the mapped value is "<< iter->second;

    //也可以这样
    for(auto& v : mp)
        print v.first and v.second
```

---

## 数论

### lucas求组合数

```c++
#include <bits/stdc++.h>
using namespace std;

typedef long long ll;

ll p;

const int maxn = 1e5 + 10;

ll qpow(ll x, ll n){
	ll res = 1;
	while(n){
		if(n & 1) res = (res * x) % p;
		x = (x * x) % p;
		n >>= 1;
	}
	
	return res;
}

ll C(ll up, ll down){
	if(up > down) return 0;
	ll res = 1;

//	for(int i = up + 1; i <= down; ++ i){
//		res = (res * i) % p;
//	}
//	for(int i = 1; i <= down - up; ++ i){
//		res = (res * qpow(i, p - 2)) % p; 
//	}

	for(int i = 1, j = down; i <= up; ++ i, -- j){
		res = (res * j) % p;
		res = (res * qpow(i, p - 2)) % p;
	}
	
	return res;
}


ll lucas(ll up, ll down){
	if(up < p && down < p) return C(up, down);
	return C(up % p, down % p) * lucas(up / p, down / p) % p; 
}

int main(){
	ios::sync_with_stdio(0); cin.tie(0); cout.tie(0);
	
	int T;
	cin >> T;
	while (T --){
		ll down, up;
		cin >> down >> up >> p;
		
		cout << lucas(up, down) % p << endl;
	}
	
	return 0;
} 
```

---

### 扩展欧几里得求逆元

```c++
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
```

