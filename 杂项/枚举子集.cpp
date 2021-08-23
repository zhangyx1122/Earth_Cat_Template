  cin >> n;
  for (int s = n; s; s = (s - 1) & n) {
      cout << bitset<8>(s) << endl;
  }
