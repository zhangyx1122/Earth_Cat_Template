#define FI(n) FastIO::read(n)
#define FO(n) FastIO::write(n)
#define Flush FastIO::Fflush()
//程序末尾写上   Flush;

namespace FastIO {
    const int SIZE = 1 << 16;
    char buf[SIZE], obuf[SIZE], str[60];
    int bi = SIZE, bn = SIZE, opt;
    double D[] = {0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001, 0.000000001, 0.0000000001};

    int read(char *s) {
        while (bn) {
            for (; bi < bn && buf[bi] <= ' '; bi++);
            if (bi < bn)
                break;
            bn = fread(buf, 1, SIZE, stdin);
            bi = 0;
        }
        int sn = 0;
        while (bn) {
            for (; bi < bn && buf[bi] > ' '; bi++)
                s[sn++] = buf[bi];
            if (bi < bn)
                break;
            bn = fread(buf, 1, SIZE, stdin);
            bi = 0;
        }
        s[sn] = 0;
        return sn;
    }

    bool read(int &x) {
        int n = read(str), bf = 0;
        if (!n)
            return 0;
        int i = 0;
        if (str[i] == '-')
            bf = 1, i++;
        else if (str[i] == '+')
            i++;
        for (x = 0; i < n; i++)
            x = x * 10 + str[i] - '0';
        if (bf)
            x = -x;
        return 1;
    }

    bool read(long long &x) {
        int n = read(str), bf;
        if (!n)
            return 0;
        int i = 0;
        if (str[i] == '-')
            bf = -1, i++;
        else
            bf = 1;
        for (x = 0; i < n; i++)
            x = x * 10 + str[i] - '0';
        if (bf < 0)
            x = -x;
        return 1;
    }

    void write(int x) {
        if (x == 0)
            obuf[opt++] = '0';
        else {
            if (x < 0)
                obuf[opt++] = '-', x = -x;
            int sn = 0;
            while (x)
                str[sn++] = x % 10 + '0', x /= 10;
            for (int i = sn - 1; i >= 0; i--)
                obuf[opt++] = str[i];
        }
        if (opt >= (SIZE >> 1)) {
            fwrite(obuf, 1, opt, stdout);
            opt = 0;
        }
    }

    void write(long long x) {
        if (x == 0)
            obuf[opt++] = '0';
        else {
            if (x < 0)
                obuf[opt++] = '-', x = -x;
            int sn = 0;
            while (x)
                str[sn++] = x % 10 + '0', x /= 10;
            for (int i = sn - 1; i >= 0; i--)
                obuf[opt++] = str[i];
        }
        if (opt >= (SIZE >> 1)) {
            fwrite(obuf, 1, opt, stdout);
            opt = 0;
        }
    }

    void write(unsigned long long x) {
        if (x == 0)
            obuf[opt++] = '0';
        else {
            int sn = 0;
            while (x)
                str[sn++] = x % 10 + '0', x /= 10;
            for (int i = sn - 1; i >= 0; i--)
                obuf[opt++] = str[i];
        }
        if (opt >= (SIZE >> 1)) {
            fwrite(obuf, 1, opt, stdout);
            opt = 0;
        }
    }

    void write(char x) {
        obuf[opt++] = x;
        if (opt >= (SIZE >> 1)) {
            fwrite(obuf, 1, opt, stdout);
            opt = 0;
        }
    }

    void Fflush() {
        if (opt)
            fwrite(obuf, 1, opt, stdout);
        opt = 0;
    }
}; // namespace FastIO
