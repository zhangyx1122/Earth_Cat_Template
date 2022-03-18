int gray_encode(int num) {
    return num ^ (num >> 1);
}

int gray_decode(int num) {
    int head;
    if (!num) return 0;
    head = 1 << int(log(num) / log(2));
    return head + gray_decode((num ^ head) ^ (head >> 1));
}

