#!usr/bin/env pyhton
# coding:utf-8


def GPS_cycleslip():
    for i in range(1, 100):
        for j in range(1, 100):
            for k in range(1, 100):
                condition_1 = abs(-6 * i + 1 * j + 7 * k) <= 1
                condition_2 = abs(3 * i + 0 * j - 4 * k) <= 1
                condition_3 = abs(4 * i + -8 * j + 3 * k) <= 1
                if condition_1 and condition_2 and condition_3:
                    print(i, j, k)


def BDS_cycleslip():
    for i in range(1, 100):
        for j in range(1, 100):
            for k in range(1, 100):
                condition_1 = abs(-4 * i + 1 * j + 4 * k) <= 1
                condition_2 = abs(-3 * i + 6 * j - 2 * k) <= 1
                condition_3 = abs(4 * i - 2 * j - 3 * k) <= 1
                if condition_1 and condition_2 and condition_3:
                    print(i, j, k)


print('BDS:')
BDS_cycleslip()
print('GPS')
GPS_cycleslip()
