#!/usr/bin/env python

import sys, csv

def read_csv(dst, path):
    with open(path, 'r') as f:
        r = csv.reader(f, delimiter='\t')
        for row in r:
            if row[0] > row[5]:
                dst.append((row[5], row[0], row[4]))
            else:
                dst.append((row[0], row[5], row[4]))

a = []
read_csv(a, sys.argv[1])

b = []
read_csv(b, sys.argv[2])

print(float(len(set(a) & set(b))) / len(set(a) | set(b)))
