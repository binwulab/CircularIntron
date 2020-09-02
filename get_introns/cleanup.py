"""Cleans up bed file with scientific notation."""
from sys import argv

filename = argv[1]

with open(filename, 'r') as f:
    for line in f:
        items = line.strip().split()
        start = str(int(float(items[1])))
        stop = str(int(float(items[2])))
        new_items = [items[0]] + [start, stop] + items[3:]
        print('\t'.join(new_items))
