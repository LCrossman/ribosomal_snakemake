#!/usr/bin/python


import sys
import re

def check_excluded_strains(inlist, finlist):
    results = []
    for fin in finlist:
        pattern = re.compile(fin)
        for match in [inl for inl in inlist if pattern.match(inl)]:
            print(match)
            results.append(match)
    return results


infile = open("cleannames.txt", 'r')
excl = open("finalnames.txt", 'r')
ids = []

inlist = [lin.rstrip() for lin in infile]
finlist = [line.rstrip() for line in excl]


results = check_excluded_strains(inlist, finlist)

res = set(inlist) - set(results)

outfile = open("excluded_strains.txt", 'w')

for exc in res:
    outfile.write("{}\n".format(exc))
        
