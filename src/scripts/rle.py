#!/usr/bin/env python

import fileinput
import re

def print_lines(current_lines):
        print("".join(current_lines))

current_lines = []

for l in fileinput.input():
        if l[0] == '>':
                if len(current_lines) > 0: print_lines(current_lines)
                print(l.strip())
                current_lines = []
                continue
        li = re.sub("(.)\\1+", "\\1", l.strip())
        if len(li) == 0: continue
        if len(current_lines) == 0:
                current_lines.append(li)
        else:
                if current_lines[-1][-1] == li[0]:
                        if len(li) == 1: continue
                        current_lines.append(li[1:])
                else:
                        current_lines.append(li)
if len(current_lines) > 0: print_lines(current_lines)
