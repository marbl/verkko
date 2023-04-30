#!/usr/bin/env python3

import sys
import re
import os

spltd =[['no_read'],['no_read']]
if not sys.stdin.isatty():
    input_stream = sys.stdin

# otherwise, read the given filename                                            
else:
    try:
        input_filename = sys.argv[1]
    except IndexError:
        message = 'need filename as first argument if stdin is not full'
        raise IndexError(message)
    else:
        input_stream = open(input_filename, 'rU')

for line in input_stream:
    spltd[0] = spltd[1]
    spltd[1] = line.split()
    if spltd[0][0] == spltd[1][0]:
        if spltd[0][2] != spltd[1][2]:
            print (f"{spltd[0][0]}\t{spltd[0][2]}\t{spltd[1][2]}")
#            print (f"{spltd[0][0]}\t{spltd[0][2]}\t{spltd[1][2]}\t{spltd[0][4]}\t{spltd[1][4]}")

