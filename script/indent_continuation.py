#!/usr/bin/env python3
# Indent continuation lines for fortran files read through stdin
# usage: 
#  cat file.f90 | indent_continuation.py
#
import re
import sys

indent_chars = [ ":", "=", "::" ]
indent_regex = re.compile(r"^\ *")
continuation_regex = re.compile(r'&\ *$')
indent_next = False
current_indent = 0
parenthesis = 0

for line in sys.stdin:
    if not indent_next:
        for char in indent_chars:
            if char in line:
                current_indent = line.find(char) + len(char) + 1

    if indent_next and parenthesis == 0:
        line = indent_regex.sub(current_indent*" ", line)
    indent_next = continuation_regex.search(line)
    print(line, end="")

    parenthesis += line.count("(") - line.count(")")\
                 + line.count("[") - line.count("]")

