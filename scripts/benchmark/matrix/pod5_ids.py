#!/usr/bin/env python3
"""Print the read ids held in a pod5 file or directory (metadata only, no signal). Usage: pod5_ids.py <pod5 file|dir>"""
import os, sys
import pod5
src = sys.argv[1]
files = [src] if os.path.isfile(src) else sorted(os.path.join(d, f) for d, _, fs in os.walk(src) for f in fs if f.endswith(".pod5"))
n = 0
for f in files:
    with pod5.Reader(f) as r:
        for rid in r.read_ids:
            sys.stdout.write(str(rid) + "\n"); n += 1
print(f"pod5_ids: {n:,} reads in {len(files)} file(s)", file=sys.stderr)
