## run as $ python flatten.py '*.txt'

import sys
from glob import glob

if len(sys.argv) > 2: filter = 0
else: filter = 1

pattern = sys.argv[1]

for f, file in enumerate( glob(pattern) ):
  for l, line in enumerate( open(file) ):
    if f == 0 and l == 0:
      sys.stdout.write("Sample\t%s" % line)
    elif l != 0:
      col = line.split("\t")
      if not filter or (col[5] in ['exonic', 'splicing'] and not col[8].startswith('syn') and not col[13].startswith('syn')):
        sys.stdout.write('%s\t%s' % (file, line))
