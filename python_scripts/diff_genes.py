# open file1 = has fewer genes
# open file2 = for genes not found in file1, print line

import sys

g1 = []
for line in open(sys.argv[1]):
  col = line.strip().split(',')
  g1.append(col[0])

for line in open(sys.argv[2]):
  col = line.strip().split(',')
  if col[0] not in g1:
    sys.stdout.write(line)