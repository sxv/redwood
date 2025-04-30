# for sample in matrix, print total reads per sample

import sys

f = sys.argv[1]

dd = {}
for l,line in enumerate(open(f)):
  col = line.strip().split(',')
  if l == 0:
    samples = col[2:]
    for s in samples:
      dd[s] = 0
  else:
    gene = col[0]
    counts = col[2:]
    for ci, c in enumerate(counts):
      s = samples[ci]
      dd[s] += int(c)

for s in samples:
  print('%s,%s' % (s,dd[s]))
