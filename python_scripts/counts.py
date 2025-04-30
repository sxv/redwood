# for sample in matrix, how many genes have at least <threshold X> reads?
# desired output:
# dd = {"sample1": [2342,289,14,2],..} showing sample1 has 2342 genes with at least 1 read, 289 with at least 10 reads..

thresholds = [1,10,100,1000]

import sys
vst = sys.argv[1]

dd = {}
for l,line in enumerate(open(vst)):
  col = line.strip().split(',')
  if l == 0:
    samples = col[2:]
    for s in samples:
      dd[s] = [ [], [], [], [] ] # initialize dd[s] to [],[],[],[] to hold gene lists
      # print(dd)
  else:
    gene = col[0]
    # print(gene)
    counts = col[2:]
    for ci, c in enumerate(counts):
      s = samples[ci]
      for i,t in enumerate(thresholds):
        # print(t)
        if float(c) >= t:
          dd[s][i].append(gene)
          # if('15365_722_StrNod2or4Epi' in s): print('appending %s to %s (threshold %s). %s >= %s' % (gene,s,i,c,t))
        # else:
          # print('%s not >= %s' % (float(c), t))

# print(dd['"15365_722_StrNod2or4Epi"'])

for s in samples:
  output = [s]
  for i,t in enumerate(thresholds):
    output.append(str(len(dd[s][i])))
  print(','.join(output))
