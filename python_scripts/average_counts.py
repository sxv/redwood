import sys

all_counts = []
for l,line in enumerate(open(sys.argv[1])):
  if l > 0:
    counts = line.strip().split(',')[1:]
    if counts[0] != 'NA':
      all_counts += counts

total_counts = 0
for c in all_counts:
  total_counts += float(c)

print('%s/%s = %s' % (total_counts,len(all_counts), total_counts/len(all_counts)))