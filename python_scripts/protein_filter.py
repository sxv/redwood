import sys

pc = []
for line in open(sys.argv[1]):
  pc.append(line.strip())

for line in open(sys.argv[2]):
  line = line.strip()
  gene = line.split(',')[0]
  if gene in pc:
    print(line)