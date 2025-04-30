import sys
from glob import glob

if len(sys.argv) > 1:
  pattern = sys.argv[1]
file = glob(pattern)[0]

genes = {}
for l, line in enumerate( open(file) ):
  if l == 0: print line
  gene_name = line.split("\t")[7]
  case = line.split("\t")[0]
  if gene_name not in genes:
    genes[gene_name] = {"lines":[], "unique":[], "recur":[]}
  gene = genes[gene_name]
  gene['lines'].append(line)
  if case in gene['unique']:
    gene['unique'].remove(case)
    gene['recur'].append(case)
  else:
    if case not in gene["recur"]:
      gene["unique"].append(case)
ordered = sorted(genes, key = lambda x: (len(genes[x]["unique"])*-1, len(genes[x]["lines"])))
for gene in ordered:
  print "%s shows exactly one mutation in %s individual samples, and appears in a total of %s calls across samples." % (gene, len(genes[gene]["unique"]), len(genes[gene]["lines"]))
  print "".join(sorted(genes[gene]["lines"], key = lambda x: x.split()[2]))
