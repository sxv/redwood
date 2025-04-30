# for sample in matrix, how many genes have at least <threshold X> reads?
# desired output:
# dd = {"sample1": [2342,289,14,2],..} showing sample1 has 2342 genes with at least 1 read, 289 with at least 10 reads..


epi100 = ['MYO6','SPINT1','KRT18','PTPRF','XBP1','TFAP2A','CLDN7','ARFGEF3','SPINT2','ERBB3','SLC38A1','FOXA1','PPDPF','DDR1','MARCKSL1','DBNDD1','TRPS1','PRSS8','MLPH','SPDEF','BICDL2','EPCAM','MALAT1','DHCR24','MUC1','CREB3L4','GALNT7','ELF3','EZR','CASZ1','KIAA1522','CLDN4','GATA3','TACSTD2','SREBF1','AZGP1','KRT19','ESRP1','DSP','ELAPOR1','PRDX2','JUP','GRHL2','CRNDE','MAL2','IGSF9','SNHG25','CLDN3','TSTD1','SMIM14','HSP90AB1','PATJ','CGN','ARSD','KCTD3','NR2F6','HDGF','GRHL1','CD24','TMBIM6','ALCAM','ATP1B1','AR','PROM2','METRN','FAM83H','TOB1','H2AJ','TPD52','CDS1','JPT2','CA12','MSI2','SOX4','IRX3','CYB561','RGL2','TSPAN13','SDC4','MAP7','SHROOM3','SLC9A3R1','IRX5','PRLR','MYB','CD9','COX6C','SERINC2','STARD10','MAGED2','ZBTB7B','CDH1','RBM47','SLC39A6','ZNF587','PKP3','TMEM184A','NECTIN4','GRB7']

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
  else:
    gene = col[0]
    gene_name = col[1]
    if gene_name in epi100:
      # print(gene)
      counts = col[2:]
      for ci, c in enumerate(counts):
        s = samples[ci]
        for i,t in enumerate(thresholds):
          # print(t)
          if float(c) >= t:
            dd[s][i].append(gene)


for s in samples:
  output = [s]
  for i,t in enumerate(thresholds):
    output.append(str(len(dd[s][i])))
  print(','.join(output))
