# for sample in matrix, how many genes have at least <threshold X> reads?
# desired output:
# dd = {"sample1": [2342,289,14,2],..} showing sample1 has 2342 genes with at least 1 read, 289 with at least 10 reads..


str100 = ['IGFBP7','C1R','COL6A3','COL1A1','CCDC80','AEBP1','COL1A2','VCAN','MMP2','DCN','C1S','SPARC','COL6A2','TCF4','TIMP2','VIM','FSTL1','UACA','GSN','LRP1','COL3A1','LUM','SFRP2','COL6A1','SERPING1','BGN','POSTN','FBN1','COL15A1','HTRA1','ZEB1','CCN2','PRRX1','LAMB1','ZEB2','EBF1','MBNL1','FBLN2','SERPINF1','C3','IGF1','HLA-E','THBS2','MEG3','PCOLCE','IFI16','CAVIN1','NID1','CD93','KCTD12','IGKC','PDGFRB','CRISPLD2','GIMAP7','LHFPL6','CELF2','CD74','AKAP12','CALD1','SFRP4','MS4A6A','MAF','ITGA1','NRP1','CDH11','CXCL12','MEF2C','FBLN1','WIPF1','HLA-DRA','THY1','MXRA5','SAMHD1','EMILIN1','CD34','STAB1','TACC1','HEG1','AQP1','MFAP4','C11orf96','COL5A2','TNS1','MMP14','COL4A1','FBLN5','ISLR','HLA-DPB1','HLA-DPA1','SLIT3','LAMA4','SPARCL1','FYB1','FN1','MYO6','CAV1','CTSK','HLA-DRB1','ELN','ADGRA2']

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
    if gene_name in str100:
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
