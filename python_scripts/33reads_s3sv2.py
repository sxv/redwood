# usage:
#   python 3reads.py /batch/dir1/ /batch/dir2/        # run all steps on two batches
#   python 3reads.py -2 -3 -keep /batch/dir    # run steps 2 and 3 and use STAR --loadAndKeep
#   python 3reads.py /batch/dir1 /batch/dir2/ -13     # run steps 1 and 3 on two batches
import sys, os, time
from glob import glob

base_dir = '/media/data/sxv/S3S/'
MAKE_FASTQ = '%s/3SEQtools/make_fastq.sh' % base_dir
ALIGN = '%s/3SEQtools/align_smart-3seq_v2.sh' % base_dir
MAKE_EXPRESSION_TABLE = '%s/3SEQtools/make_expression_table.R' % base_dir

batch_dir = sys.argv[1]
sys.stdout.write('processing batch: %s. ' % batch_dir)
library_name = batch_dir.split('/')[-1] or batch_dir.split('/')[-2]
save_dir = '/media/stroma-s3s.results/%s/' % batch_dir.replace('data','')
os.system("bash -c 'mkdir -p %s'" % save_dir)

def make_fastq():
  if not glob('%s/fastq/*gz' % batch_dir):
    os.system("bash -c 'cp %s %s'" % (MAKE_FASTQ, batch_dir))
    os.system("bash -c 'cd %s && bash %s .'" % (batch_dir, MAKE_FASTQ))
    os.system("bash -c 'cd %s && mkdir -p fastq && mv *gz fastq'" % batch_dir)
    os.system("bash -c 'cp -nR %s/fastq %s &'" % (batch_dir, save_dir))
    sys.stdout.write('make_fastq complete.')
  else: sys.stdout.write('fastq exists.')

def align():
  if glob('%s/fastq/*gz' % batch_dir) and not glob('%s/bam/*bam' % batch_dir):
    os.system("bash -c 'cp %s %s/fastq/'" % (ALIGN, batch_dir))
    os.system("bash -c 'cd %s/fastq && bash %s -d /media/data/sxv/star7_gencode33-68/ *fastq.gz'" % (batch_dir, ALIGN))
    os.system("bash -c 'cd %s && mkdir bam && mv fastq/*ba? fastq/*log bam/'" % batch_dir)
    sys.stdout.write('align complete.')
  else: sys.stdout.write('bam exists.')

def qc():
  if glob('%s/bam/*bam' % batch_dir) and not glob('%s/bam/*multiqc*' % batch_dir):
    os.system("bash -c 'cd %s/bam && for f in *align.log ; do cp -L $f ${f//align.log/Log.final.out} ; done && multiqc *out'" % batch_dir)
  else:
    print 'multiqc already exists.'

def make_expression_table():
  if glob('%s/bam/*bam' % batch_dir) and not glob('%s/bam/gene_expression.xlsx' % batch_dir):
    os.system("bash -c 'cp %s %s/bam'" % (MAKE_EXPRESSION_TABLE, batch_dir))
    os.system("bash -c 'cd %s/bam && Rscript %s --no-rlog /media/data/sxv/gencode.v33.annotation.gtf *bam'" % (batch_dir, MAKE_EXPRESSION_TABLE))
    os.system("bash -c 'cd %s/bam && xlsx2csv -s 1 gene_expression.xlsx raw_reads.csv'" % batch_dir)
    os.system("bash -c 'cd %s/bam && cp raw_reads.csv %s.raw_reads.csv'" % (batch_dir, library_name))
  else:
    sys.stdout.write('make expression table complete.')
  sys.stdout.write('read table exists.')

def save_results():
  os.system("bash -c 'mkdir -p %s'" % save_dir)
  os.system("bash -c 'cp -nR %s/fastq %s'" % (batch_dir, save_dir))
  os.system("bash -c 'cp -nR %s/bam %s'" % (batch_dir, save_dir))
  os.system("bash -c 'cp %s/bam/*.raw_reads.csv %s/bam/*report*html %s'" % (batch_dir, batch_dir, save_dir))
  sys.stdout.write('results saved to %s. ' % save_dir)

starts = time.time()
make_fastq()
align()
qc()
make_expression_table()
save_results()
sys.stdout.write('completed in %s seconds\n.' % (time.time() - starts))
