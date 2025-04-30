#!/usr/bin/env python3

import sys
import os
import time
import subprocess
from glob import glob

# Constants
BASE_DIR = '/media/data/S3S'
MAKE_FASTQ = os.path.join(BASE_DIR, '3SEQtools/make_fastq.sh')
ALIGN = os.path.join(BASE_DIR, '3SEQtools/align_smart-3seq_v2.sh')
MAKE_EXPRESSION_TABLE = os.path.join(BASE_DIR, '3SEQtools/make_expression_table.R')
GENOME_DIR = '/media/data/sxv/star7_gencode33-68/'
GTF_FILE = '/media/data/sxv/gencode.v33.annotation.gtf'

def run(cmd, cwd=None):
    subprocess.run(cmd, shell=True, cwd=cwd, check=True)

def make_fastq(batch_dir):
    fastq_dir = os.path.join(batch_dir, 'fastq')
    if not glob(os.path.join(fastq_dir, '*.gz')):
        run(f'cp {MAKE_FASTQ} {batch_dir}')
        run(f'bash {os.path.basename(MAKE_FASTQ)} .', cwd=batch_dir)
        os.makedirs(fastq_dir, exist_ok=True)
        run('mv *.gz fastq', cwd=batch_dir)
        print('make_fastq complete.')
    else:
        print('fastq exists.')

def align(batch_dir):
    fastq_dir = os.path.join(batch_dir, 'fastq')
    bam_dir = os.path.join(batch_dir, 'bam')
    if glob(os.path.join(fastq_dir, '*.gz')) and not glob(os.path.join(bam_dir, '*.bam')):
        run(f'cp {ALIGN} {fastq_dir}')
        run(f'bash {os.path.basename(ALIGN)} -d {GENOME_DIR} *fastq.gz', cwd=fastq_dir)
        os.makedirs(bam_dir, exist_ok=True)
        run('mv *.bam *.log ../bam', cwd=fastq_dir)
        print('align complete.')
    else:
        print('bam exists.')

def qc(batch_dir):
    bam_dir = os.path.join(batch_dir, 'bam')
    if glob(os.path.join(bam_dir, '*.bam')) and not glob(os.path.join(bam_dir, '*multiqc*')):
        run('for f in *align.log; do cp -L "$f" "${f//align.log/Log.final.out}"; done && multiqc *out', cwd=bam_dir)
        print('multiqc complete.')
    else:
        print('multiqc already exists.')

def make_expression_table(batch_dir, library_name):
    bam_dir = os.path.join(batch_dir, 'bam')
    if glob(os.path.join(bam_dir, '*.bam')) and not os.path.exists(os.path.join(bam_dir, 'gene_expression.xlsx')):
        run(f'cp {MAKE_EXPRESSION_TABLE} {bam_dir}')
        run(f'Rscript {os.path.basename(MAKE_EXPRESSION_TABLE)} --no-rlog {GTF_FILE} *.bam', cwd=bam_dir)
        run('xlsx2csv -s 1 gene_expression.xlsx raw_reads.csv', cwd=bam_dir)
        run(f'cp raw_reads.csv {library_name}.raw_reads.csv', cwd=bam_dir)
        print('make expression table complete.')
    else:
        print('read table exists.')

def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py /<batchname>/<batchid>/")
        sys.exit(1)

    batch_dir = sys.argv[1].rstrip('/')
    library_name = os.path.basename(batch_dir) or os.path.basename(os.path.dirname(batch_dir))

    print(f'Processing batch: {batch_dir}')
    start_time = time.time()

    make_fastq(batch_dir)
    align(batch_dir)
    qc(batch_dir)
    make_expression_table(batch_dir, library_name)

    print(f'Completed in {round(time.time() - start_time, 2)} seconds.')

if __name__ == '__main__':
    main()
