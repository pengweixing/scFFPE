#################################################
#  File Name:debacode.sh
#  Author: Pengwei.Xing
#  Mail: xingwei421@qq.com,pengwei.xing@igp.uu.se,xpw1992@gmail.com
#  Created Time: Thu 31 Jul 2025 03:12:09 PM CEST
#################################################

python fastq_debarcoding.allow_mismatch.v3.1.py -r1 1m.R1_001.fastq.gz -r2 1m.R2_001.fastq.gz -b barcode.txt -o 01.debarcoding -p 20 -ml 40
