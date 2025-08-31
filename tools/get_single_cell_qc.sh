#################################################
#  File Name:get_single_cell_qc.sh
#  Author: Pengwei.Xing
#  Mail: xingwei421@qq.com,pengwei.xing@igp.uu.se,xpw1992@gmail.com
#  Created Time: Fri 01 Aug 2025 03:58:59 PM CEST
#################################################

python get_singcell_qc.py -b1 02.mapping/Spleen.sort.markdup.merge.bam -b2 02.mapping/{samplename}.sort.rmdup.merge.pair.withinpeaks.bam -b3 02.mapping/{samplename}.sort.rmdup.merge.pair.withintss.bam -o {samplename}.singlecell.qc.txt -q 2
