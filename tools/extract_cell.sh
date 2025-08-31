#################################################
#  File Name:extract_cell.sh
#  Author: Pengwei.Xing
#  Mail: xingwei421@qq.com,pengwei.xing@igp.uu.se,xpw1992@gmail.com
#  Created Time: Fri 01 Aug 2025 03:35:51 PM CEST
#################################################


Rscript plot_qc.r -i 06.QC/{sample}.singlecell.qc.txt -p 15 -t 15 -n 1000 -o {sample}_selected.cell.names}
#params: -p = FRiP, -t = FRiT, -n = The number of frags

python extract_cell.py  -n 06.QC/{sample}_selected.cell.names  -f 02.mapping/{sample}.rmdup.merge.pair.sort.fragments.gz -b 02.mapping/{sample}.rmdup.merge.pair.singlcell.sort.bam -ob {output.selectedbam} -of {output.selectedfrags} && bgzip -c {output.selectedfrags} > {output.selectedfrags2}
