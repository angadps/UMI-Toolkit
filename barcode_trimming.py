"""
fastq1=$1/*R1*gz
fastq2=$1/*R2*gz
prefix=$2
odir=$3

zcat $fastq1 | paste - - - - | tr '\t' ' ' |
awk -F" " '{bc=substr($3,1,12);qt=substr($5,1,12); n=match(bc,"TTAT");
if(n>0) {bc=substr(bc,1,n+3);qt=substr(qt,1,n+3);} printf "%s BC:Z:%s\n%s\n%s\n%s\n",$1,bc,substr($3,13),$4,substr($5,13);}' > ${odir}/${prefix}-1.fq.bk



zcat $fastq2 | paste - - - - | tr '\t' ' ' | awk -F" " '{bc=substr($3,1,12);qt=substr($5,1,12); n=match(bc,"TTAT"); if(n>0) {bc=substr(bc,1,n+3);qt=substr(qt,1,n+3);} printf "%s BC:Z:%s\n%s\n%s\n%s\n",$1,bc,substr($3,13),$4,substr($5,13);}' > ${odir}/${prefix}-2.fq.bk
#cat ${odir}/${prefix}-1.fq.bk ${odir}/${prefix}-2.fq.bk | grep ^@NRUSCA | cut -d' ' -f 2 | sort | uniq -c > ${odir}/${prefix}.count

rm ${odir}/${prefix}-1.fq ${odir}/${prefix}-2.fq
paste ${odir}/${prefix}-1.fq.bk ${odir}/${prefix}-2.fq.bk | paste - - - - | tr '\t' ' ' | awk -F" " -v fastq1=${odir}/${prefix}-1.fq -v fastq2=${odir}/${prefix}-2.fq '{bc=$2"-"substr($4,6); printf "%s %s\n%s\n%s\n%s\n", $1,bc,$5,$7,$9 >> fastq1; printf "%s %s\n%s\n%s\n%s\n", $3,bc,$6,$8,$10 >> fastq2;}'

"""
from sys import argv
import gzip

script, fastq_file1, fastq_file2, out1, out2 = argv


def get_barcode(seq):
    barcode_index = seq.rfind("TTAT",0, 12)
    if barcode_index != -1:
        barcode = seq[0:barcode_index+4]
        new_sequence = seq[12::]
    else:
        barcode = seq[0:12]
        new_sequence = seq[12::]
    return (barcode, new_sequence)

output1 = gzip.open(out1, 'wt')
output2 = gzip.open(out2, 'wt')
line_index = 1
fastq1 = []
fastq2 = []
with gzip.open(fastq_file1,'rt') as fq1, gzip.open(fastq_file2,'rt') as fq2:
    for line1, line2 in zip(fq1, fq2):
        line1 = line1.strip()
        line2 = line2.strip()
        if line_index == 1:
            ##This is the seqid
            fastq1 = []
            fastq2 = []
            seq_id1 = line1.split()[0]
            seq_id2 = line2.split()[0]
        elif line_index == 2:
            #This is the sequence. Get the barcode and adjust the sequence. Add barcode to the seq_id
            #barcode1_result = get_barcode(line1)
            #barcode2_result = get_barcode(line2)
            #if barcode1_result and barcode2_result:
            barcode1, seq1 = get_barcode(line1)
            barcode2, seq2 = get_barcode(line2)
            seq_id1 = seq_id1 + " BC:Z:" + barcode1 + "-" + barcode2
            seq_id2 = seq_id2 + " BC:Z:" + barcode1 + "-" + barcode2
            fastq1.append(seq_id1)
            fastq1.append(seq1)
            fastq2.append(seq_id2)
            fastq2.append(seq2)
        elif line_index == 3:
            #This is '+'
            fastq1.append(line1)
            fastq2.append(line2)
        elif line_index == 4:
            #This is the quality score. Adjust based on barcode length and print out the 4 lines.
            fastq1.append(line1[12::])
            fastq2.append(line2[12::])
            output1.write("\n".join(fastq1))
            output1.write("\n")
            output2.write("\n".join(fastq2))
            output2.write("\n")
            line_index = 0
        line_index += 1

