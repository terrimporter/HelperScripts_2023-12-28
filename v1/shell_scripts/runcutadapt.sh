#!/bin/zsh
#Script to grab all paired end files, loop through, run cutadapt multiple times to grab sequences with certain sets of trimmed primers
#USAGE zsh runcutadapt.sh


for r1 in *R1.fastq
do

echo $r1

base=$r1[1,4]
suffix=_R2.fastq

r2=$base$suffix

echo $r2

#CO1 trim B and E primers, 5'-3' do not revcom reverse primer, keeping trimmed seqs only, run seqprep later
#CO1 BE amplicon, B, E
cutadapt -g CCDGAYATRGCDTTYCCDCG -G GTRATDGCDCCDGCDARDAC --trimmed-only -m 200 -M 300 -q 20,20 --max-n=3 -o ${r1}.Btrimmed -p ${r2}.Etrimmed $r1 $r2

#CO1 F230 amplicon, LCO1490, 230_R
cutadapt -g GGTCAACAAATCATAAAGATATTGG -G CTTATRTTRTTTATDCGDGGRAADGC --trimmed-only -m 200 -M 300 -q 20,20 --max-n=3 -o ${r1}.Ltrimmed -p ${r2}.2trimmed $r1 $r2

#16S V3-V4 amplicon, 16S_V3_F, 16S_V4_R
cutadapt -g ACTCCTACGGGAGGCAGCAG -G GGACTACARGGTATCTAAT --trimmed-only -m 200 -M 300 -q 20,20 --max-n=3 -o ${r1}.V3trimmed -p ${r2}.V4trimmed $r1 $r2

#Fish Mini_SH-E amplicon, Fish_miniE_F, Fish_miniE_R
cutadapt -g ACYAADCAYAAAGAYATDGGCAC -G CTTATRTTRTTTATDCGDGGRAADGC --trimmed-only -m 200 -M 300 -q 20,20 --max-n=3 -o ${r1}.EFtrimmed -p ${r2}.ERtrimmed $r1 $r2

#Fish Mini_SH-F amplicon, Fish_miniF_F, Fish_miniF_R
cutadapt -g GGDACDGGDTGRACDGTDTAYCCYCC -G CTTCAGGGTGDCCGAARAATC --trimmed-only -m 200 -M 300 -q 20,20 --max-n=3 -o ${r1}.FFtrimmed -p ${r2}.FRtrimmed $r1 $r2

done
	

