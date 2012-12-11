#!/bin/bash
#$ -S /bin/bash -cwd
#$ -o TAXON.out -j y
#$ -M k8hertweck@gmail.com -m be
#$ -l highprio
#$ -N TAXON

##MUST HAVE PROGRAMS INSTALLED:
##RepeatMasker
##THIS FILE IN TAXON FOLDER; RESULTS FILE ALREADY PRESENT: CHANGE TAXON

cd /Users/kate/Desktop/REdata/TAXON/annotate

#PREP
echo 'PREP' | tee -a annotate.out
mkdir RM DNA retroelements other retroelements/LTRs

#RUN REPEATMASKER
echo 'REPEATMASKER' | tee -a annotate.out
cp ../results/MSR1/genome.scf.fasta RM/scafNuc.fas
cd RM
repeatmasker -gff -species liliopsida scafNuc.fas | tee -a annotate.out
cd ..
cp RM/scafNuc.fas.out .

#PARSING INTO REPEAT CLASSES	
echo 'TOTAL REPEATS' 
tail +4 scafNuc.fas.out | wc -l | tee -a annotate.out
	
	#RETROELEMENTS
	echo 'RETROELEMENTS' | tee -a annotate.out
	echo 'LTRs' | tee -a annotate.out
		awk '$11 ~ /LTR/ {print $0}' scafNuc.fas.out > retroelements/LTRs/LTR.out
		#SORT LTR FILE
		tail +4 retroelements/LTRs/LTR.out | sort -k11d -o retroelements/LTRs/LTRsort.out
		wc -l retroelements/LTRs/LTRsort.out | tee -a annotate.out
	echo 'Copia' | tee -a annotate.out
		awk '$11 ~ /Copia/ {print $0}' scafNuc.fas.out > retroelements/LTRs/Copia.out
		wc -l retroelements/LTRs/Copia.out | tee -a annotate.out
	echo 'Gypsy' | tee -a annotate.out
		awk '$11 ~ /Gypsy/ {print $0}' scafNuc.fas.out > retroelements/LTRs/Gypsy.out
		wc -l retroelements/LTRs/Gypsy.out | tee -a annotate.out
	echo 'SINEs' | tee -a annotate.out
		awk '$11 ~ /SINE/ {print $0}' scafNuc.fas.out > retroelements/SINEs.out
		wc -l retroelements/SINEs.out | tee -a annotate.out
	echo 'LINEs' | tee -a annotate.out
		awk '$11 ~ /LINE/ {print $0}' scafNuc.fas.out > retroelements/LINEs.out
		wc -l retroelements/LINEs.out | tee -a annotate.out
	
	#DNA
	echo 'DNA' | tee -a annotate.out
	awk '$11 ~ /DNA/ {print $0}' scafNuc.fas.out > DNA/DNA.out
	#SORT DNA FILE
	tail +4 DNA/DNA.out | sort -k11d -o DNA/DNAsort.out
	wc -l DNA/DNAsort.out | tee -a annotate.out	
	#OTHER
	echo 'OTHER' | tee -a annotate.out
	echo 'ROLLING CIRCLES' | tee -a annotate.out
		awk '$11 ~ /RC/ {print $0}' scafNuc.fas.out > other/RC.out
		wc -l other/RC.out | tee -a annotate.out
	echo 'SIMPLE REPEATS' | tee -a annotate.out
		awk '$11 == "Simple_repeat" {print $0}' scafNuc.fas.out > other/simple.out
		wc -l other/simple.out | tee -a annotate.out
	echo 'SATELLITES' | tee -a annotate.out
		awk '$11 ~ /Satellite/ {print $0}' scafNuc.fas.out > other/satellite.out
		wc -l other/satellite.out | tee -a annotate.out
	echo 'rRNA' | tee -a annotate.out
		awk '$11 ~ /rRNA/ {print $0}' scafNuc.fas.out > other/rRNA.out
		wc -l other/rRNA.out | tee -a annotate.out
	echo 'LOW COMPLEXITY' | tee -a annotate.out
		awk '$11 ~ /Low_complexity/ {print $0}' scafNuc.fas.out > other/lowcompl.out
		wc -l other/lowcompl.out | tee -a annotate.out