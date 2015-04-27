#!/bin/bash
#$ -S /bin/bash -cwd
#$ -o TAXON.out -j y
#$ -M k8hertweck@gmail.com -m be
#$ -l highprio
#$ -N TAXON

#classify DNA TEs by superfamily

#go to directory where input file of taxon name is located
cd /Users/kate/Desktop/REdata/Asparagales/$1/pipeline

#DNA
cd DNA
echo 'TOTAL DNA TRANSPOSONS' | tee ../$1superfam.out
cat DNA.out | awk '{print $5}' | sort | uniq | wc -l | tee -a ../$1superfam.out

echo 'EnSpm' | tee -a ../$1superfam.out
	cat DNA.out | awk '$11 ~ /EnSpm/ {print $0}' DNA.out > DNAEnSpm.out
	cat DNAEnSpm.out | awk '{print $5}' | uniq | sort > DNAEnSpm.lst		
	wc -l DNAEnSpm.lst | tee -a ../$1superfam.out
	grep -f DNAEnSpm.lst ../scaf/scafreads.lst | awk '{s+=$3}END{print s}' | tee -a ../$1superfam.out
	
echo 'hAT' | tee -a ../$1superfam.out
	cat DNA.out | awk '$11 ~ /hAT/ {print $0}' DNA.out > DNAhAT.out
	cat DNAhAT.out | awk '{print $5}' | uniq | sort > DNAhAT.lst		
	wc -l DNAhAT.lst | tee -a ../$1superfam.out
	grep -f DNAhAT.lst ../scaf/scafreads.lst | awk '{s+=$3}END{print s}' | tee ../$1superfam.out

echo 'MuDR' | tee -a ../$1superfam.out
	cat DNA.out | awk '$11 ~ /MuDR/ {print $0}' DNA.out > DNAMuDR.out
	cat DNAMuDR.out | awk '{print $5}' | uniq | sort > DNAMuDR.lst		
	wc -l DNAMuDR.lst | tee -a ../$1superfam.out
	grep -f DNAMuDR.lst ../scaf/scafreads.lst | awk '{s+=$3}END{print s}' | tee ../$1superfam.out
	
echo 'PIF' | tee -a ../$1superfam.out
	cat DNA.out | awk '$11 ~ /PIF/ {print $0}' DNA.out > DNAPIF.out
	cat DNAPIF.out | awk '{print $5}' | uniq | sort > DNAPIF.lst		
	wc -l DNAPIF.lst | tee -a ../$1superfam.out
	grep -f DNAPIF.lst ../scaf/scafreads.lst | awk '{s+=$3}END{print s}' | tee ../$1superfam.out

echo 'TcMar' | tee -a ../$1superfam.out
	cat DNA.out | awk '$11 ~ /TcMar/ {print $0}' DNA.out > DNATcMar.out
	cat DNATcMar.out | awk '{print $5}' | uniq | sort > DNATcMar.lst		
	wc -l DNATcMar.lst | tee -a ../$1superfam.out	
	grep -f DNATcMar.lst ../scaf/scafreads.lst | awk '{s+=$3}END{print s}' | tee ../$1superfam.out
	
