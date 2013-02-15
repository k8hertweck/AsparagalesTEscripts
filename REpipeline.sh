#!/bin/bash
#$ -S /bin/bash -cwd
#$ -o EL01TE.out -j y
#$ -M k8hertweck@gmail.com -m e
#$ -l highprio
#$ -N EL01TE

##REPLACE TAXON WITH NAME
##TAXON FOLDER IN /home/nescent/kh200/repeats
##READS IN /netscratch/kh200/data
##REF INDEX FILES IN /home/nescent/kh200/data/refs
##THREE FILES: .sh, sr
##sr_config.txt: READ FILE, PATH, SET READ PARAMETERS
##THINGS TO CHANGE IN THIS FILE?
##MUST HAVE PROGRAMS INSTALLED: MSR, repeatmasker, smalt, samtools, seqtk, cdbyank
#need to add: sending echoed output to file

cd /Users/kate/Desktop/REdata/Poaceae/EL01/pipeline

##MSR ON ALL RAW READS (sr_config.txt is in taxon folder)
	#echo 'MSR' 
	#mkdir temp MSR 
	#cd temp
	#perl /opt/apps/MSR-CA/bin/runSRCA.pl ../sr_config.txt
	#./assemble.sh

	#cd .. 
	#mv sr_config.txt temp/CA/9-terminator/* MSR/
	#rm -r temp
	mkdir scaf
	cp MSR/genome.scf.fasta scaf/scaf.fas
	
#QUANTIFICATION
	
	echo 'QUANTIFICATION AND PREP'
	
	#list scaf.fas names
	cd scaf
	grep ">" scaf.fas | sed 's/>//' > scaf.lst
	wc -l scaf.lst
	#make index
	cdbfasta scaf.fas
	
	#MAP READS TO SCAFFOLDS
	smalt index scaf scaf.fas 
	smalt map -f sam -o scaf.sam scaf ~/data/EL01TRIM.fastq 
	
	##convert from SAM to BAM
	samtools view -bS -o scaf.bam scaf.sam 
	##sort BAM
	samtools sort scaf.bam scafsort
	##summary stats: how many sequences for each reference
	samtools index scafsort.bam
	samtools idxstats scafsort.bam > scafreads.lst
	
	#SUM COLUMNS AND REPORT STATS	
	echo 'TOTAL READS MAPPED' 
	awk '{s+=$3}END{print s}' scafreads.lst 
	echo 'READS UNMAPPED'
	awk '{s+=$4}END{print s}' scafreads.lst 

##REMOVE CP/MT FROM SCAFFOLDS (START IN TAXON FOLDER) and reads

	echo 'REMOVE CP/MT SCAFFOLDS'
	
	##blast and parse mtcp scaffolds
	blastall -p blastn -W 50 -m 9 -d /Users/kate/Desktop/REdata/referenceSeq/organelles/ZeaOrg.fasta -i scaf.fas -o cpmtblast.out
	awk '$0 ~ /MT/ {print $1}' cpmtblast.out | sort | uniq > mtscafs.lst
	wc -l mtscafs.lst
	awk '$0 ~ /CP/ {print $1}' cpmtblast.out | sort | uniq > cpscafs.lst
	wc -l cpscafs.lst
	
	#pull out mtcp scaffolds
	cat mtscafs.lst cpscafs.lst | sort | uniq > scafCPMT.lst
	wc -l scafCPMT.lst
	#quantify reads for cpmt
	echo 'READS MAPPED TO CPMT'
	grep -f scafCPMT.lst scafreads.lst | awk '{s+=$3}END{print s}'

#RUN REPEATMASKER ON ALL SCAFFOLDS (start in taxon folder)
	echo 'REPEATMASKER' 
	repeatmasker -gff -species liliopsida scaf.fas 

#UNIQUE REPEAT SCAFFOLDS
	#FIND UNIQUE REPEAT SCAFFOLDS
	echo 'UNIQUE REPEAT SCAFFOLDS'
	tail -n+4  scaf.fas.out | awk '{print $5}' | sort | uniq > scafRE.lst
	wc -l scafRE.lst 
	
	#TOTAL REPEAT READS MAPPED
	echo 'TOTAL REPEAT READS MAPPED'
	grep -f scafRE.lst scafreads.lst| awk '{s+=$3}END{print s}'
	
#UNKNOWN SCAFFOLDS
	#LIST UNKNOWN SCAFFOLDS
	echo 'UNKNOWN SCAFFOLDS'
	diff scaf.lst scafCPMT.lst > nuc.temp 
	diff nuc.temp scafRE.lst | awk '$2 ~ /scf/ {print $2}' > scafUnknown.lst
	wc -l scafUnknown.lst
	rm nuc.temp
	
	#fasta unknown scafs to BLAST later
	cat scafUnknown.lst | cdbyank scaf.fas.cidx -o scafUnknown.fas
	
	cd ..

#PARSING INTO REPEAT CLASSES (start in taxon folder)
echo 'PARSING REPEATS' 
mkdir retroelements DNA retroelements/LINE retroelements/SINE retroelements/LTRs other
	
	#RETROELEMENTS
	echo 'RETROELEMENTS' 
	echo 'LTRs' 
		#find and count scaffolds annotated as category
		awk '$11 ~ /LTR/ {print $0}' scaf/scaf.fas.out > retroelements/LTRs/LTR.out
		awk '{print $5}' retroelements/LTRs/LTR.out | sort | uniq > retroelements/LTRs/LTRscafs.lst
		wc -l retroelements/LTRs/LTRscafs.lst 
		#pull out scaffolds annotated as category
		cat retroelements/LTRs/LTRscafs.lst | cdbyank scaf/scaf.fas.cidx -o retroelements/LTRs/LTRscafs.fas 
		#count reads mapping to category scaffolds (from original smalt)
		echo 'TOTAL READS MAPPED'
		grep -f retroelements/LTRs/LTRscafs.lst scaf/scafreads.lst | awk '{s+=$3}END{print s}'		
	echo 'Copia' 
		#find and count scaffolds annotated as category
		awk '$11 ~ /Copia/ {print $0}' scaf/scaf.fas.out > retroelements/LTRs/Copia.out
		awk '{print $5}' retroelements/LTRs/Copia.out | sort | uniq > retroelements/LTRs/Copiascafs.lst
		wc -l retroelements/LTRs/Copiascafs.lst 
		#pull out scaffolds annotated as category
		cat retroelements/LTRs/Copiascafs.lst | cdbyank scaf/scaf.fas.cidx -o retroelements/LTRs/Copiascafs.fas 
		#count reads mapping to category scaffolds (from original smalt)
		echo 'TOTAL READS MAPPED'
		grep -f retroelements/LTRs/Copiascafs.lst scaf/scafreads.lst | awk '{s+=$3}END{print s}'		 
	echo 'Gypsy' 
		#find and count scaffolds annotated as category
		awk '$11 ~ /Gypsy/ {print $0}' scaf/scaf.fas.out > retroelements/LTRs/Gypsy.out
		awk '{print $5}' retroelements/LTRs/Gypsy.out | sort | uniq > retroelements/LTRs/Gypsyscafs.lst
		wc -l retroelements/LTRs/Gypsyscafs.lst 
		#pull out scaffolds annotated as category
		cat retroelements/LTRs/Gypsyscafs.lst | cdbyank scaf/scaf.fas.cidx -o retroelements/LTRs/Gypsyscafs.fas 
		#count reads mapping to category scaffolds (from original smalt)
		echo 'TOTAL READS MAPPED'
		grep -f retroelements/LTRs/Gypsyscafs.lst scaf/scafreads.lst | awk '{s+=$3}END{print s}'	
	echo 'SINEs'
		#find and count scaffolds annotated as category
		awk '$11 ~ /SINE/ {print $0}' scaf/scaf.fas.out > retroelements/SINE.out
		awk '{print $5}' retroelements/SINE.out | sort | uniq > retroelements/SINEscafs.lst
		wc -l retroelements/SINEscafs.lst 
		#pull out scaffolds annotated as category
		cat retroelements/SINEscafs.lst | cdbyank scaf/scaf.fas.cidx -o retroelements/SINEscafs.fas 
		#count reads mapping to category scaffolds (from original smalt)
		echo 'TOTAL READS MAPPED'
		grep -f retroelements/SINEscafs.lst scaf/scafreads.lst| awk '{s+=$3}END{print s}'	 
	echo 'LINEs' 
		#find and count scaffolds annotated as category
		awk '$11 ~ /LINE/ {print $0}' scaf/scaf.fas.out > retroelements/LINE.out
		awk '{print $5}' retroelements/LINE.out | sort | uniq > retroelements/LINEscafs.lst
		wc -l retroelements/LINEscafs.lst 
		#pull out scaffolds annotated as category
		cat retroelements/LINEscafs.lst | cdbyank scaf/scaf.fas.cidx -o retroelements/LINEscafs.fas 
		#count reads mapping to category scaffolds (from original smalt)
		echo 'TOTAL READS MAPPED'
		grep -f retroelements/LINEscafs.lst scaf/scafreads.lst| awk '{s+=$3}END{print s}'
	echo 'DNA'
		#find and count scaffolds annotated as category
		awk '$11 ~ /DNA/ {print $0}' scaf/scaf.fas.out > DNA/DNA.out
		awk '{print $5}' DNA/DNA.out | sort | uniq > DNA/DNAscafs.lst
		wc -l DNA/DNAscafs.lst 
		#pull out scaffolds annotated as category
		cat DNA/DNAscafs.lst | cdbyank scaf/scaf.fas.cidx -o DNA/DNAscafs.fas 
		#count reads mapping to category scaffolds (from original smalt)
		echo 'TOTAL READS MAPPED'
		grep -f DNA/DNAscafs.lst scaf/scafreads.lst| awk '{s+=$3}END{print s}'
	echo 'RC'
		#find and count scaffolds annotated as category
		awk '$11 ~ /RC/ {print $0}' scaf/scaf.fas.out > other/RC.out
		awk '{print $5}' other/RC.out | sort | uniq > other/RCscafs.lst
		wc -l other/RCscafs.lst 
		#pull out scaffolds annotated as category
		cat other/RCscafs.lst | cdbyank scaf/scaf.fas.cidx -o other/RCscafs.fas 
		#count reads mapping to category scaffolds (from original smalt)
		echo 'TOTAL READS MAPPED'
		grep -f other/RCscafs.lst scaf/scafreads.lst| awk '{s+=$3}END{print s}'
	echo 'Satellite'
		#find and count scaffolds annotated as category
		awk '$11 ~ /Satellite/ {print $0}' scaf/scaf.fas.out > other/Satellite.out
		awk '{print $5}' other/Satellite.out | sort | uniq > other/Satellitescafs.lst
		wc -l other/Satellitescafs.lst 
		#pull out scaffolds annotated as category
		cat other/Satellitescafs.lst | cdbyank scaf/scaf.fas.cidx -o other/Satellitescafs.fas 
		#count reads mapping to category scaffolds (from original smalt)
		echo 'TOTAL READS MAPPED'
		grep -f other/Satellitescafs.lst scaf/scafreads.lst| awk '{s+=$3}END{print s}'	
	echo 'rRNA'
		#find and count scaffolds annotated as category
		awk '$11 ~ /rRNA/ {print $0}' scaf/scaf.fas.out > other/rRNA.out
		awk '{print $5}' other/rRNA.out | sort | uniq > other/rRNAscafs.lst
		wc -l other/rRNAscafs.lst 
		#pull out scaffolds annotated as category
		cat other/rRNAscafs.lst | cdbyank scaf/scaf.fas.cidx -o other/rRNAscafs.fas 
		#count reads mapping to category scaffolds (from original smalt)
		echo 'TOTAL READS MAPPED'
		grep -f other/rRNAscafs.lst scaf/scafreads.lst| awk '{s+=$3}END{print s}'	