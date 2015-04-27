#!/bin/bash
#$ -S /bin/bash -cwd
#$ -o TAXONTE.out -j y
#$ -M k8hertweck@gmail.com -m e
#$ -l highprio
#$ -N TAXONTE

##REPLACE "TAXON" WITH NAME OF TAXON
##TAXON FOLDER IN /home/nescent/kh200/repeats
##READS IN /netscratch/kh200/data
##REF INDEX FILES IN /home/nescent/kh200/data/refs
##sr_config.txt: CHANGE READ FILE, PATH, SET READ PARAMETERS
##DEPENDENCIES: MSR, repeatmasker, smalt, samtools, seqtk, cdbyank

##MSR ON ALL RAW READS (sr_config.txt is in taxon folder)
	echo 'MSR' 
	mkdir temp MSR 
	cd temp
	perl /opt/apps/MSR-CA/bin/runSRCA.pl ../sr_config.txt
	./assemble.sh

	cd .. 
	mv sr_config.txt temp/CA/9-terminator/* MSR/
	rm -r temp
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
	smalt map -f sam -o scaf.sam scaf ~/data/TAXONTRIM.fastq 
	
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
	#REMOVE CPMT SCAFFOLDS
	diff scaf.lst scafCPMT.lst | awk '$2 ~ /scf/ {print $2}' > nuc.lst
	cat nuc.lst | cdbyank scaf.fas.cidx -o nuc.fas

#RUN REPEATMASKER ON ALL SCAFFOLDS (start in taxon folder)
	echo 'REPEATMASKER' 
	repeatmasker -nolow -species liliopsida nuc.fas 

#REPEAT SCAFFOLDS
	
	#REMOVE JUNK HITS FROM REPEATMASKER OUTPUT
	echo 'REMOVE JUNK'
	awk '$11 ~ /RC/ {print $5}' nuc.fas.out > RC.lst
	wc -l RC.lst
	awk '$11 ~ /Low_complexity/ {print $5}' nuc.fas.out > lowcomp.lst
	wc -l lowcomp.lst
	awk '$11 ~ /Simple_repeat/ {print $5}' nuc.fas.out > simple.lst
	wc -l simple.lst
	grep -v "Helitron" nuc.fas.out | grep -v "Low_complexity" | grep -v "Simple_repeat"  > nojunk.fas.out
	
	#TOTAL REPEAT READS MAPPED (AMBIGUOUS + ANNOTATED)
	echo 'TOTAL REPEAT READS MAPPED'
	tail -n+4  nojunk.fas.out | awk '{print $5}' | sort | uniq > scafRE.lst
	wc -l scafRE.lst
	grep -f scafRE.lst scafreads.lst | awk '{s+=$3}END{print s}'	
	
	#FIND SCAFFOLDS WITH AMBIGUOUS ANNOTATIONS
	tail -n+4  nojunk.fas.out | awk '{print $5,$11}' | sort | uniq | awk '{print $1}' | uniq -d > ambtemp.lst
	awk '$11 ~ /Unknown/ {print $5}' nojunk.fas.out > unknown.lst
	wc -l unknown.lst
	cat unknown.lst ambtemp.lst > ambiguous.lst
	wc -l ambiguous.lst
	echo 'AMBIGUOUS READS MAPPED'
	grep -f ambiguous.lst scafreads.lst	| awk '{s+=$3}END{print s}'
	
	#FIND ANNOTATED REPEAT SCAFFOLDS
	echo 'ANNOTATED REPEAT SCAFFOLDS'
	grep -v -f ambiguous.lst nojunk.fas.out > annotate.fas.out
	tail -n+4  annotate.fas.out | awk '{print $5}' | sort | uniq > scafannotate.lst
	wc -l scafannotate.lst 
	#rm RE.temp
	
	#UNKNOWN SCAFFOLDS (INCLUDING JUNK)
	#LIST UNKNOWN SCAFFOLDS
	echo 'UNKNOWN SCAFFOLDS'
	grep -v -f scafRE.lst nuc.lst > scafUnknown.lst
	wc -l scafUnknown.lst
	#fasta unknown scafs to BLAST later
	cat scafUnknown.lst | cdbyank scaf.fas.cidx -o scafUnknown.fas

cd ..	

#PARSING INTO REPEAT CLASSES (start in taxon folder)
echo 'PARSING REPEATS' 
mkdir retroelements DNA retroelements/LINE retroelements/SINE retroelements/LTRs other
	
	#RETROELEMENTS
	echo 'RETROELEMENTS' 
	echo 'LTRs' 
		awk '$11 ~ /LTR/ {print $0}' scaf/annotate.fas.out > retroelements/LTRs/LTR.out
		awk '{print $5}' retroelements/LTRs/LTR.out | sort | uniq > retroelements/LTRs/LTRscafs.lst
		wc -l retroelements/LTRs/LTRscafs.lst 
		cat retroelements/LTRs/LTRscafs.lst | cdbyank scaf/scaf.fas.cidx -o retroelements/LTRs/LTRscafs.fas 
		echo 'TOTAL READS MAPPED'
		grep -f retroelements/LTRs/LTRscafs.lst scaf/scafreads.lst | awk '{s+=$3}END{print s}'		
	echo 'Copia' 
		awk '$11 ~ /Copia/ {print $0}' scaf/annotate.fas.out > retroelements/LTRs/Copia.out
		awk '{print $5}' retroelements/LTRs/Copia.out | sort | uniq > retroelements/LTRs/Copiascafs.lst
		wc -l retroelements/LTRs/Copiascafs.lst 
		cat retroelements/LTRs/Copiascafs.lst | cdbyank scaf/scaf.fas.cidx -o retroelements/LTRs/Copiascafs.fas 
		echo 'TOTAL READS MAPPED'
		grep -f retroelements/LTRs/Copiascafs.lst scaf/scafreads.lst | awk '{s+=$3}END{print s}'		 
	echo 'Gypsy' 
		awk '$11 ~ /Gypsy/ {print $0}' scaf/annotate.fas.out > retroelements/LTRs/Gypsy.out
		awk '{print $5}' retroelements/LTRs/Gypsy.out | sort | uniq > retroelements/LTRs/Gypsyscafs.lst
		wc -l retroelements/LTRs/Gypsyscafs.lst 
		cat retroelements/LTRs/Gypsyscafs.lst | cdbyank scaf/scaf.fas.cidx -o retroelements/LTRs/Gypsyscafs.fas 
		echo 'TOTAL READS MAPPED'
		grep -f retroelements/LTRs/Gypsyscafs.lst scaf/scafreads.lst | awk '{s+=$3}END{print s}'	
	echo 'SINEs'
		awk '$11 ~ /SINE/ {print $0}' scaf/annotate.fas.out > retroelements/SINE.out
		awk '{print $5}' retroelements/SINE.out | sort | uniq > retroelements/SINEscafs.lst
		wc -l retroelements/SINEscafs.lst 
		cat retroelements/SINEscafs.lst | cdbyank scaf/scaf.fas.cidx -o retroelements/SINEscafs.fas 
		echo 'TOTAL READS MAPPED'
		grep -f retroelements/SINEscafs.lst scaf/scafreads.lst| awk '{s+=$3}END{print s}'	 
	echo 'LINEs' 
		awk '$11 ~ /LINE/ {print $0}' scaf/annotate.fas.out > retroelements/LINE.out
		awk '{print $5}' retroelements/LINE.out | sort | uniq > retroelements/LINEscafs.lst
		wc -l retroelements/LINEscafs.lst 
		cat retroelements/LINEscafs.lst | cdbyank scaf/scaf.fas.cidx -o retroelements/LINEscafs.fas 
		echo 'TOTAL READS MAPPED'
		grep -f retroelements/LINEscafs.lst scaf/scafreads.lst| awk '{s+=$3}END{print s}'
	echo 'DNA'
		awk '$11 ~ /DNA/ {print $0}' scaf/annotate.fas.out > DNA/DNA.out
		awk '{print $5}' DNA/DNA.out | sort | uniq > DNA/DNAscafs.lst
		wc -l DNA/DNAscafs.lst 
		cat DNA/DNAscafs.lst | cdbyank scaf/scaf.fas.cidx -o DNA/DNAscafs.fas 
		#echo 'TOTAL READS MAPPED'
		grep -f DNA/DNAscafs.lst scaf/scafreads.lst| awk '{s+=$3}END{print s}'
	echo 'Satellite'
		awk '$11 ~ /Satellite/ {print $0}' scaf/annotate.fas.out > other/Satellite.out
		awk '{print $5}' other/Satellite.out | sort | uniq > other/Satellitescafs.lst
		wc -l other/Satellitescafs.lst 
		cat other/Satellitescafs.lst | cdbyank scaf/scaf.fas.cidx -o other/Satellitescafs.fas 
		#echo 'TOTAL READS MAPPED'
		grep -f other/Satellitescafs.lst scaf/scafreads.lst| awk '{s+=$3}END{print s}'	
	echo 'rRNA'
		awk '$11 ~ /rRNA/ {print $0}' scaf/annotate.fas.out > other/rRNA.out
		awk '{print $5}' other/rRNA.out | sort | uniq > other/rRNAscafs.lst
		wc -l other/rRNAscafs.lst 
		cat other/rRNAscafs.lst | cdbyank scaf/scaf.fas.cidx -o other/rRNAscafs.fas 
		#echo 'TOTAL READS MAPPED'
		grep -f other/rRNAscafs.lst scaf/scafreads.lst| awk '{s+=$3}END{print s}'	