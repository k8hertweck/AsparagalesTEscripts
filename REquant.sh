#!/bin/bash
#$ -S /bin/bash -cwd
#$ -o TAXON.out -j y
#$ -M k8hertweck@gmail.com -m be
#$ -l highprio
#$ -N TAXON

##MUST HAVE PROGRAMS INSTALLED:
##RepeatMasker [smalt, samtools, cdbyank]
##THIS FILE IN TAXON FOLDER; RESULTS FILE ALREADY PRESENT: CHANGE TAXON AND READS

cd /Users/kate/Desktop/REdata/TAXON/annotate

#UNIQUE SCAFFOLDS
	echo 'UNIQUE SCAFFOLDS' | tee -a quant.out
	mkdir uniquescafs
	cd uniquescafs
	tail +4 ../scafNuc.fas.out | awk '{print $5}' | sort | uniq > uniquescafs.lst
	wc -l uniquescafs.lst | tee -a ../quant.out
	
	#PULL OUT UNIQUE SCAFFOLDS WITH SEQUENCES FROM scafNuc.fas
	cdbfasta ../RM/scafNuc.fas
	cat uniquescafs.lst | cdbyank ../RM/scafNuc.fas.cidx -o uniquescafs.fas 
	
	#MAP READS TO UNIQUE SEQUENCES
	smalt index uniquescafs uniquescafs.fas | tee -a ../quant.out
	smalt map -f sam -o uniquescafs.sam uniquescafs ~/data/TAXONTRIM.fastq | tee -a ../quant.out
	
	##convert from SAM to BAM
	samtools view -bS -o uniquescafs.bam uniquescafs.sam 
	##sort BAM
	samtools sort uniquescafs.bam uniquescafssort
	##summary stats: how many sequences for each reference
	samtools index uniquescafssort.bam
	samtools idxstats uniquescafssort.bam > uniquescafsreads.lst
	
	#SUM COLUMNS AND REPORT STATS	
	echo 'TOTAL SCAFFOLD LENGTH' uniquescafsreads.lst | tee -a ../quant.out
	awk '{s+=$2}END{print s}' uniquescafsreads.lst | tee -a ../quant.out
	echo 'TOTAL READS MAPPED' uniquescafsreads.lst | tee -a ../quant.out
	awk '{s+=$3}END{print s}' uniquescafsreads.lst | tee -a ../quant.out
	echo 'READS UNMAPPED' uniquescafsreads.lst | tee -a ../quant.out
	awk '{s+=$4}END{print s}' uniquescafsreads.lst | tee -a ../quant.out
	
	cd ..
	
#RETROELEMENTS
echo 'RETROELEMENTS' | tee -a quant.out
	echo 'LTRs' | tee -a quant.out		
		cd retroelements/LTRs
		awk '{print $5}' LTR.out | sort | uniq > LTRscafs.lst
		wc -l LTRscafs.lst | tee -a ../../quant.out
	
		#PULL OUT UNIQUE SCAFFOLDS WITH SEQUENCES FROM scafNuc.fas
		cat LTRscafs.lst | cdbyank ../../RM/scafNuc.fas.cidx -o LTRscafs.fas 
		
		#MAP READS TO UNIQUE SEQUENCES
		smalt index LTRscafs LTRscafs.fas | tee -a ../../quant.out
		smalt map -f sam -o LTRscafs.sam LTRscafs ~/data/TAXONTRIM.fastq | tee -a ../../quant.out
		
		##convert from SAM to BAM
		samtools view -bS -o LTRscafs.bam LTRscafs.sam 
		##sort BAM
		samtools sort LTRscafs.bam LTRscafssort
		##summary stats: how many sequences for each reference
		samtools index LTRscafssort.bam
		samtools idxstats LTRscafssort.bam > LTRscafsreads.lst
		
		#SUM COLUMNS AND REPORT STATS	
		echo 'TOTAL SCAFFOLD LENGTH' LTRscafsreads.lst | tee -a ../../quant.out
		awk '{s+=$2}END{print s}' LTRscafsreads.lst | tee -a ../../quant.out
		echo 'TOTAL READS MAPPED' LTRscafsreads.lst | tee -a ../../quant.out
		awk '{s+=$3}END{print s}' LTRscafsreads.lst | tee -a ../../quant.out
		echo 'READS UNMAPPED' LTRscafsreads.lst | tee -a ../../quant.out
		awk '{s+=$4}END{print s}' LTRscafsreads.lst | tee -a ../../quant.out

	echo 'Copia' | tee -a ../../quant.out
		awk '{print $5}' Copia.out | sort | uniq > Copiascafs.lst
		wc -l Copiascafs.lst | tee -a ../../quant.out
	
		#PULL OUT UNIQUE SCAFFOLDS WITH SEQUENCES FROM scafNuc.fas
		cat Copiascafs.lst | cdbyank ../../RM/scafNuc.fas.cidx -o Copiascafs.fas 
		
		#MAP READS TO UNIQUE SEQUENCES
		smalt index Copiascafs Copiascafs.fas | tee -a ../../quant.out
		smalt map -f sam -o Copiascafs.sam Copiascafs ~/data/TAXONTRIM.fastq | tee -a ../../quant.out
		
		##convert from SAM to BAM
		samtools view -bS -o Copiascafs.bam Copiascafs.sam 
		##sort BAM
		samtools sort Copiascafs.bam Copiascafssort
		##summary stats: how many sequences for each reference
		samtools index Copiascafssort.bam
		samtools idxstats Copiascafssort.bam > Copiascafsreads.lst
		
		#SUM COLUMNS AND REPORT STATS	
		echo 'TOTAL SCAFFOLD LENGTH' Copiascafsreads.lst | tee -a ../../quant.out
		awk '{s+=$2}END{print s}' Copiascafsreads.lst | tee -a ../../quant.out
		echo 'TOTAL READS MAPPED' Copiascafsreads.lst | tee -a ../../quant.out
		awk '{s+=$3}END{print s}' Copiascafsreads.lst | tee -a ../../quant.out
		echo 'READS UNMAPPED' Copiascafsreads.lst | tee -a ../../quant.out
		awk '{s+=$4}END{print s}' Copiascafsreads.lst | tee -a ../../quant.out

	echo 'Gypsy' | tee -a ../../quant.out
		awk '{print $5}' Gypsy.out | sort | uniq > Gypsyscafs.lst
		wc -l Gypsyscafs.lst | tee -a ../../quant.out
	
		#PULL OUT UNIQUE SCAFFOLDS WITH SEQUENCES FROM scafNuc.fas
		cat Gypsyscafs.lst | cdbyank ../../RM/scafNuc.fas.cidx -o Gypsyscafs.fas 
		
		#MAP READS TO UNIQUE SEQUENCES
		smalt index Gypsyscafs Gypsyscafs.fas | tee -a ../../quant.out
		smalt map -f sam -o Gypsyscafs.sam Gypsyscafs ~/data/TAXONTRIM.fastq | tee -a ../../quant.out
		
		##convert from SAM to BAM
		samtools view -bS -o Gypsyscafs.bam Gypsyscafs.sam 
		##sort BAM
		samtools sort Gypsyscafs.bam Gypsyscafssort
		##summary stats: how many sequences for each reference
		samtools index Gypsyscafssort.bam
		samtools idxstats Gypsyscafssort.bam > Gypsyscafsreads.lst
		
		#SUM COLUMNS AND REPORT STATS	
		echo 'TOTAL SCAFFOLD LENGTH' Gypsyscafsreads.lst | tee -a ../../quant.out
		awk '{s+=$2}END{print s}' Gypsyscafsreads.lst | tee -a ../../quant.out
		echo 'TOTAL READS MAPPED' Gypsyscafsreads.lst | tee -a ../../quant.out
		awk '{s+=$3}END{print s}' Gypsyscafsreads.lst | tee -a ../../quant.out
		echo 'READS UNMAPPED' Gypsyscafsreads.lst | tee -a ../../quant.out
		awk '{s+=$4}END{print s}' Gypsyscafsreads.lst | tee -a ../../quant.out
	
		cd ..

	echo 'SINE' | tee -a ../quant.out
		awk '{print $5}' SINEs.out | sort | uniq > SINEscafs.lst
		wc -l SINEscafs.lst | tee -a ../quant.out
	
		#PULL OUT UNIQUE SCAFFOLDS WITH SEQUENCES FROM scafNuc.fas
		cat SINEscafs.lst | cdbyank ../RM/scafNuc.fas.cidx -o SINEscafs.fas 
		
		#MAP READS TO UNIQUE SEQUENCES
		smalt index SINEscafs SINEscafs.fas | tee -a ../quant.out
		smalt map -f sam -o SINEscafs.sam SINEscafs ~/data/TAXONTRIM.fastq | tee -a ../quant.out
		
		##convert from SAM to BAM
		samtools view -bS -o SINEscafs.bam SINEscafs.sam 
		##sort BAM
		samtools sort SINEscafs.bam SINEscafssort
		##summary stats: how many sequences for each reference
		samtools index SINEscafssort.bam
		samtools idxstats SINEscafssort.bam > SINEscafsreads.lst
		
		#SUM COLUMNS AND REPORT STATS	
		echo 'TOTAL SCAFFOLD LENGTH' SINEscafsreads.lst | tee -a ../quant.out
		awk '{s+=$2}END{print s}' SINEscafsreads.lst | tee -a ../quant.out
		echo 'TOTAL READS MAPPED' SINEscafsreads.lst | tee -a ../quant.out
		awk '{s+=$3}END{print s}' SINEscafsreads.lst | tee -a ../quant.out
		echo 'READS UNMAPPED' SINEscafsreads.lst | tee -a ../quant.out
		awk '{s+=$4}END{print s}' SINEscafsreads.lst | tee -a ../quant.out

	echo 'LINE' | tee -a ../quant.out
		awk '{print $5}' LINEs.out | sort | uniq > LINEscafs.lst
		wc -l LINEscafs.lst | tee -a ../quant.out
	
		#PULL OUT UNIQUE SCAFFOLDS WITH SEQUENCES FROM scafNuc.fas
		cat LINEscafs.lst | cdbyank ../RM/scafNuc.fas.cidx -o LINEscafs.fas 
		
		#MAP READS TO UNIQUE SEQUENCES
		smalt index LINEscafs LINEscafs.fas | tee -a ../quant.out
		smalt map -f sam -o LINEscafs.sam LINEscafs ~/data/TAXONTRIM.fastq | tee -a ../quant.out
		
		##convert from SAM to BAM
		samtools view -bS -o LINEscafs.bam LINEscafs.sam 
		##sort BAM
		samtools sort LINEscafs.bam LINEscafssort
		##summary stats: how many sequences for each reference
		samtools index LINEscafssort.bam
		samtools idxstats LINEscafssort.bam > LINEscafsreads.lst
		
		#SUM COLUMNS AND REPORT STATS	
		echo 'TOTAL SCAFFOLD LENGTH' LINEscafsreads.lst | tee -a ../quant.out
		awk '{s+=$2}END{print s}' LINEscafsreads.lst | tee -a ../quant.out
		echo 'TOTAL READS MAPPED' LINEscafsreads.lst | tee -a ../quant.out
		awk '{s+=$3}END{print s}' LINEscafsreads.lst | tee -a ../quant.out
		echo 'READS UNMAPPED' LINEscafsreads.lst | tee -a ../quant.out
		awk '{s+=$4}END{print s}' LINEscafsreads.lst | tee -a ../quant.out
			
		cd ..

#DNA TRANSPOSONS
	echo 'DNA' | tee -a quant.out
	cd DNA	
	awk '{print $5}' DNA.out | sort | uniq > DNAscafs.lst
	wc -l DNAscafs.lst | tee -a ../quant.out

	#PULL OUT UNIQUE SCAFFOLDS WITH SEQUENCES FROM scafNuc.fas
	cat DNAscafs.lst | cdbyank ../RM/scafNuc.fas.cidx -o DNAscafs.fas 
	
	#MAP READS TO UNIQUE SEQUENCES
	smalt index DNAscafs DNAscafs.fas | tee -a ../quant.out
	smalt map -f sam -o DNAscafs.sam DNAscafs ~/data/TAXONTRIM.fastq | tee -a ../quant.out
	
	##convert from SAM to BAM
	samtools view -bS -o DNAscafs.bam DNAscafs.sam 
	##sort BAM
	samtools sort DNAscafs.bam DNAscafssort
	##summary stats: how many sequences for each reference
	samtools index DNAscafssort.bam
	samtools idxstats DNAscafssort.bam > DNAscafsreads.lst
	
	#SUM COLUMNS AND REPORT STATS	
	echo 'TOTAL SCAFFOLD LENGTH' DNAscafsreads.lst | tee -a ../quant.out
	awk '{s+=$2}END{print s}' DNAscafsreads.lst | tee -a ../quant.out
	echo 'TOTAL READS MAPPED' DNAscafsreads.lst | tee -a ../quant.out
	awk '{s+=$3}END{print s}' DNAscafsreads.lst | tee -a ../quant.out
	echo 'READS UNMAPPED' DNAscafsreads.lst | tee -a ../quant.out
	awk '{s+=$4}END{print s}' DNAscafsreads.lst | tee -a ../quant.out
	
	cd ..
	
#OTHER
echo 'OTHER' | tee -a quant.out
	cd other

echo 'ROLLING CIRCLES' | tee -a ../quant.out
	awk '{print $5}' RC.out | sort | uniq > RCscafs.lst
	wc -l RCscafs.lst | tee -a ../quant.out

	#PULL OUT UNIQUE SCAFFOLDS WITH SEQUENCES FROM scafNuc.fas
	cat RCscafs.lst | cdbyank ../RM/scafNuc.fas.cidx -o RCscafs.fas 
	
	#MAP READS TO UNIQUE SEQUENCES
	smalt index RCscafs RCscafs.fas | tee -a ../quant.out
	smalt map -f sam -o RCscafs.sam RCscafs ~/data/TAXONTRIM.fastq | tee -a ../quant.out
	
	##convert from SAM to BAM
	samtools view -bS -o RCscafs.bam RCscafs.sam 
	##sort BAM
	samtools sort RCscafs.bam RCscafssort
	##summary stats: how many sequences for each reference
	samtools index RCscafssort.bam
	samtools idxstats RCscafssort.bam > RCscafsreads.lst
	
	#SUM COLUMNS AND REPORT STATS	
	echo 'TOTAL SCAFFOLD LENGTH' RCscafsreads.lst | tee -a ../quant.out
	awk '{s+=$2}END{print s}' RCscafsreads.lst | tee -a ../quant.out
	echo 'TOTAL READS MAPPED' RCscafsreads.lst | tee -a ../quant.out
	awk '{s+=$3}END{print s}' RCscafsreads.lst | tee -a ../quant.out
	echo 'READS UNMAPPED' RCscafsreads.lst | tee -a ../quant.out
	awk '{s+=$4}END{print s}' RCscafsreads.lst | tee -a ../quant.out

	
echo 'SIMPLE REPEATS EXEMPT' 
echo 'SATELLITES EXEMPT' 
echo 'rRNA EXEMPT' 
echo 'LOW COMPLEXITY EXEMPT' 