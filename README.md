AsparagalesTEscripts
====================

Methods from characterization of TEs from low-coverage genome survey sequences in exemplar Asparagales taxa, from [Hertweck 2013](http://www.nrcresearchpress.com/doi/abs/10.1139/gen-2013-0042). This repository is kept as a part of this publication; see other repos for updates to these methods.

**Dependencies:**
* MSR-CA (Maryland Short Read Assembler, now MaSuRCA)
* smalt
* samtools
* seqtk
* cdbyank/cdbfasta
* RepeatMasker

**REpipeline.sh** runs MSR-CA, maps reads to contigs, filters out organellar contigs, runs RepeatMasker on remaining contigs, and summarizes results.

**REDNAsuperfam.sh** parses results further into DNA transposon superfamilies.
