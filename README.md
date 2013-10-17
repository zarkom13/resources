Mutational-signature-analysis
=============================

Makes the data and plot similar to Fig 2 of this paper: 
http://www.nature.com/nature/journal/vaop/ncurrent/full/nature12477.htm

==== Installation

Install prerequisites: 

	       R (tested on R version 2.14.1)
	       Perl (tested on 5.14)
	       Statistics::R
	       Getopt::Long
	       Bio::SeqIO
	       ggplot2

1) install R - follow instructions here:
http://cran.r-project.org/
Linux: 
http://cran.r-project.org/bin/linux/

Mac: 
http://cran.r-project.org/bin/macosx/

Windows: 
http://cran.r-project.org/bin/windows/

2) install Perl modules: 
Using these instructions: 
http://www.cpan.org/modules/INSTALL.html

install these Perl modules: 
Statistics::R Getopt::Long Bio::SeqIO

3) install ggplot2, using these instructions: 
http://math.usask.ca/~longhai/software/installrpkg.html

You can probably just start R then type this: 
install.packages("ggplot2")

==== Running

./generate_signature_data.pl --mutations [mutations file] --reference [reference fasta file]

The mutations file should be a tab delimited file such as that produced by Varscan: 

chr1	13418	G	A	11	0	0%	G	8	4	33.33%	R	Somatic	1.0	0.055900621118012056	4	4	4	0	6	5	0	0

But this script only cares about the first 4 columsn: 
chromosome\tposition\treference nucleotide\tmutated nucleotide

Run TEST.sh to make sure things are working. 