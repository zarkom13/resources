#!/usr/bin/perl -w 

# given a set of somatic mutations and a fasta file with the reference genome 
# generate data suitable for an R plot of the mutational signature data, as in 
# this paper
# http://www.nature.com/nature/journal/vaop/ncurrent/full/nature12477.html
# and Fig S2 here:
# http://www.nature.com/nature/journal/vaop/ncurrent/extref/nature12477-s1.pdf

use strict; 
use Getopt::Long;
use Bio::SeqIO; 

my $mutations;
my $reference;
my $debug;

my $usage = "$0 
--mutations [mutations file, a tab-delimited file with chromosome\\tposition\treference residue\tmutated residue data]
--reference [reference fasta file]
";

GetOptions(
    'mutations=s' => \$mutations,
    'reference=s' => \$reference,
    'debug=s' => \$debug
    );

# check the things...
unless (defined $mutations){
    die "No --mutations argument specified\n$usage\n";
}
unless (defined $reference){
    die "No --reference argument specified\n$usage\n";
}
unless (-r $mutations){
    die "Can't read the file you specified as --mutations argument: $mutations";
}
unless (-r $reference){
    die "Can't read the file you specified as --reference argument: $reference";
}

open (my $MUTATIONS, "$mutations") || die "Can't open --mutations file $mutations: $!";

print STDERR "parsing reference sequence...\n";
my $refIO = Bio::SeqIO->new('-file' => $reference);
my $refMap; # a hash with each reference sequence ID (key) and its sequence (value)

while ( my $ref = $refIO->next_seq ){
    $refMap->{$ref->primary_id} = $ref->seq;
}
undef;
print STDERR "done.";

# the following is a hash, organized thusly:
# $signatureHistogramData->{mutation}->{3' residue}->{5' residue}
my $signatureHistogramData;
print STDERR "parsing somatic mutations to make histogram data for signature plot...\n";
while ( my $line = <$MUTATIONS> ){
    next if $line =~ m/^\#/;
    my @fields = split /\t/, $line; 
    unless ( defined $fields[0] && defined $fields[1] ){
	warn "field 1 undefined! skipping this line.\n";
	next;
    }
    # look up this mutated residue, and the 3' and 5' base...
    unless ( exists $refMap->{$fields[0]} ){
	warn "Can't find reference sequence $fields[0]!\n";
	next;
    }

    my $mutated = uc($fields[3]);
    my $thisRef = uc(substr( $refMap->{$fields[0]}, $fields[1] - 1, 1));
    my $threePrime = uc(substr( $refMap->{$fields[0]}, $fields[1] - 2, 1));
    my $fivePrime = uc(substr( $refMap->{$fields[0]}, $fields[1], 1));

    if ( $thisRef ne uc($fields[2]) ){
	warn "the item in the third column ($fields[0]) is not the same as what's in the reference sequence ($thisRef)!\n";
    }

    # first, if $thisRef is a G or A, we count this as a C or T, respectively. 
    # This is because the data in that paper are organized according to the pyrimidine of the
    # base-pair that was mutated
    if ( $thisRef eq 'A' or $thisRef eq 'G' ){ # count this as a T, reverse complement $thisRef and $mutated
	warn "before: reference is $thisRef, mutated is $mutated\n" if $debug;
	$mutated =~ tr/ACGT/TGCA/;
	$thisRef =~ tr/ACGT/TGCA/;
	warn "after: reference is $thisRef, mutated is $mutated\n" if $debug;
    }

    # increment the relevant histogram datum
    my $thisMutation = $thisRef . "->" . $mutated;
    $signatureHistogramData->{$thisMutation}->{$threePrime}->{$fivePrime}++;    
}
print STDERR "done.";

print STDERR "printing out histogram data...\n";
foreach my $thisMutation ( sort keys %{$signatureHistogramData} ){
    foreach my $threePrime ( "A", "C", "G", "T" ){
	foreach my $fivePrime ( "A", "C", "G", "T" ){
	    print $thisMutation . "\t" . $threePrime . "\t" . $fivePrime . "\t";
	    if ( exists $signatureHistogramData->{$thisMutation}->{$threePrime}->{$fivePrime} ){
		print $signatureHistogramData->{$thisMutation}->{$threePrime}->{$fivePrime};
	    } else {
		print 0;
	    }
	    print "\n";
	}
    }
}


print STDERR "done.";
