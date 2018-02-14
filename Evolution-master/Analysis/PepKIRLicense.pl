#!/usr/bin/perl -w

use Getopt::Long;
use Data::Dumper;
use warnings;
use strict;


#This script calculates the average lisencing per KIR  over time

my $dir;
my $poolfile;
my $outfile; 
my @mhc;
my @matrix;
my %map;
my $loci = 1;
my $mhc_loci = 1;
my $population=0;
my %genefiles;
my $assocfile="Association.data";
my @final;
my $mapfile="Sequence.map";
my %assoc;

my $argc=scalar(@ARGV);  # The arguments

for(my $i=0; $i<$argc; $i++)
{
  my $entry=shift @ARGV;
  if($entry=~/^-dir/) {$i++; $dir=shift @ARGV;}
  elsif($entry=~/^-loci_kir/) {$i++; $loci=shift @ARGV;} 
  elsif($entry=~/^-number_loci_mhc/) {$i++; $mhc_loci=shift @ARGV;}
  else {warn "\nunknown option: $entry\n"; exit;}
}

if ($argc==0)
{
    ShowHelp();
}

#Getting the data from the directory

$mapfile="./" . $dir . "/Sequence.map";
$assocfile="./" . $dir . "/Association.data";


#print $assocfile, "\n";

# Now load all the MHC's and the KIR-MHC association from the Association.data file


open(FILE3, "$assocfile") or die "Could not open $assocfile: $!";

while (<FILE3>)
	{	
		chomp;
		next if($_ =~ /#/);
		my ($kirmole, $mhcmole) = split /\s/, $_;
		#print $kirmole, "\t";
		#print $mhcmole, "\n";
		push(@mhc, $mhcmole);
		push (@{$assoc{$kirmole}},$mhcmole);
	}

#Remove the duplicates
#print @mhc, "\n";
@mhc= grep defined,@mhc;
foreach my $kk (keys %assoc)
{
@{$assoc{$kk}}=Unique(@{$assoc{$kk}});
}

@mhc=Unique(@mhc);
#print Dumper(@mhc), "\n";
#print @kir, "\n"; 
#print @kir, "\n";

#exit;
#@kir=Unique(@kir); 

#print Dumper(%assoc), "ok", "\n";
#exit;

#if ((scalar @mhc)<28) {
#print scalar @mhc," That's a bummer, sorry \n"; 
#exit;
#}
close (FILE3);


# First of all load the map
open(FILE, $mapfile) or die "Could not open $mapfile $!";
while (<FILE>) {

            next if($_ =~ /#/);
            my @line = split /\s/, $_;
            $map{$line[0]}=$line[1];
}

close(FILE);

# Open the outfile

$outfile="AVKIRLicensing_L=" . $dir . ".txt";

open(OUT, ">$outfile") or die "could not open $outfile: $!";


# Open the directory 
opendir DIR, $dir or die "Cannot open the $dir $!";

# Select all the files

    while (my $file = readdir DIR)
    {
        next unless $file =~/\.Genes.txt/;
        my @name = split /\./, $file;
	my $time = $name[0];
	$genefiles{$time}=$file;
    }
    closedir DIR;


foreach my $f (sort {$a<=>$b} keys %genefiles)
{

#Start looking at all the files

my @name = split /\./, ${genefiles{$f}};
my $time = $name[0];

my $fname="./" . $dir . "/" . ${genefiles{$f}};
# Read the contents of each file


#print "$time \n";

open(FILE, $fname) or die "Could not open $f $!";

# Local scope for the maps (will be used again and again)

my %association;

# First fill the association map with zeros
for (my $i=1;$i<15;$i++)
{
$association{$i}=0;
}


# First column Host ID (0), Second(1) to Fifth MHC(4), Sixth(5) to Fifteenth(14) Haplotype one, Sixteenth(15) to Twentyfifth(24) the Haplotype two NOTE TO SELF: These can be made variable for future use in the next project 

while (<FILE>) {

            next if($_ =~ /#/);
            my @line = split /\t/, $_;
            
# The first haplotype combined into one string of 80 Characters, this will be used as a key for the first map
# The first map will save computation since the same procedure will not take place for the second time          
# First map (%haplo): Haplotype and recognition
# Second map (%association): How many it recognizes - how many are there in the population
            
            my $molecule;
            my $rec=0;
	    my $column_m = $mhc_loci*2+1;
	    my $column = $loci*2+$mhc_loci*2+1;	
#The key is the KIR sequence

            for (my $i = $column_m; $i<$column; $i++)
            {
             $molecule=$map{$line[$i]}; 
    #Check if it already exists, if not do the process, else do not do the process again
		$rec=scalar(@{$assoc{$molecule}});			       
		#print $molecule, "\t", $rec, "\n"; 
		#exit;
            $association{$rec}+=1;
        #print $molecule;
        #print "\n";
    $population=$line[0] unless $population>=$line[0];  
$rec=0;

}
}


#print Dumper(%association), "\n";
#exit;
# Now write at the outfile
my $av;
print OUT $time, "\t", $population, "\t";
$population=$population*2*$loci;
for (my $i=1;$i<15;$i++)
{
my $test=(($association{$i}/$population)*$i);
#print $test ,"\t";
$av+=$test;

}
#print "\n";
print OUT $av, "\n";

# Close the file and go to the next one
close(FILE);
$population=0;
$av=0;
#print "end \n";


}


close(OUT);

sub Unique
{
    my %hash =map{ $_, 1 } @_;
    @_ = keys %hash;
    return @_;
}



sub ShowHelp
{
die("This script calculates the frequency of different KIR moplecules
grouped by the number of MHC's they recognize.
The syntax is :./KIRfrequency -dir <directory name>\n");
}
