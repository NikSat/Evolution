#!/usr/bin/perl -w


use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;


my $number;
my $Dir;
my $Years;
my $outfile; 

# I will use a hash and a reference to sort the results 

my %storage;
my %secondstorage;


my $argc=scalar(@ARGV);  # The arguments


for(my $i=0; $i<$argc; $i++)
	{
  	my $entry=shift @ARGV;
	if($entry=~/^-dir/) {$i++; $Dir=shift @ARGV;}
	elsif($entry=~/^-number/) {$i++; $number=shift @ARGV;}
  	else {warn "\nunknown option: $entry\n"; exit;}
	}

if($argc==0)
	{
	showHelp();
	}


$outfile="PresentedPeptidesAverage_L=".$number. ".txt";
# open the directory
opendir DIR, $Dir or die "Cannot open $Dir $!";




# now open the subdirectories in the directory

while (my $subdir= readdir(DIR)) {
	next if($subdir =~ /\D/);
	my $name=$subdir;
	$name=~ s/\D//g;
	my $directory=$Dir. '/' .$subdir;
	opendir SUB_DIR, $directory or die "Cannot open subdirectory $subdir $!";
	
	my $filename=$directory.'/AveragePepPresentation_L='.$number.'.txt';;
	open(FILE, $filename) or die "Could not open ", $filename ," in directory $subdir $!";
	while (<FILE>) {

			next if($_ =~ /#/);
			my($time, $self,$nonself) = split /\s/, $_;
			{
			push(@{$storage{$time}},$self);
			push(@{$secondstorage{$time}},$nonself);
			}
					}
	close(FILE);


close(SUB_DIR);
print "$subdir Ok \n";

}


#print Dumper(%storage), "\n";
#exit;

#calculate the SI 


#open the output file
	open(OUT, ">$outfile") or die "could not open $outfile: $!";
	print OUT "# Time\t Average Self\t SD\t Average NonSelf\t SD\n";
	
# print the contents of the hash in the correct order

foreach my  $key (sort {$a<=>$b} keys %storage)
{
my $aveSI=CalcAverage(@{%storage->{$key}});
my $sdSI=CalcSD(@{%storage->{$key}});
my $secondaveSI=CalcAverage(@{%secondstorage->{$key}});
my $secondsdSI=CalcSD(@{%secondstorage->{$key}});
print OUT $key, "\t", $aveSI, "\t" ,$sdSI, "\t", $secondaveSI, "\t" ,$secondsdSI,"\n";
}

close(OUT);


sub showHelp
{
	die("This script Calculates the mean of Presented Peptides in each Directory
This directories must have a distinct numbers refering to the mutant number\n")
}

sub CalcAverage
{
my (@pool)=@_;
my $N=scalar(@pool);
my $sum=0;
foreach my $member (@pool){
$sum+=$member;
}
my $av=$sum/$N;
return $av;
}


sub CalcSD
{
my (@pool)=@_;
my $N=scalar(@pool);
my $sum=0;
foreach my $member (@pool){
$sum+=$member;
}
my $av=$sum/$N;
foreach (@pool){
@_=(@_-$av)**2;
}
foreach my $member (@pool){
$sum+=$member;
}
$sum=($sum)**0.5;
my $SD=$sum/$N;
return $SD;
}
