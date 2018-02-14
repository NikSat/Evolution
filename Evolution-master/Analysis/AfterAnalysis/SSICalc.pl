#!/usr/bin/perl -w


use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my $Dir;
my $type="kir";
my $Years;


# I will use a hash and a reference to sort the results 



my $argc=scalar(@ARGV);  # The arguments


for(my $i=0; $i<$argc; $i++)
	{
  	my $entry=shift @ARGV;
	if($entry=~/^-years/) {$i++; $Years=shift @ARGV;}
	elsif($entry=~/^-hap/) {$i++; $type="hap";}
	elsif($entry=~/^-kir/) {$i++; $type="kir";}
	elsif($entry=~/^-dir/) {$i++; $Dir=shift @ARGV;}
  else {warn "\nunknown option: $entry\n"; exit;}
	}

if($argc==0)
	{
	showHelp();
	}

# open the directory
opendir DIR, $Dir or die "Cannot open $Dir $!";

# now open the subdirectories in the directory

while (my $subdir= readdir(DIR)) {
	next if($subdir =~ /\D/);
	my $name=$subdir;
	$name=~ s/\D//g;
	my $directory=$Dir. '/' .$subdir;
	opendir SUB_DIR, $directory or die "Cannot open subdirectory $subdir $!";
	my %storage;
	for (my $m=1;$m<11;$m++)
	{ 
	my $filename=$directory.'/Diversity_'.$type.'_L='.$m.'.txt';
	open(FILE, $filename) or die "Could not open $filename in directory $subdir $!";
	my @SI;
	while (<FILE>) {

			next if($_ =~ /#/);
			chomp;
			my @line = split /\s/, $_;
			#print $line[0], "\t";
			#print $line[1], "\n";
			unless ($line[0]<$Years) 
						{
						push(@SI,$line[1]);
						}
					}
	close(FILE);
#print Dumper(@SI);
#exit;
	my $Average=CalcAverage(@SI);
	my $SISD=CalcSD(@SI);


$storage{$m}=[$Average,$SISD];
#print Dumper(%storage);
#exit;
}
my $outfile='SI'.$type.'Average'.$subdir.'.txt';

#open the output file
	open(OUT, ">$outfile") or die "could not open $outfile: $!";
	print OUT "# Run\tAverage SI\tSD\n";
	
# print the contents of the hash in the correct order
foreach my  $name (sort {$a<=>$b} keys %storage)
{

print OUT $name, "\t", $storage{$name}[0], "\t",$storage{$name}[1], "\n";

}

close(OUT);

print "$subdir Ok \n";


close(SUB_DIR);


}

sub showHelp
{

die("This script Calculates the average SI for the last 
variable number of years. The number of (Total years - years of calculation must be given).
A Directory containing the Directories with the results of the runs should be given. 
This directories must have a distinct number from 1 to 10 \n");

}

sub CalcAverage
{
my (@pool)=@_;
my $N=scalar(@pool);
my $sum;
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
my $sum;
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

