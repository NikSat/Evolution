#!/usr/bin/perl -w


use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my $Dir;
my $type="kir";
my $outfile;
my $file;

# I will use a hash and a reference to sort the results 

my $storage;

my $argc=scalar(@ARGV);  # The single argument


for(my $i=0; $i<$argc; $i++)
	{
  	my $entry=shift @ARGV;
	if($entry=~/^-dir/) {$i++; $Dir=shift @ARGV;}
	elsif($entry=~/^-hap/) {$i++; $type="hap";}
	elsif($entry=~/^-kir/) {$i++; $type="kir";}
  else {warn "\nunknown option: $entry\n"; exit;}
	}

if($argc==0)
	{
	showHelp();
	}

if ($Dir=~ /\/$/)
{
$Dir=substr($Dir, 0, -1); 
}

$outfile = "RSIprep".$type. $Dir .".txt";

# open the directory
opendir DIR, $Dir or die "Cannot open $Dir $!";

# now open the subdirectories in the directory

while (my $subdir= readdir(DIR)) {
	next if($subdir =~ /\D/);
	my $nname=$subdir;
	$nname=~ s/\D//g;
#	my $directory=$Dir. '/' .$subdir. '/';
#	opendir SUB_DIR, $directory or die "Cannot open subdirectory $subdir $!";
#	print $directory , "\n";
#	while (my $file = readdir(SUB_DIR))
#	{
#	next unless ($file =~ /SIAverage/);
#	$file=$directory . $file;
	$file="SI".$type."Average".$nname.".txt";
	print $file, "\n";
	open(FILE, $file) or die "Could not open $file $!";
	while (<FILE>) {
			next if($_ =~ /#/);
			my($spec, $pop, $sd, $chronic, $sd_chronic) = split /\t/, $_;
			push(@{$storage->{$spec}->{"pop"}},$pop);
			push(@{$storage->{$spec}->{"sd"}},$sd);
			push(@{$storage->{$spec}->{"chr"}},$chronic);
			push(@{$storage->{$spec}->{"sd chr"}},$sd_chronic);
		
			}
	close(FILE);

#}
#close(SUB_DIR);
#print Dumper($storage);
#print "$subdir Ok \n";
}


# open the output file
	open(OUT, ">$outfile") or die "could not open $outfile: $!";
	print OUT "# Run\tPopulation\n";
	
#print the contents of the hash in the correct order
foreach my  $name (sort {$a<=>$b} keys %{$storage})
{
#my $average=CalcAverage(@{$storage->{$name}->{"pop"}});
#my $poolSD=CalcSD(@{$storage->{$name}->{"pop"}});
#my $chroav=CalcAverage(@{$storage->{$name}->{"chr"}});
#my $chrpoolSD=CalcSD(@{$storage->{$name}->{"chr"}});


# $average, "\t",$poolSD, "\t", $chroav, "\t", $chrpoolSD, "\n";
foreach my $member (@{$storage->{$name}->{"pop"}})
{

print OUT $name,"\t" ,$member, "\n";

}

}

close(OUT);


sub showHelp
{
	die("This script collects all the Average population results and creates
A Directory containing the Directories with the results of the runs should be
given.This directories must have a distinct number from 1 to 10
Syntax: ./MultiAverage -dir DIRNAME \n")
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
$sum=($sum)**0.5;
my $SD=$sum/$N;
return $SD;
}


sub CalcPoolSD
{
my (@pool)=@_;
my $N=scalar(@pool);
my $sum=0;
foreach my $member (@pool){
$member=($member)**2;
$sum+=$member;
}
my $SD=$sum/$N;
$SD=($SD)**0.5;
return $SD;
}


