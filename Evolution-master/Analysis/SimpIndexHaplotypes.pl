#!/usr/bin/perl -w

#use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
            
# This script calculates the SRI for KIRs, MHCs or KIR haplotypes 
 
my $freq=0;
my $sum;
my $gene_type = "mhc";
my $loci = 1;
my $mhc_loci = 1;
my $outfile;
my $argc=scalar(@ARGV);
my $dir;
my $tempoutfilename;
#print Dumper(@ARGV);
for(my $i=0; $i<$argc; $i++)
{
  my $entry=shift @ARGV;
  if($entry=~/^-dir/) {$i++; $dir=shift @ARGV;}
  elsif($entry=~/^-out/) {$i++; $tempoutfilename=shift @ARGV;}
  elsif($entry=~/^-mhc/) {$i++; $gene_type="mhc";}
  elsif($entry=~/^-kir/) {$i++; $gene_type="kir";}
  elsif($entry=~/^-haplotypes/) {$i++; $gene_type="hap";}
  elsif($entry=~/^-decoy/) {$i++; $gene_type="decoy";}
  elsif($entry=~/^-loci_kir/) {$i++; $loci=shift @ARGV;} 
  elsif($entry=~/^-number_loci_mhc/) {$i++; $mhc_loci=shift @ARGV;} 
  else {warn "\nunknown option: $entry\n"; exit;}
}

if ($tempoutfilename) {$outfile=$tempoutfilename;}
else {$outfile="Diversity_" . $gene_type . "_L=" . $dir . ".txt";}



if($argc==0)
{
	showHelp();
}



my $hash;
#my $gene_file=1;
my $histogram;
opendir DIR, $dir or die "Cannot open $dir $!";

if($gene_type eq 'mhc' || $gene_type eq 'kir')
{
	while (my $file = readdir DIR)
	{
		next unless $file =~/\.Genes.txt/;
		$hash->{$file}=1;
	}
	closedir DIR;

	open(OUT, ">$outfile") or die "could not open $outfile: $!";
	print OUT "# time Simpson's diversity\n";
	foreach my $f (sort {$a<=>$b} keys %{$hash})
	{
		$sum = 0;
		my @genes;
		my $total_pop = 0;
		my $gene_file = "$dir/$f";
		my @name = split /\./, $f;
		my $time = $name[0];
		open(FILE, "$gene_file") or die "Could not open $gene_file: $!";
	
	
	#extract the mhc info of the host population (columns 2 and 3) or the kir info (columns >3)
		while (<FILE>)
		{
			next if($_ =~ /#/);
			my @line = split /\t/, $_;
			my $column_m = $mhc_loci*2+1;
			
			if($gene_type eq "mhc")
			{
				for(my $m = 1; $m<$column_m; $m++)
				{
					push(@genes, $line[$m]);
				}	
			}
			if($gene_type eq "kir")
			{
				my $column = $loci*2+$mhc_loci*2+1;				
				for(my $i = $column_m; $i<$column; $i++)
				{
					push(@genes, $line[$i]);
				}
			}
		}
		close(FILE);
		$total_pop = scalar(@genes);
		#print $total_pop, "\n";
		$histogram->{$time}->{$_}++ for @genes;
		foreach my $nr (keys %{$histogram->{$time}})
		{
			$histogram->{$time}->{$nr}->{$freq} = $histogram->{$time}->{$nr}/$total_pop;
			$sum+= ($histogram->{$time}->{$nr}->{$freq})**2;
			#print $time, " ", $total_pop, " ",$nr, " ", $histogram->{$time}->{$nr}->{$freq}, " ",($histogram->{$time}->{$nr}->{$freq})**2, " ", $sum, "\n";
		}
	
		my $diversity_index = 0;
		$diversity_index = 1/$sum;
		print OUT $time, " ", $diversity_index,"\n";
		$#genes = -1;
	}
	close OUT;
#	print Dumper($histogram);
}

if($gene_type eq 'decoy')
{
	while (my $file = readdir DIR)
	{
		next if $file =~/Genes\.txt/;
		next if $file =~/Age\.txt/;
		next if $file =~/notBornChildren\.txt/;
		next if $file =~/PopulationSize\.txt/;
		next if $file =~/\.data/;
		$hash->{$file}=1;
	}
	closedir DIR;
#	print Dumper($hash);
	open(OUT, ">$outfile") or die "could not open $outfile: $!";
	print OUT "# time Simpson's diversity\n";
	foreach my $f (sort {$a<=>$b} keys %{$hash})
	{
		$sum = 0;
		my $total_pop = 0;
		my $gene_file = "$dir/$f";
		my @name = split /\./, $f;
		my $time = $name[0];
		my @genes;
		open(FILE, "$gene_file") or die "Could not open $gene_file: $!";
	
	#extract the mhc info of the host population (columns 2 and 3) or the kir info (columns >3)
		while (<FILE>)
		{
			next if($_ =~ /#/);
			my @line = split /\t/, $_;
			push(@genes, $line[10]) unless $line[10]== 0;
		}
		close(FILE);
		$total_pop = scalar(@genes);
		#print $total_pop, "\n";
		$histogram->{$time}->{$_}++ for @genes;
		#print join ",", @genes;
		foreach my $nr (keys %{$histogram->{$time}})
		{
			$histogram->{$time}->{$nr}->{$freq} = $histogram->{$time}->{$nr}/$total_pop;
			$sum+= ($histogram->{$time}->{$nr}->{$freq})**2;
	#		print $time, " ", $total_pop, " ",$nr, " ", $histogram->{$time}->{$nr}->{$freq}, " ",($histogram->{$time}->{$nr}->{$freq})**2, " ", $sum, "\n";
		}
#		print $sum , "\n";
		my $diversity_index = 0;
		if($sum == 0)
		{
			$diversity_index = 0;
		}
		else
		{
			$diversity_index = 1/$sum;
		}
		print OUT $time, " ", $diversity_index,"\n";
		$#genes = -1;
	}
	close OUT;

}


if($gene_type eq 'hap')
{
	while (my $file = readdir DIR)
	{
		next unless $file =~/\.Genes.txt/;
		$hash->{$file}=1;
	}
	closedir DIR;

	open(OUT, ">$outfile") or die "could not open $outfile: $!";
	print OUT "# time Simpson's diversity\n";
	foreach my $f (sort {$a<=>$b} keys %{$hash})
	{
		$sum = 0;
		my @genes;
		my $total_pop = 0;
		my $gene_file = "$dir/$f";
		my @name = split /\./, $f;
		my $time = $name[0];
		open(FILE, "$gene_file") or die "Could not open $gene_file: $!";
	
	#extract the mhc info of the host population (columns 2 and 3) or the kir info (columns >3)
		while (<FILE>)
		{
			next if($_ =~ /#/);
			my @line = split /\t/, $_;
			my $column_m = $mhc_loci*2+1;
			my $column = $loci+$mhc_loci*2+1;	
			my $secondcolumn=$loci*2+$mhc_loci*2+1;
			my $haplotype;
			my $secondhaplotype; 			
			for(my $i = $column_m; $i<$column; $i++)
				{
					$haplotype=$haplotype . $line[$i];
				}
			push(@genes, $haplotype);
			for(my $i = $column; $i<$secondcolumn; $i++)
			{
				$secondhaplotype=$secondhaplotype . $line[$i];
			}
			push(@genes, $secondhaplotype);
		}
		close(FILE);
		$total_pop = scalar(@genes);
		#print $total_pop, "\n";
		$histogram->{$time}->{$_}++ for @genes;
		foreach my $nr (keys %{$histogram->{$time}})
		{
			$histogram->{$time}->{$nr}->{$freq} = $histogram->{$time}->{$nr}/$total_pop;
			$sum+= ($histogram->{$time}->{$nr}->{$freq})**2;
			#print $time, " ", $total_pop, " ",$nr, " ", $histogram->{$time}->{$nr}->{$freq}, " ",($histogram->{$time}->{$nr}->{$freq})**2, " ", $sum, "\n";
		}
	
		my $diversity_index = 0;
		$diversity_index = 1/$sum;
		print OUT $time, " ", $diversity_index,"\n";
		$#genes = -1;
	}
	close OUT;
#	print Dumper($histogram);
}



		#print Dumper ($histogram);
sub showHelp
{
	die("This tool calculates the Simpson's Diversity Index. \nUsage: simsons_index.pl [options]\n
	Options\n:
	-dir (directory where the *Genes.txt files are to be found)
	-out (outfile, optional, default is Diversity_gene_type_directory.txt)
	-mhc (calculate MHC SI)
	-kir (calculate KIR SI)
	-hap (calculate the Haplotypes SI)
	-loci_kir
	-number_loci_mhc
	-decoy (calculate decoy SI)
	\n");
}
