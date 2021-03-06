#!/usr/bin/perl -w


use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;


#This script calculates the average recognized self and viral pMHC complexes over time
             
my $file;
my $sum;
my $nsum;
my $loci = 1;
my $mhc_loci = 1;
my $outfile;
my $argc=scalar(@ARGV);
my $thresshold;
my $dir;
my $tempoutfilename;
my $mapfile;
my %map;
my @matrix;
my $matrixfile;
my $poolfile;
my @mhc;
my @kir;
my @self;
my @nonself;
my $selfpMHC;
my $nonselfpMHC;
my $selfcheck;
my $nonselfcheck;
my $library="CMFILVWYAGTSQNEDHRKP";
my $outheatmap;
my $totalnonselfcheck;
my $totalselfcheck;
my $skir;
my $nskir;



#print Dumper(@ARGV);

for(my $i=0; $i<$argc; $i++)
{
  my $entry=shift @ARGV;
  if($entry=~/^-dir/) {$i++; $dir=shift @ARGV;}
  elsif($entry=~/^-out/) {$i++; $tempoutfilename=shift @ARGV;}
  elsif($entry=~/^-loci_kir/) {$i++; $loci=shift @ARGV;} 
  elsif($entry=~/^-number_loci_mhc/) {$i++; $mhc_loci=shift @ARGV;} 
  else {warn "\nunknown option: $entry\n"; exit;}
}

if ($tempoutfilename) {$outfile=$tempoutfilename;}
else {$outfile="Collected_L=" . $dir . ".txt";}

$mapfile="./" . $dir . "/Sequence.map";
$poolfile="./" . $dir . "/GenePool.data";
$matrixfile="./" . $dir . "/JMMATRIXCM.csv";
$file="./" . $dir . "/200000.Backup.data";
$outfile="AveragePepPresentation_L=".$dir.".txt";

if($argc==0)
{
	showHelp();
}


# Creating the multidimensional array

for (my $i=0;$i<20;$i++)
{
    $matrix[$i]=[];
}


my %dirthress =(
1 => -620,
2 => -256,
3 => -106,
4 => -12,
5 => 61,
6 => 117,
7 => 173,
8 => 204,
9 => 351,
10 => 373,
);

$thresshold=$dirthress{$dir};

#print $thresshold;
#exit;


# First of all load the map
open(FILE, $mapfile) or die "Could not open $mapfile $!";
while (<FILE>) {

            next if($_ =~ /#/);
            my @line = split /\s/, $_;
            $map{$line[0]}=$line[1];
}

close(FILE);

# Now to load the matrix for use in the second part of the script

open(FILE, $matrixfile) or die "Could not open $matrixfile $!";
while (<FILE>) {
    
    # Get the library from the first line
            my $number=$.-2;
            if($.==1) 
            {
                my @line= split /\s/, $_;
                for (my $k=0;$k<20;$k++)
                {
                    $library=$library . $line[$k];
                }
            }                       
            my @line= split /\s/, $_;
    # Throw away the first element 
            @line=@line[1..$#line];
            for (my $k=0;$k<20;$k++)
            {
                $matrix[$number][$k]=$line[$k];
            }
}

close(FILE);

#print Dumper(@matrix);
#exit;

# Now get all the genes, self and non self peptides

open(FILE, "$file") or die "Could not open $file: $!";

#get the first two lines 


my $line1= <FILE>;
my $line2= <FILE>;
my $line3= <FILE>;	

#the first line contains the MHCs and KIRs 

my @spline1=split /\t/, $line1;


my $kirnum=$spline1[1];
my $mhcnum=$kirnum/10;

# now collect the MHCs
for (my $in=2; $in<($mhcnum+2); $in++)
{
push (@mhc, $map{$spline1[$in]});  
}
# now collect the KIRs
for (my $im=($mhcnum+2); $im< ($mhcnum+$kirnum+2) ; $im++)
{
push (@kir, $map{$spline1[$im]});
}

# Now get the self peptides in the second line

my @spline2=split /\t/, $line2;

for (my $is=0; $is<1000; $is++)
{
push (@self, $map{$spline2[$is]});  
}

#Now collect the non self
my @spline3=split /\t/, $line3;

for (my $is=1; $is<51; $is++)
{
push (@nonself, $map{$spline3[$is]});  
}


close(FILE);


#print Dumper(@mhc);
#print Dumper(@kir);
#print Dumper(@self);
#print Dumper(@nonself);
#exit; 

# First create for all MHCs all the pMHC complexes

foreach my $tempmhc (@mhc)
{
 	foreach my $temppep (@self)
		{	
		#print $tempmhc, "\n";
		#print $temppep, "\n";
			my $bindscore=givescore($temppep,$tempmhc);
		#print $bindscore, "\n";	
			if ($bindscore>=121) {
						my $presented=present($tempmhc,$temppep); 
		#print $presented, "\n";
		#exit;
						push (@{$selfpMHC->{$tempmhc}},$presented);
					}
		}
 	foreach my $tempnonpep (@nonself)
		{	
			my $bindscore=givescore($tempnonpep,$tempmhc);
			if ($bindscore>=121) {
						my $presented=present($tempmhc,$tempnonpep); 
						push (@{$nonselfpMHC->{$tempmhc}},$presented);
					}
		}

} 

#Find how many pMHCs each KIR recognizes (store numbers, not ratios)

foreach my $totest (@kir)
	{
		my $spepsum=0;
		my $nspepsum=0;
		foreach my $tobetested (keys %{$selfpMHC})
		{
		my $subsum=scalar (@{$selfpMHC->{$tobetested}});
		foreach my $toagainbetested (@{$selfpMHC->{$tobetested}})
				{
				my $index=getscore($totest,$toagainbetested);
				if ($index>=$thresshold) 
					{
					$totalselfcheck->{$totest}->{$tobetested}++
					} 
				}

		$totalselfcheck->{$totest}->{$tobetested}=0.0 unless ($totalselfcheck->{$totest}->{$tobetested});
		$spepsum+=$totalselfcheck->{$totest}->{$tobetested};
		}
		$skir->{$totest}=$spepsum;
		foreach my $ntobetested (keys %{$nonselfpMHC})
		{
		my $nsubsum=scalar (@{$nonselfpMHC->{$ntobetested}});
		foreach my $ntoagainbetested (@{$nonselfpMHC->{$ntobetested}})
				{
				my $nindex=getscore($totest,$ntoagainbetested);
				if ($nindex>=$thresshold) 
					{
					$totalnonselfcheck->{$totest}->{$ntobetested}++
					} 
				}
		$totalnonselfcheck->{$totest}->{$ntobetested}=0.0 unless ($totalnonselfcheck->{$totest}->{$ntobetested});
		$nspepsum+=$totalnonselfcheck->{$totest}->{$ntobetested};
		}
		$nskir->{$totest}=$nspepsum;
	}




#print Dumper($selfpMHC);
#print Dumper(@mhc);
#print Dumper($nonselfpMHC);
#print Dumper($totalselfcheck);
#print Dumper($totalnonselfcheck);
#print Dumper($skir);
#print Dumper($nskir);
#exit;

#Select all genes from the gene files and calculate their frequences

my $hash;
#my $gene_file=1;
my $histogram;
my $haphistogram;
opendir DIR, $dir or die "Cannot open $dir $!";

while (my $file = readdir DIR)
	{
		next unless $file =~/\.Genes.txt/;
		my $tempname=(split /\./, $file)[0];
		$hash->{$tempname}=1;
	}
#print Dumper($hash);
#exit;
	closedir DIR;


	open(OUT, ">$outfile") or die "could not open $outfile: $!";
	print OUT "#Time\tSelf\tNonSelf\n";
	foreach my $f (sort {$a<=>$b} keys %{$hash})
	{
		$sum = 0;
		$nsum= 0;
		my @genes;
		my $total_pop = 0;
		my $gene_file = $dir.'/'.$f.'.Genes.txt';
		#print $gene_file, "\n";
		my $time = $f;
		#print $time, "\n";
		open(FILE, "$gene_file") or die "Could not open $gene_file: $!";
		
		#exit;	
	#extract the mhc info of the host population (columns 2 and 3) or the kir info (columns >3)
		while (<FILE>)
		{
			next if($_ =~ /#/);
			my @line = split /\t/, $_;
			my $column_m = $mhc_loci*2+1;
			my $column = $loci+$mhc_loci*2+1;	
			my $secondcolumn=$loci*2+$mhc_loci*2+1;
			for(my $i = $column_m; $i<$column; $i++)
				{
					push(@genes, $map{$line[$i]});
				}
			for(my $i = $column; $i<$secondcolumn; $i++)
				{
					push(@genes, $map{$line[$i]});
				}				
		}
		close(FILE);
		$total_pop = scalar(@genes);
		#print $total_pop, "\n";
		$histogram->{$time}->{$_}++ for @genes;


		foreach my $nr (keys %{$histogram->{$time}})
		{
			$histogram->{$time}->{$nr}= $histogram->{$time}->{$nr}/$total_pop;
			$sum+= ($histogram->{$time}->{$nr})*($skir->{$nr});
			$nsum+= ($histogram->{$time}->{$nr})*($nskir->{$nr});
			#print $time, " ", $total_pop, " ",$nr, " ", $histogram->{$time}->{$nr}->{$freq}, " ",($histogram->{$time}->{$nr}->{$freq})**2, " ", $sum, "\n";
		}
		print OUT $time,"\t", $sum,"\t",$nsum ,"\n";
		#print $time,"\t", $sum,"\t",$nsum ,"\n";

	}
close OUT;

#print Dumper(@selectedkir);
#exit;



sub Unique
{
    my %hash =map{ $_, 1 } @_;
    @_ = keys %hash;
    return @_;
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

sub getscore
{
    my $score=0;
    my $first=shift;
    my $second=shift;
     #$first=$map{$first};
     #$second=$map{$second};
    my @firstchar=split("",$first);
    #print Dumper(@firstchar);
    #print "\n";
    my @secondchar=split("",$second);
    #print @secondchar;
    #print "\n";
    for (my $l=0;$l<16;$l++)
    {
    my $one=index($library,$firstchar[$l]);
    my $two=index($library,$secondchar[$l]);
        $score=$score + $matrix[$one][$two];    
    }
    #print $score;
    #print "\n";
    return $score;
}


sub givescore
{
    my $score=0;
    my $first=shift;
    my $second=shift;
     #$first=$map{$first};
     #$second=$map{$second};
    my @firstchar=split("",$first);
    #print @firstchar;
    #print "\n";
    my @secondchar=split("",$second);
    #print @secondchar;
    #print "\n";
    for (my $l=0;$l<8;$l++)
    {
    my $one=index($library,$firstchar[$l]);
    my $two=index($library,$secondchar[$l+4]);
        $score=$score + $matrix[$one][$two];    
    }
    #print $score;
    #print "\n";
    return $score;
}

sub present
{
    my $score=0;
    my $first=shift;
    my $second=shift;
     #$first=$map{$first};
     #$second=$map{$second};
    my @firstchar=split("",$first);
    #print @firstchar;
    #print "\n";
    my @secondchar=split("",$second);
    #print @secondchar;
    #print "\n";
    for (my $l=0;$l<8;$l++)
    {
    $firstchar[$l+4]=$secondchar[$l];
    }
    #print @firstchar;
    #print "\n";
    $score=join("",@firstchar);
    #print "$score";
    #print "\n";
    return $score;
}




sub showHelp
{
	die("This tool calculates the average recognized self and non self peptides. Usage ./CheckEvol.pl [options]\n
	Options\n:
	-dir (directory where the *Genes.txt files are to be found)
	-out (outfile, optional, default is Collected_L=<directory name>.txt)
	-loci_kir
	-number_loci_mhc
	\n");
}
