#!/usr/bin/perl

=head
 pathwayTool.pl -- tools for pathway analysis

 author: Yi Zheng
 2017-01-30 finish
=cut

use strict;
use warnings;
use IO::File;
use Getopt::Std;
#use Statistics::Multtest;
use FindBin;

my %options;
getopts('t:x:p:', \%options);
pwy_enrich(\%options, \@ARGV); # enrichment analysis

# ==================
# kentnf: subroutine
# ==================

=head2 
 pwy_enrich: pathway enrichment analysis
=cut
sub pwy_enrich
{
	my ($options, $files) = @_;
	my $usage = qq'
USAGE: $0 -t idlink -x dblink input_gene_list pathway_file output

	-t link for gene id
	-x link for pathway id
    -p cutoff of pvalue (default 0.05)

* the input gene list should be changed gene in DE analysis
* the pathway file should be output of pathways tools

';
	print $usage and exit unless @$files == 3;
	my ($gene_list,$pathway_file, $output_file) = @$files;
	print $usage and exit unless -s $gene_list;
	print $usage and exit unless -s $pathway_file;
	print $usage and exit unless defined $output_file;
	print $usage and exit unless defined $$options{'t'};
	print $usage and exit unless defined $$options{'x'};
	my $idlink = $$options{'t'};
	my $dblink = $$options{'x'};
	my $pcutoff = 0.05;
	$pcutoff = $$options{'p'} if defined $$options{'p'} && $$options{'p'} > 0;

	# load gene list (changed gene) to hash
	my %changed_gene;
	my $fh1 = IO::File->new($gene_list) || die $!;
	while(<$fh1>)
	{
		chomp;
		$changed_gene{$_} = 1;
	}
	$fh1->close;

	# load pathway to hash
	# key: pwy_id
	# value: pwy_name
	#
	# key: pwy_id
	# value: gene1 \t gene2 \t ... \t geneN
	# 
	# check pwy id and name uniq at same time
	my %pwy_name;
	my %pwy_gene;
	my %all_pwy_gene;
	my %all_pwy_changed_gene;

	my $fh2 = IO::File->new($pathway_file) || die $!;
	while(<$fh2>)
	{
		chomp;
		next if $_ =~ m/^#/;
		#gene_ID        gene_description        pathway_ID      pathway_name
		my @a = split(/\t/, $_);
		die "Error in line $_\n" unless scalar @a == 4;
		my ($gid, $g_desc, $pid, $p_name) = @a;

		# check pwy id and pwy name
		if (defined $pwy_name{$pid} )
		{
			die "Error in pwy $pid\n" if $p_name ne $pwy_name{$pid};
		}
		else
		{
			$pwy_name{$pid} = $p_name;
		}

		if (defined $pwy_gene{$pid})
		{
			$pwy_gene{$pid}.= "\t".$gid;
		}
		else
		{
			$pwy_gene{$pid} = $gid;
		}

		$all_pwy_gene{$a[0]} = 1;

		if ( defined $changed_gene{$a[0]} )
		{
			$all_pwy_changed_gene{$a[0]} = 1;
		}

	}
	$fh2->close;

	my $N = scalar(keys %all_pwy_gene);		# N: gene in all pathways
	my $n = scalar(keys %all_pwy_changed_gene);	# n: changed gene in pathways

	# uniq the gene in each pwy, then get pvalue of changed pathways
	my %uniq_gene;
	my @uniq_gene;

	my %pwy_clsgene; # key pwy id; value: array of cluster gene in pwy
    my $output_num = 0;
	my %output_hash; # key: pvalue, value, array of content line

	foreach my $pid (sort keys %pwy_gene)
	{
		%uniq_gene = ();
		my @gene = split(/\t/, $pwy_gene{$pid});

		foreach my $gid (@gene) { $uniq_gene{$gid} = 1; }
		@uniq_gene = sort keys %uniq_gene;

		my $M = scalar(@uniq_gene);		# M: gene in particular pathways
		my $x = 0;				# x: changed gene in particular pathways
		foreach my $gid (@uniq_gene) {
			if ( defined $changed_gene{$gid} ) {
				$x++;
				push(@{$pwy_clsgene{$pid}}, $gid);
			}
		}
	
		my $p_name = $pwy_name{$pid};

		# compute pvalue
		#################################################################################
		#  input format
		#  hypergeometric(N, n, M, x);
		#  N: gene in all pathways
		#  n: changed gene in pathways
		#  M: gene in particular pathways
		#  x: changed gene in particular pathways
		#
		#  check this link for confirmation:
		#  http://www.geneprof.org/GeneProf/tools/hypergeometric.jsp
		################################################################################

		if ($x > 0) {
			# the order should be N M n x
			# my $pvalue = hypergeometric($N, $n, $M, $x);
			my $pvalue = hypergeometric($N, $M, $n, $x);
			my $background = "$M out of $N genes";
			my $cluster = "$x out of $n genes";
			my $genes = join(", ", @{$pwy_clsgene{$pid}});
			my $output_line = "$pid\t$p_name\t$cluster\t$background\t$pvalue\t$genes\n";
			push (@{$output_hash{$pvalue}}, $output_line);
			$output_num++;
		}
	}

    if ( $output_num > 0) { 
    	my $out1 = IO::File->new(">".$output_file) || die $!;
		print $out1 "Pathway ID\tPathway name\tCluster frequency\tGenome frequency\tP-value\tGenes annotated to the pathway\n";
		foreach my $p (sort {$a<=>$b} keys %output_hash) {
			foreach my $m (@{$output_hash{$p}}) {
				print $out1 $m;
			}
		}
		$out1->close;
	}

	# adjust p value to q value
	if ( $output_num > 10) {
		my $output_dat = $output_file.".qvalue.txt";

		my $R_CODE =<< "END";
library(qvalue)
data<-read.delim(file="$output_file", header=TRUE)
p<-data[,5]
qobj<-qvalue(p)
alldata<-cbind(data,qobj\$qvalue)
write.table(file="$output_dat", sep="\\t", alldata);
END

		#print $R_CODE; exit;
		open R,"|/usr/bin/R --vanilla --slave" or die $!;
		print R $R_CODE;
		close R;
	}

	# append Qvalue to output and apply p cutoff
    # load q value to hash
	my %pwy_qvalue;
	if (-s $output_file.".qvalue.txt") {
		my $ap0 = IO::File->new($output_file.".qvalue.txt") || die $!;
		<$ap0>;
		while(<$ap0>) {
			chomp;
			$_ =~ s/"//ig;
			my @a = split(/\t/, $_);
			$pwy_qvalue{$a[1]} = $a[7];
		}
		$ap0->close;
	}

	my $add_qv_content;

	if (scalar keys %pwy_qvalue > 0) {
		$add_qv_content = "Pathway ID\tPathway Name\tCulster Frequency\tGenome Frequency\tP-value\tadj P\tGenes annotated to the pathway\n";
	} else {
		$add_qv_content = "Pathway ID\tPathway Name\tCulster Frequency\tGenome Frequency\tP-value\tGenes annotated to the pathway\n";
	}
	my $ap1 = IO::File->new($output_file) || die $!;
	<$ap1>;
	while(<$ap1>) {
		chomp;
		my @a = split(/\t/, $_);
		if (defined $pwy_qvalue{$a[0]}) {
			my $q = parse_pvalue($pwy_qvalue{$a[0]});
			my $p_return = parse_pvalue($a[4]);
			next if ($a[4] > $pcutoff && $q > $pcutoff);
			$add_qv_content.="$a[0]\t$a[1]\t$a[2]\t$a[3]\t$p_return\t$q\t$a[5]\n";
		}
		else
		{
			my $p_return = parse_pvalue($a[4]);
			next if $a[4] > $pcutoff;
			$add_qv_content.="$a[0]\t$a[1]\t$a[2]\t$a[3]\t$p_return\t$a[5]\n";
		}
	}
	$ap1->close;

    $ap1 = IO::File->new(">".$output_file) || die $!;
	print $ap1 $add_qv_content;
	$ap1->close;
	# append id link and db link to output

    my $ap2 = IO::File->new(">>". $output_file) || die $!;
    print $ap2 "#idlink\t$idlink\n#dblink\t$dblink\n";
    $ap2->close;
	# better to parse the output file format
}

#####################
###  Subroutines  ###
#####################

sub hypergeometric {
    my $n = $_[0]; # N Total number of genes in all the pathways
    my $np = $_[1];# M Total number of genes in a particular pathway
    my $k = $_[2]; # n Total number of changed genes (in the input list) from all the pathways
    my $r = $_[3]; # x total number of changed genes (in the input list) from the particular pathway
    my $nq;
    my $top;
    
    $nq = $n - $np;

    my $log_n_choose_k = lNchooseK( $n, $k );

    $top = $k;
    if ( $np < $k ) {
        $top = $np;
    }

    my $lfoo = lNchooseK($np, $top) + lNchooseK($nq, $k-$top);
    my $sum = 0;

    for (my $i = $top; $i >= $r; $i-- ) {
        $sum = $sum + exp($lfoo - $log_n_choose_k);

        if ( $i > $r) {
            $lfoo = $lfoo + log($i / ($np-$i+1)) +  log( ($nq - $k + $i) / ($k-$i+1)  )  ;
        }
    }
    return $sum;
}

sub lNchooseK {
    my $n = $_[0];
    my $k = $_[1];
    my $answer = 0;

    if( $k > ($n-$k) ){
        $k = ($n-$k);
    }

    for(my $i=$n; $i>($n-$k); $i-- ) {
        $answer = $answer + log($i);
    }

    $answer = $answer - lFactorial($k);
    return $answer;
}

sub lFactorial {
    my $returnValue = 0;
    my $number = $_[0];
    for(my $i = 2; $i <= $number; $i++) {
        $returnValue = $returnValue + log($i);
    }
    return $returnValue;
}

sub parse_pvalue {
	my $pvalue = shift;

	my $pvalue_return;
	if ($pvalue =~ m/(\S+)e(\S+)/) {
		my ($p1, $p2) = ($1, $2);
		$p1 = sprintf("%.2f", $p1);
		$pvalue_return = $p1."e".$p2;
	}
	else
	{
		my @np = split(//, $pvalue);
		my $dig = '';
		my $f = 0;
		foreach my $n (@np) {
			$dig.= $n;
			$f++ if ($n ne '.' && $n > 0);
			last if $f > 2;
		}
		$pvalue_return = $dig;
	}

	return $pvalue_return;
}
