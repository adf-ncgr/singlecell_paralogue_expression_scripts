#!/usr/bin/env perl
use strict;
use Getopt::Long;
my $bias_margin = 0;
my $min_cell_counts = 0;
my $min_counts = 0;
my $sum_column = "count_sum";
my $cluster_matrix_output;
my $cell_counts_file;
my $subgenomes;
my $unbiased_value="NA";
my $underpowered_value="NA";
my $respect_given_pair_order=0;
GetOptions(
    "bias_margin=f" => \$bias_margin,
    "min_cell_counts=i" => \$min_cell_counts,
    "min_counts=f" => \$min_counts,
    "sum_column=s" => \$sum_column,
    "cluster_matrix_output=s" => \$cluster_matrix_output,
    "subgenomes=s" => \$subgenomes,
    "cell_counts_file=s" => \$cell_counts_file,
    "unbiased_value:s" => \$unbiased_value,
    "underpowered_value:s" => \$underpowered_value,
    "respect_given_pair_order" => \$respect_given_pair_order,
);
my $counts_file = shift;
my $homeolog_pairs_file = shift;

my %cell_counts;
if (defined $cell_counts_file) {
    open(CCF, $cell_counts_file) || die $!;
    while (<CCF>) {
	chomp;
	my ($gene,$remainder) = /^([^\t]+)\t(.*)/;
	my @counts = split /\t/, $remainder;
	$cell_counts{$gene} = \@counts;
    }
    close CCF;
}


open(CF, $counts_file) || die $!;
my $header=<CF>;
chomp $header;
my ($gene_header, @headers)=split /\t/, $header;
my $sum_column_idx;
if (defined $sum_column) {
    for (my $i=0; $i<@headers; $i++) {
        if ($headers[$i] eq $sum_column) {
            $sum_column_idx = $i;
            last;
        }
    }
}
my %gene_counts;
my @cluster_total_counts;
my @clusters_with_counts;
my @cluster_similarity_counts;
my @no_data;
while (<CF>) {
    chomp;
    my ($gene,$remainder) = /^([^\t]+)\t(.*)/;
    my @counts = split /\t/, $remainder;
    #NB: depends on these files being similarly ordered!
    if ($min_cell_counts) {
        my $cell_counts = $cell_counts{$gene};
        for (my $i=0; $i < @$cell_counts; $i++) {
            if ($cell_counts->[$i] < $min_cell_counts) {
                $counts[$i] = 0;
            }
        }
    }
    $gene_counts{$gene} = \@counts;
    for (my $i = 0; $i < @counts; $i++) {
        $cluster_total_counts[$i] += $counts[$i];
    }
    if (! @no_data) {
        my $count = split /\t/, $remainder;
        @no_data=split(//,"0"x$count);
    }
}
close CF;

my %subgenome_assignments;
if (defined $subgenomes) {
    open(S, $subgenomes) || die $!;
    while (<S>) {
        chomp;
        my ($g1,$s1,$g2,$s2) = split /\t/;
        $subgenome_assignments{$g1} = $s1;
        $subgenome_assignments{$g2} = $s2;
    }
    close S;
}
#print $header,"\tbiased1\tbiased2\tunbiased\tbiased1*biased2\n";
print $header,"\tbiased_count\tbiased_1\tbiased_2\tF_ex\tB_fix";
if (defined $subgenomes) {
    print "\tsubgenomes";
}
print "\n";
open(HP, $homeolog_pairs_file) || die $!;
my $pair = 1;
while (<HP>) {
    chomp;
    my ($g1,$g2) = split /\s+/;
    my $gc1 = $gene_counts{$g1};
    $gc1 = \@no_data unless $gc1;
    my $gc2 = $gene_counts{$g2};
    $gc2 = \@no_data unless $gc2;
    my @ordered_pair;
    #my $name_pair = "$pair:";
    my $name_pair = "";
    my $given_name_pair = "$g1/$g2";
    my $subgenome_pair;
    my $bias_direction = &compare($gc1,$gc2);
    if ($bias_direction >= 0) {
        @ordered_pair = ($gc1,$gc2);
        $name_pair .= "$g1/$g2";
        my $s1 = $subgenome_assignments{$g1};
        $s1 = "?" unless defined $s1;
        my $s2 = $subgenome_assignments{$g2};
        $s2 = "?" unless defined $s2;
        if ($bias_direction > 0) {
            $subgenome_pair = $s1.">".$s2;
        }
        else {
            $subgenome_pair = $s1."=".$s2;
        }
    }
    else {
        @ordered_pair = ($gc2,$gc1);
        $name_pair .= "$g2/$g1";
        my $s1 = $subgenome_assignments{$g1};
        $s1 = "?" unless defined $s1;
        my $s2 = $subgenome_assignments{$g2};
        $s2 = "?" unless defined $s2;
        $subgenome_pair = $s2.">".$s1;
    }
    #print $pair,"_",$g1,"\t",join("\t",@$gc1),"\n";
    #print $pair,"_",$g2,"\t",join("\t",@$gc2),"\n";
    #print ($respect_given_pair_order ? $given_name_pair : $name_pair),"\t";
    if ($respect_given_pair_order) {
        print $given_name_pair;
    }
    else {
        print $name_pair;
    }
    my $num_biased1 = 0;
    my $num_biased2 = 0;
    my $num_unbiased = 0;
    my $num_underpowered = 0;
    my @biases;
    for (my $i=0; $i < @$gc1; $i++) {
        my $sum = $gc1->[$i] + $gc2->[$i];
        #NA causes trouble for heatmaps in R, but only in terms of not getting dendograms. for now, we'll live without them
        my $ratio = "NA";
        $biases[$i] = 0;
        if ($sum >= $min_counts && $sum > 0) {
            $ratio = $ordered_pair[0]->[$i] / $sum;
            #exclude column representing the sums from being included in index calculations
            if ((!defined $sum_column_idx) || $i != $sum_column_idx) {
                my $bias = &is_biased($ordered_pair[0]->[$i], $ordered_pair[1]->[$i]);
                if ($bias == 1) {
                    $num_biased1++;
                    $biases[$i] = 1;
                }
                elsif ($bias == -1) {
                    $num_biased2++;
                    $biases[$i] = -1;
                }
                else {
		    if (defined $bias && $bias == 0) {
			$num_unbiased++;
			#per Jeremy's request, as I am currently interpresting it
			if ( length($unbiased_value) ) {
				$ratio=$unbiased_value;
			}
		    }
		    else  {
			$num_underpowered++;
			if ( length($underpowered_value) ) {
				$ratio=$underpowered_value;
			}
		    }
                }
            }
        }
	if ($respect_given_pair_order && !($name_pair eq $given_name_pair)) {
		print "\t", ($ratio =~ /^NA$/i ? $ratio : sprintf("%0.2f",1-$ratio));
	}
	else {
		print "\t", ($ratio =~ /^NA$/i ? $ratio : sprintf("%0.2f",$ratio));
	}
    }
    #print "\t$num_biased1\t$num_biased2\t$num_unbiased\t",$num_biased1*$num_biased2;
    print "\t",$num_biased1+$num_biased2,"\t",$num_biased1,"\t",$num_biased2;
    #F_ex = expression Fixation index
    #NB: before we distinguished between unbiased and underpowered, this would have included the latter in the denominator
    #print "\t",($num_biased1+$num_biased2+$num_unbiased ? sprintf("%0.2f",($num_biased1+$num_biased2)/($num_biased1+$num_biased2+$num_unbiased)) : "NA");
    print "\t",($num_biased1+$num_biased2+$num_unbiased+$num_underpowered ? sprintf("%0.2f",($num_biased1+$num_biased2)/($num_biased1+$num_biased2+$num_unbiased+$num_underpowered)) : "NA");
    #B_ex = expression Balance index
    my $balance=0;
    if ($num_biased1+$num_biased2) {
        $balance=(2*$num_biased2)/($num_biased1+$num_biased2);
    }
    print "\t",($num_biased1+$num_biased2 ? sprintf("%0.2f",$balance) : "NA");
    if (defined $subgenomes) {
        print "\t",$subgenome_pair;
    }
    print "\n";
    for (my $i=0; $i < @biases; $i++) {
        for (my $j=$i+1; $j < @biases; $j++) {
            my $similarity = $biases[$i]*$biases[$j];
            if ($similarity > 0) {
                #$cluster_similarity_counts[$i][$j]+=$balance;
                if ($balance) {
                $cluster_similarity_counts[$i][$j]+=1;
                }
                $clusters_with_counts[$i] = 1;
                $clusters_with_counts[$j] = 1;
            }
            elsif ($similarity < 0) {
                #$cluster_similarity_counts[$j][$i]+=$balance;
                if ($balance) {
                $cluster_similarity_counts[$j][$i]+=1;
                }
                $clusters_with_counts[$i] = 1;
                $clusters_with_counts[$j] = 1;
            }
        }
    }
    $pair++;
}
if ($cluster_matrix_output) {
    open(CS, ">$cluster_matrix_output") || die $!;
    print CS "cluster";
    for (my $i=0; $i < @headers; $i++) {
        next if defined $sum_column_idx && $i==$sum_column_idx;
        next unless $clusters_with_counts[$i];
        print CS "\t$headers[$i]";
    }
    print CS "\n";
    for (my $i=0; $i < @headers; $i++) {
        next if defined $sum_column_idx && $i==$sum_column_idx;
        next unless $clusters_with_counts[$i];
        print CS "$headers[$i]";
        for (my $j=0; $j < @headers; $j++) {
            next if defined $sum_column_idx && $j==$sum_column_idx;
            next unless $clusters_with_counts[$j];
            if ($i==$j) {
                print CS "\tNA";
            }
            else {
                print CS "\t$cluster_similarity_counts[$i][$j]";
            }
        }
        print CS "\n";
    }
}

sub is_biased() {
    my ($count1, $count2) = @_;
    my $sum = $count1+$count2;
    if ($sum < $min_counts || $sum == 0) {
        return 0;
    }
    #my $ratio = $count1 / $sum;
    #if ($ratio >= 1-$bias_margin) {
    my ($center, $dev) = &wilson_test($count1, $sum, 1.96);
    if ($center-$dev >= 1-$bias_margin) {
        return 1;
    }
    #elsif ($ratio <= $bias_margin) {
    elsif ($center+$dev <= $bias_margin) {
        return -1;
    }
    elsif ($center-$dev > $bias_margin && $center+$dev < 1-$bias_margin) {
        return 0;
    }
    return undef;
}

sub compare() {
    my ($counts1, $counts2) = @_;
    my $size = scalar(@$counts1);
    if (! $size) {
        $size = scalar(@$counts2);
    }
    my $compare=0;
    for (my $i=0; $i < $size; $i++) {
        next if defined $sum_column_idx && $i==$sum_column_idx;
        if (&is_biased($counts1->[$i], $counts2->[$i])) {
            if ($counts1->[$i] >= $counts2->[$i]) {
                $compare += 1;
            }
            else {
                $compare -= 1;
            }
        }
    }
    return $compare;
}

sub wilson_test {
    my ($ns, $n, $z) = @_;
    my $nf = $n - $ns;
    my $center = ($ns + ($z*$z/2))/($n+($z*$z));
    my $dev = ($z/($n+($z*$z)))*sqrt((($ns*$nf)/$n)+(($z*$z)/4));
    return ($center, $dev);
}
