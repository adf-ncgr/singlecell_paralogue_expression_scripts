#!/usr/bin/env perl
use strict;
use Getopt::Long;
use Text::NSP::Measures::2D::Fisher::right;
use Text::NSP::Measures::2D::Fisher::left;

my $feature_pair_file;
my %matrix_dir2cluster_file;
my $doublets;
my $keepers_enum;
my %cluster2mtx_barcode_suffixes;
#if set off, we will ignore bc that don't have clusters rather than considering it a fatal error
my $all_clustered=1;
my $summarize=1;
my $min_count_sum=0;
my $cluster_map;
my %cluster_map;
my %cluster_counts;
my $fisher=1;
my $wilson=0;
my $bias_margin=0;
GetOptions(
    "feature_pair_file=s" => \$feature_pair_file,
    "matrix_dir2cluster_file=s" => \%matrix_dir2cluster_file,
    "doublets=s" => \$doublets,
    "keepers_enum=s" => \$keepers_enum,
    "cluster_map=s" => \$cluster_map,
    "cluster2mtx_suffix=s" => \%cluster2mtx_barcode_suffixes,
    "all_clustered!" => \$all_clustered,
    "summarize!" => \$summarize,
    "min_count_sum=i" => \$min_count_sum,
    "fisher!" => \$fisher,
    "wilson!" => \$wilson,
    "bias_margin=f" => \$bias_margin,
);
die "must supply --feature_pair_file\n" unless defined $feature_pair_file;
die "must supply at least one pair to --matrix_dir2cluster_file\n" unless keys %matrix_dir2cluster_file;
if (defined $cluster_map) {
    open(CM, $cluster_map) || die $!;
    while (<CM>) {
        chomp;
        my ($assigned_cluster, $mapped_cluster) = split /\t/;
        $cluster_map{$assigned_cluster} = $mapped_cluster;
    }
    close CM;
}
open(FPF, $feature_pair_file) || die $!;
my %feature_pairs;
while (<FPF>) {
    chomp;
    my ($f1,$f2) = split /\t/;
    $feature_pairs{$f1}->{rank} = 1;
    $feature_pairs{$f1}->{paired_with} = $f2;
    #HACK! allow me to use this script for single genes by letting them pair with themselves
    if ($f2 cmp $f1) {
        $feature_pairs{$f2}->{rank} = 2;
        $feature_pairs{$f2}->{paired_with} = $f1;
    }
}
close FPF;

my %bc2cluster;
my %bcidx2cluster;
foreach my $matrix_dir (keys %matrix_dir2cluster_file) {
    my $clusters = $matrix_dir2cluster_file{$matrix_dir};
    open(FT, "zcat $matrix_dir/features.tsv.gz |") || die $!;
    my %ftidx2id;
    my $ftidx = 1;
    while (<FT>) {
        chomp;
        my ($ftid, $ftname, undef) = split /\s+/;
        $ftidx2id{$ftidx} = $ftid;
        $ftidx++;
    }
    close FT;

    my %cluster_ids;
    open(C, $clusters) || die;
    <C>;
    while (<C>) {
            chomp;
            my ($bc, $cluster) = split /,/;
            $bc = $matrix_dir.":".$bc;
            if (defined $cluster_map{$cluster}) {
                $cluster = $cluster_map{$cluster}
            }
            my ($cluster_barcode_suffix) = ($bc =~ /[ATGC]+(\S*)/);
            my $mtx_barcode_suffix = $cluster2mtx_barcode_suffixes{$cluster_barcode_suffix};
            if (keys %cluster2mtx_barcode_suffixes) { 
                #die "no mtx suffix corresponds to $cluster_barcode_suffix\n" unless defined $mtx_barcode_suffix;
                #next unless defined $mtx_barcode_suffix;
                $bc =~ s/$cluster_barcode_suffix$/$mtx_barcode_suffix/;
            }
            $bc2cluster{$bc} = $cluster;
            $cluster_counts{$cluster}++;
            $cluster_ids{$cluster} = 1;
    }
    close C;
    my %doublets;
    if (defined $doublets) {
        open(D, $doublets) || die $!;
        while (<D>) {
            chomp;
            $doublets{$_} = 1;
        }
        close D;
    }
    my %keepers;
    if (defined $keepers_enum) {
        open(K, $keepers_enum) || die $!;
        while (<K>) {
            chomp;
            $keepers{$_} = 1;
        }
    }

    open(BC, "zcat $matrix_dir/barcodes.tsv.gz |") || die $!;
    my $bcidx = 1;
    while (<BC>) {
        if ($doublets{$bcidx}) {
            $bcidx++;
            next;
        }
        chomp;
        my $bc = $matrix_dir.":".$_;
        if (defined $keepers_enum && ! $keepers{$bc}) {
            $bcidx++;
            next;
        }
        my $cluster = $bc2cluster{$bc};
        if ($all_clustered) {
            die "could not determine cluster for $_\n" unless defined $cluster;
        }
        $bcidx2cluster{$matrix_dir.":".$bcidx} = $cluster;
        $bcidx++;
    }
    close BC;


    open(MTX, "zcat $matrix_dir/matrix.mtx.gz |") || die $!;
    <MTX>;
    <MTX>;
    <MTX>;
    while (<MTX>) {
        chomp;
        my ($ftidx, $bcidx, $count) = split /\s+/;
        $bcidx = $matrix_dir.":".$bcidx;
        next if $doublets{$bcidx};
        my $ftid = $ftidx2id{$ftidx};
        if ($feature_pairs{$ftid}) {
            $feature_pairs{$ftid}->{cells}->{$bcidx} = $count;
        }
    }
    close MTX;
}

$, = "\t";
$\ = "\n";
if ($summarize) {
    if ($fisher) {
        print "overlap_pval", "nonoverlap_pval", "gene1", "gene2", "cluster", "both_on", "gene1_only", "gene2_only", "neither", "underpowered";
    }
    else {
        print "gene1", "gene2", "cluster", "both_on", "gene1_only", "gene2_only", "neither", "underpowered";
    }
}
else {
    print "gene1", "gene2", "cluster", "gene1+gene2", "gene1", "gene2";
}
foreach my $f (sort keys %feature_pairs) {
    next unless $feature_pairs{$f}->{rank} == 1;
    my $f2 = $feature_pairs{$f}->{paired_with};
    my %bc = map {$_ => 1} keys %{$feature_pairs{$f}->{cells}}, keys %{$feature_pairs{$f2}->{cells}};
    my %cluster_summaries;
    foreach my $bc (keys %bc) {
        #next if $doublets{$bcidx};
        next unless defined $bcidx2cluster{$bc};
        my $count1 = $feature_pairs{$f}->{cells}->{$bc};
        $count1 = 0 unless defined $count1;
        my $count2 = $feature_pairs{$f2}->{cells}->{$bc};
        $count2 = 0 unless defined $count2;
        my $cluster = $bcidx2cluster{$bc};
        my $count_sum = $count1+$count2;
        if ($count_sum >= $min_count_sum) {
            my $bias;
            if ($wilson) {
                $bias = &is_biased($count1, $count2);
                print STDERR $bc,$cluster,$count1,$count2,(defined $bias ? $bias : "NA");
            }
            if (! $summarize) {
                print $f, $f2, $cluster, $count1+$count2, $count1, $count2;
            }
            else {
                if (! $wilson) {
                    if ($count1 && $count2) {
                        $cluster_summaries{$bcidx2cluster{$bc}}->{count_both_on}++;
                    }
                    elsif ($count1) {
                        $cluster_summaries{$bcidx2cluster{$bc}}->{count_1_on}++;
                    }
                    elsif ($count2) {
                        $cluster_summaries{$bcidx2cluster{$bc}}->{count_2_on}++;
                    }
                }
                else {
                    if (! defined $bias) {
                        $cluster_summaries{$bcidx2cluster{$bc}}->{count_underpowered}++;
                    }
                    elsif ($bias == 0) {
                        $cluster_summaries{$bcidx2cluster{$bc}}->{count_both_on}++;
                    }
                    elsif ($bias == 1) {
                        $cluster_summaries{$bcidx2cluster{$bc}}->{count_1_on}++;
                    }
                    elsif ($bias == -1) {
                        $cluster_summaries{$bcidx2cluster{$bc}}->{count_2_on}++;
                    }
                }
            }
        }
    }
    if ($summarize) {
        foreach my $cluster (sort keys %cluster_summaries) {
            my $count_both_on = $cluster_summaries{$cluster}->{count_both_on};
            $count_both_on = 0 unless defined $count_both_on;
            my $count_1_on = $cluster_summaries{$cluster}->{count_1_on};
            $count_1_on = 0 unless defined $count_1_on;
            my $count_2_on = $cluster_summaries{$cluster}->{count_2_on};
            $count_2_on = 0 unless defined $count_2_on;
            my $count_underpowered = $cluster_summaries{$cluster}->{count_underpowered};
            $count_underpowered = 0 unless defined $count_underpowered;
            my $count_neither_on = $cluster_counts{$cluster} - ($count_both_on + $count_1_on + $count_2_on) - $count_underpowered;
            if ($fisher) {
                my $overlap_pval = Text::NSP::Measures::2D::Fisher::right::calculateStatistic(n11=>$count_both_on, n1p=>($count_both_on+$count_1_on), np1=>($count_2_on+$count_both_on), npp=>($count_both_on+$count_1_on+$count_2_on+$count_neither_on));
                my $nonoverlap_pval = Text::NSP::Measures::2D::Fisher::left::calculateStatistic(n11=>$count_both_on, n1p=>($count_both_on+$count_1_on), np1=>($count_2_on+$count_both_on), npp=>($count_both_on+$count_1_on+$count_2_on+$count_neither_on));
                print $overlap_pval, $nonoverlap_pval, $f, $f2, $cluster, $count_both_on, $count_1_on, $count_2_on, $count_neither_on, $count_underpowered;
            }
            else {
                print $f, $f2, $cluster, $count_both_on, $count_1_on, $count_2_on, $count_neither_on, $count_underpowered;
            }
        }
    }
}

#code below is copy-paste from jjd-homeologs.pl + modification to test for unbiased vs underpowered
sub is_biased() {
    my ($count1, $count2) = @_;
    my $sum = $count1+$count2;
    if ($sum < $min_count_sum || $sum == 0) {
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

sub wilson_test {
    my ($ns, $n, $z) = @_;
    my $nf = $n - $ns;
    my $center = ($ns + ($z*$z/2))/($n+($z*$z));
    my $dev = ($z/($n+($z*$z)))*sqrt((($ns*$nf)/$n)+(($z*$z)/4));
    return ($center, $dev);
}
