#!/usr/bin/perl

use strict;

use Getopt::Long;


my ($inputFile, $outputGFF, $outputTab, $outputConfig, $sample, $profileSetName);

&GetOptions(
    "inputFile=s"       => \$inputFile,
    "outputGFF=s"       => \$outputGFF,
    "outputTab=s"       => \$outputTab,
    "outputConfig=s"    => \$outputConfig,
    "sample=s"          => \$sample,
    "profileSetName=s"  => \$profileSetName,
    
    );

open(IN, $inputFile) or die "Cannot open $inputFile for reading: $!";

open(OUT, ">$outputTab") or die "Cannot open $outputTab for writing: $!";

open(GFF, ">$outputGFF") or die "Cannot open $outputGFF for writing: $!";

open(CONFIG, ">$outputConfig") or die "Cannot open $outputConfig for writing: $!";

print OUT  "sequence_source_id\tsegment_start\tsegment_end\tscore1\tscore2\tp_value\n";

my @columns = ('Name', 'File Name', 'Source Id Type', 'Input ProtocolAppNodes', 'Protocol', 'ProtocolParams', 'ProfileSet');
print CONFIG join("\t", @columns) . "\n";


while (my $line=<IN>) {
    if ($line=~/^#/) {
        next;
    }
    chomp($line);
    my @arr = split(/\t/, $line);
    # Storing Normalized Tag Count in score1 and Fold Change vs Control in score 2.
    print OUT "$arr[1]\t$arr[2]\t$arr[3]\t$arr[5]\t$arr[10]\t$arr[11]\n";

    print GFF join("\t", $arr[1], "HOMER", "enriched_region", $arr[2], $arr[3], $arr[5], $arr[4], ".", "ID=${sample}_$arr[0]; foldChangeVsInput=$arr[10]; pvalueVsInput=$arr[11]"), "\n";
}

print CONFIG $sample. "_peaks (ChIP-Seq)\t$outputTab\tsegment\t\tHOMER peak calls\tstyle|histone;fdr|0.001\t$profileSetName\n";

close IN;
close OUT;
close GFF;
close CONFIG;
