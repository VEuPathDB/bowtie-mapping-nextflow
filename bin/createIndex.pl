#!/usr/bin/perl

use strict;
use Getopt::Long;

my($sampleName,$bowtieIndex,$isColorspace,$tool);

&GetOptions( 
    "isColorspace:s" => \$isColorspace,
    "bowtieIndex|x=s" => \$bowtieIndex,
    "sampleName|s=s" => \$sampleName,
            );
if ($isColorspace eq "true"){
    $tool = "bowtie-build";
}
else {
    $tool = "bowtie2-build";
}   

system("$tool $bowtieIndex $sampleName");
