#!/usr/bin/perl

#use strict;
use Statistics::LineFit;

my @xValues = ();
my @yValues = ();

while(<>)
{
    my @fields = split(/,/);
    push @xValues, $fields[0];
    push @yValues, $fields[1];
}

#for (my $i=0; $i<100; $i++)
#{
#    push @xValues, $i;
#    push @yValues, 2*$i +3;
#}



$lineFit = Statistics::LineFit->new();
$lineFit->setData (\@xValues, \@yValues) or die "Invalid data";
($intercept, $slope) = $lineFit->coefficients();
#defined $intercept or die "Can't fit line if x values are all equal";
#$rSquared = $lineFit->rSquared();
#$meanSquaredError = $lineFit->meanSqError();
#$durbinWatson = $lineFit->durbinWatson();
#$sigma = $lineFit->sigma();
#($tStatIntercept, $tStatSlope) = $lineFit->tStatistics();
#@predictedYs = $lineFit->predictedYs();
#@residuals = $lineFit->residuals();
#($varianceIntercept, $varianceSlope) = $lineFit->varianceOfEstimates();


print "M: $slope\n";
print "C: $intercept\n";
