#!/usr/local/bin/perl -w

# TODO change minright to a positive number, call it extright.
# needed for multiplicative criterion.

use Getopt::Long;
use strict;

my @ARGV_COPY  = @ARGV;
my $n_args = @ARGV;

$::debug    =  0;
$::test     =  0;
my $help    =  0;
my $fh      =  \*STDIN;
my $annot   =  "";
my $tag     =  "tag";
my $progname = 'syloscope';
$::limit    =  40;
$::mode     =  '';
$::geometric = 0;
my $fnsylbase = "";
my $fnsylbaze = "";
my $fnsyl6  =  "sylresult-6";
my $fnsyl7  =  "sylresult-7";
my $fnsyl8  =  "sylresult-8";
my $fntable =  "syltable";
my $cut_entropy = 0.0;
my $fraction=  0;
my $range   =  0;

# --debug                 debug
# --test                  test

sub help {
   print <<EOH;
Usage:
   $progname [options]
Options:
--help                  this
EOH
}
#--over                  only consider overrepresentation
#--under                 only consider underrepresentation

if
(! GetOptions
   (  "help"            =>    \$help
   ,  "test"            =>    \$::test
   ,  "debug=i"         =>    \$::debug
   ,  "annot=s"         =>    \$annot
   ,  "limit=i"         =>    \$::limit
   ,  "syl6=s"          =>    \$fnsyl6
   ,  "syl7=s"          =>    \$fnsyl7
   ,  "syl8=s"          =>    \$fnsyl8
   ,  "syl=s"           =>    \$fnsylbase
   ,  "sylgz=s"         =>    \$fnsylbaze
   ,  "table=s"         =>    \$fntable
   ,  "fraction=f"      =>    \$fraction
   ,  "range=i"         =>    \$range
   ,  "entropy=f"       =>    \$cut_entropy
   ,  "mode=s"          =>    \$::mode
   ,  "geometric"       =>    \$::geometric
   ,  "tag=s"           =>    \$tag
   )
)
   {  print STDERR "option processing failed\n";
      exit(1);
   }

if ($help) {
   help();
   exit(0);
}

if ($fraction) {
}
elsif ($range) {
   $fraction = 1/$range;
}
else {
   die "don't know what range to take\n";
   $fraction = 1;
}

if ($fnsylbase) {
   ($fnsyl6, $fnsyl7, $fnsyl8) = map { "$fnsylbase.$_" } (6, 7, 8);
}
elsif ($fnsylbaze) {
   ($fnsyl6, $fnsyl7, $fnsyl8) = map { "$fnsylbaze.$_.gz" } (6, 7, 8);
}

open(FH6, "gzip -dcf $fnsyl6|") || die "can't open file $fnsyl6";
open(FH7, "gzip -dcf $fnsyl7|") || die "can't open file $fnsyl7";
open(FH8, "gzip -dcf $fnsyl8|") || die "can't open file $fnsyl8";
open(FHT, "<$fntable") || die "can't open file $fntable";

my @lead6 = map { chomp; split; } scalar <FH6>;
shift @lead6;
my @lead7 = map { chomp; split; } scalar <FH7>;
shift @lead7;
my @lead8 = map { chomp; split; } scalar <FH8>;
shift @lead8;

my %table = ();
my %sylresult = ();

my $leadt = <FHT>;
chomp $leadt;
my @lead = split "\t", $leadt;
my %colmap = ();
for (my $i=0;$i<@lead;$i++) {
   $colmap{$lead[$i]} = $i;
}
my $n_fields = keys %colmap;
my @missing = grep { !defined($colmap{$_}) } qw(handle six2 six1A six3 seven1A seven2 eight1A type);
die "missing fields @missing in table" if @missing;

for (<FHT>) {
   chomp;
   my @fields = split "\t", $_;
   die "field count error in line $. [$_]" unless scalar @fields == $n_fields;
   my $syl81 = $fields[$colmap{eight1A}];
   die "eightmer error line $. [$_]" unless defined($syl81);
   for (keys %colmap) {
      $table{$syl81}{$_} = $fields[$colmap{$_}];
   }
}
close FHT;

for (<FH6>) {
   chomp; my @F = split; my $sixmer = shift @F; $sylresult{$sixmer} = [@F];
} close FH6;
for (<FH7>) {
   chomp; my @F = split; my $sevmer = shift @F; $sylresult{$sevmer} = [@F];
} close FH7;
for (<FH8>) {
   chomp; my @F = split; my $eightmer = shift @F; $sylresult{$eightmer} = [@F];
} close FH8;

print "eightmer\ttype\tmaxleft\tminright\tareadiff\tlocation\tmaxabs\n";

sub max {
   return $_[0] > $_[1] ? $_[0] : $_[1];
}

my @weights = (0, 0, 0, 0, 0, 0);
if ($::mode eq 'SUM5') {
   @weights = (1, 1, 0, 1, 1, 1);
}
elsif ($::mode eq 'SUM4') {
   @weights = (1, 0, 0, 1, 1, 1);
}
elsif ($::mode eq 'MAX4') {
   @weights = (1, 0, 0, 1, 1, 1);
}


for my $eightmer (sort keys %table) {

   my @words = map { $table{$eightmer}{$_} } qw(six2 six1A six3 seven1A seven2 eight1A);

   my @unify  = (0) x @{$sylresult{$words[0]}};
   my @window = ();
   my $mitype = $table{$eightmer}{type};

die "$words[0]" unless defined($sylresult{$words[0]});


   for (my $i=0; $i<@{$sylresult{$words[0]}};$i++) {        # loop over sylamer offsets

      my $sum = 0;
      my $max = 0;
      my $min = 0;

      for (my $w=0;$w<6;$w++) {
         my $score = $weights[$w] * $sylresult{$words[$w]}[$i];
         $sum += $score;
         $max = $score if $score > $max;
         $min = $score if $score < $min;
      }
      
      if ($::mode =~ /^SUM/) {
         $unify[$i] = $sum;
      }
      elsif ($::mode =~ /^MAX/) {
         $unify[$i] = $max > -1 * $min ? $max : $min;
      }
   }

   my ($maxleft, $minright) = (-1000, 1000);
   my ($areal, $arear) = (0, 0);
   my ($minpos, $maxpos) = (0, 0);

   for (my $j=0; $j < @lead6 && $lead6[$j] <= $fraction * $lead6[-1]; $j++) {
      if ($unify[$j] > $maxleft) {
         $maxleft = $unify[$j];
         $maxpos = $lead6[$j];
      }
   }
      # the last column are zeros, we ignore them.
      # this means that the maximum can be negative.
   for (my $j=@unify-1; $j >= 0 && $lead6[$j] >= (1 - $fraction) * $lead6[-1]; $j--) {
      if ($unify[$j] < $minright) {
         $minright = $unify[$j];
         $minpos = $lead6[$j];
      }
   }

   my $location = int(100 * ($maxleft > -$minright ? $maxpos : $minpos) / $lead6[-1]);
         #  get_location($eightmer, \@lead6, $sumsmooth, $maxleft > -$minright);

   for (my $i=0;$i<@unify/2;$i++) {
      my $j = @unify - $i -1;
      $areal += $unify[$i];
      $arear += $unify[$j];
   }
   $areal *= -1 if $areal < 0;
   $arear *= -1 if $arear < 0;
   $areal =   1 if $areal == 0;
   $arear =   1 if $arear == 0;
   my $areadiff = log($areal/$arear)/log(2);

   my $maxabs = abs($maxleft);
   $maxabs = abs($minright) if $maxabs < abs($minright);

   printf
      "%s\t%s\t%g\t%6g\t%6g\t%d\t%6g\n"
      ,  $eightmer
      ,  $mitype
      ,  $maxleft
      ,  -$minright
      ,  $areadiff
      ,  $location
      ,  $maxabs
   ;
}





sub get_entropy {

   my $word = $_[0];
   my %bases = ();
   for my $b (split "", $word) {
      $bases{$b}++;
   }
   my $entropy = 0;
   for (values %bases) {
      my $p = $_/length($word);
      $entropy -= log($p) * $p / log(4);
   }
   return $entropy;
}


sub get_smooth {

   my ($unify) = @_;
   my @sumsmooth = (0) x @$unify;

   my @window = ($unify->[0]);
   $sumsmooth[0] = $unify->[0];

   for (my $i=1;$i<@$unify;$i++) {
      if ($i < 3) {
         push @window, @$unify[2*$i-1, 2*$i];
      }
      elsif ($i+3 <= @$unify) {
         push @window, $unify->[$i];
         shift @window;
      }
      else {
         shift @window;
         shift @window;
      }
      # print STDERR "@window\n";
      my $smooth = 0;
      $smooth += $_ for @window;
      $sumsmooth[$i] = $smooth / @window;
   }
   return \@sumsmooth;
}


sub get_location {
   my ($word, $lead, $smooth, $lookleft) = @_;
   @$smooth = map { -1 * $_; } reverse @$smooth if !$lookleft;
   die "location mismash\n" unless @$lead == @$smooth;
   my ($max, $ofs) = (0, 0);
   for (my $i=0;$i<@$smooth;$i++) {
      if ($smooth->[$i] > $max) {
         $ofs = $lead->[$i];
         $max = $smooth->[$i];
# print STDERR "$ofs $max\n" if $word eq 'CTTTGTAA';
      }
      elsif ($ofs + 2000 <  $lead->[$i]) {
         last;
      }
   }
   $ofs = $lead->[-1] - $ofs if !$lookleft;
   return int(100*($ofs/$lead->[-1]));
}


