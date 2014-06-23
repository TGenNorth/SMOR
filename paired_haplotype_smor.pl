#!/usr/bin/env perl

use strict;
use warnings;
#use Data::Dumper;

my $positionstocheck =
    [ [ "eisPlus1000", 1319 ],
    [ "eisPlus1000", 1272 ],
    [ "eisPlus1000", 1270 ],
    [ "eisPlus1000", 1268 ],
    [ "gyrAPlus1000", 1059 ],
    [ "gyrAPlus1000", 1062 ],
    [ "gyrAPlus1000", 1101 ],
    [ "gyrAPlus1000", 1108 ],
    [ "gyrAPlus1000", 1110 ],
    [ "gyrAPlus1000", 1119 ],
    [ "gyrAPlus1000", 1120 ],
    [ "gyrAPlus1000", 1123 ],
    [ "gyrAPlus1000", 1175 ],
    [ "gyrAPlus1000", 1199 ],
    [ "gyrAPlus1000", 1214 ],
    [ "inhAPlus1000", 1099 ],
    [ "inhAPlus1000", 1124 ],
    [ "inhAPlus1000", 1141 ],
    [ "inhAPlus1000", 1143 ],
    [ "inhAPlus1000", 1150 ],
    [ "katGPlus1000", 1094 ],
    [ "katGPlus1000", 1095 ],
    [ "rpoBPlus1000", 1077 ],
    [ "rpoBPlus1000", 1084 ],
    [ "rpoBPlus1000", 1087 ],
    [ "rpoBPlus1000", 1089 ],
    [ "rpoBPlus1000", 1097 ],
    [ "rpoBPlus1000", 1098 ],
    [ "rpoBPlus1000", 1099 ],
    [ "rpoBPlus1000", 1117 ],
    [ "rpoBPlus1000", 1128 ],
    [ "rpoBPlus1000", 1129 ],
    [ "rpoBPlus1000", 1144 ],
    [ "rpoBPlus1000", 1150 ],
    [ "rrsPlus1000", 1063 ],
    [ "rrsPlus1000", 1124 ],
    [ "rrsPlus1000", 1146 ],
    [ "rrsPlus1000", 1148 ] ];

if( @ARGV != 1 )
{
  print <<EOF;
Usage:
paired_haplotype_smor.pl <inputfile.bam>
EOF
  exit();
}

my $inputfilename = shift();
print "Chromosome\tPos\t|\t#AA\t#CC\t#GG\t#TT\t#Hom\t|\t%AA\t%CC\t%GG\t%TT\t|\t#A\t#C\t#G\t#T\t#Cov\t|\t%A\t%C\t%G\t%T\t|\tRaw\t\n";
foreach my $positiontocheck (@{$positionstocheck})
{
  my $chromosomename = $positiontocheck->[0];
  $positiontocheck = $positiontocheck->[1];
  my $readdata = {};
  if( open( my $samfiledata, '-|', "samtools view '$inputfilename'"  ) )
  {

    my $linefromfile = '';
    while( $linefromfile = <$samfiledata> )
    {
      my @linedata = split( "\t", $linefromfile );
      my $readname = $linedata[0];
      my $readchromosome = $linedata[2];
      my $startposition = $linedata[3];
      my $manglingtype = $linedata[5];
      my $wholeread = $linedata[9];
      if( $readchromosome eq $chromosomename )
      {
        # This is an oversimplification of a CIGAR parser
        if( $manglingtype =~ /^((?:\d{1,8}S)?)(\d{1,8})M/ )
        {
          my $cutfromfront = $1;
          my $lengthofread = $2;
          if( $cutfromfront =~ /^(\d{1,8})S$/ ){ $wholeread = substr( $wholeread, $1 ); }
          $wholeread = substr( $wholeread, 0, $lengthofread );
          if( ( $startposition <= $positiontocheck ) && ( ( $startposition + $lengthofread - 1 ) >= $positiontocheck ) )
          {
            if( !defined( $readdata->{$readname} ) ){ $readdata->{$readname} = []; }
            push( @{$readdata->{$readname}}, uc( substr( $wholeread, ( $positiontocheck - $startposition ), 1 ) ) );
          }
        }
      }
    }
    close( $samfiledata );

    my $basedistributions = {};
    my $definedcounts = 
    {
      'AA' => 0,
      'CC' => 0,
      'GG' => 0,
      'TT' => 0,
      'hom' => 0,
      'A' => 0,
      'C' => 0,
      'G' => 0,
      'T' => 0,
      'cov' => 0
    };
    foreach my $currentread ( keys( %{$readdata} ) )
    {
      foreach my $currentbase ( @{$readdata->{$currentread}} )
      {
        $definedcounts->{$currentbase}++;
        $definedcounts->{'cov'}++;
      }
      my $basedistribution = join( '', sort( @{$readdata->{$currentread}} ) );
      if( length( $basedistribution ) == 1 ){ $basedistribution .= '-'; }
      if( defined( $definedcounts->{$basedistribution} ) )
      {
        $definedcounts->{$basedistribution}++;
        $definedcounts->{'hom'}++;
      }
      if( !defined( $basedistributions->{$basedistribution} ) ){ $basedistributions->{$basedistribution} = 0; }
      $basedistributions->{$basedistribution}++;
    }
    print "$chromosomename\t$positiontocheck\t|\t";
    print "$definedcounts->{'AA'}\t$definedcounts->{'CC'}\t$definedcounts->{'GG'}\t$definedcounts->{'TT'}\t$definedcounts->{'hom'}\t|\t";
    if( $definedcounts->{'hom'} > 0 )
    {
      print sprintf( "%.4e", ( $definedcounts->{'AA'} * 100 / $definedcounts->{'hom'} ) ) . "%\t" . sprintf( "%.4e", ( $definedcounts->{'CC'} * 100 / $definedcounts->{'hom'} ) ) . "%\t" . sprintf( "%.4e", ( $definedcounts->{'GG'} * 100 / $definedcounts->{'hom'} ) ) . "%\t" . sprintf( "%.4e", ( $definedcounts->{'TT'} * 100 / $definedcounts->{'hom'} ) ) . "%\t|\t";
    } else { print "--\t--\t--\t--\t|\t"; }
    print "$definedcounts->{'A'}\t$definedcounts->{'C'}\t$definedcounts->{'G'}\t$definedcounts->{'T'}\t$definedcounts->{'cov'}\t|\t";
    if( $definedcounts->{'cov'} > 0 )
    {
      print sprintf( "%.4e", ( $definedcounts->{'A'} * 100 / $definedcounts->{'cov'} ) ) . "%\t" . sprintf( "%.4e", ( $definedcounts->{'C'} * 100 / $definedcounts->{'cov'} ) ) . "%\t" . sprintf( "%.4e", ( $definedcounts->{'G'} * 100 / $definedcounts->{'cov'} ) ) . "%\t" . sprintf( "%.4e", ( $definedcounts->{'T'} * 100 / $definedcounts->{'cov'} ) ) . "%\t|\t";
    } else { print "--\t--\t--\t--\t|\t"; }
    foreach my $basedistribution ( sort keys( %$basedistributions ) ){ print "$basedistribution: " . $basedistributions->{$basedistribution} . ", "; }
    print "\t\n";

  } else { print STDERR "Error reading $inputfilename!\n"; }
}


