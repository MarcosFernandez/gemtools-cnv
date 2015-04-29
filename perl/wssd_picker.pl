#!/usr/bin/perl -w

use strict;
use Getopt::Std;

my $dbug = 1;
my $unit = 1000;
my $hide = 0;
my %opts;

if ( !defined $ARGV[0]) 
  {
    print "usage: $0
         given a file with chrom, window_start, window_end, reads_count, reads_bases
         count how many windows in the -w contiguous windows are of value >= cutoff (-c)
         if over -s, signal as false negative duplication region
         output the chrom, start and end of such regions.

         NEW: Refine the boundary with a second file (-k) that has chrom, window_start, window_end, reads_count, reads_bases
              using the refined window size (-i), cut at the 1st/last new window of 1st/last boundary window with >=
              new_win_size/old_win_size * threshold (-c)


        -f input file name
        -w contiguous window number to look for high-value windows
        -s window number cutoff

        -c value cutoff
        -b number of column in -f file used to compare with cutoff given by -c, starting from 0
        -m switch. should pick windows whose [-b] value is < -c

        -k second file with -i as window size
        -n first round window size
           (given this number is to see how many windows to cut for 1st/last window for each wssd region)
        -i refined window size
        -t new threshold for refined window size
        -o output file\n\n";

    exit;
  }


getopts('mb:f:w:s:c:o:k:i:n:t:', \%opts);
die "please give contiguous window size and window number cutoff\n"
  if( !defined $opts{'w'} || $opts{'w'} !~ /^\d+$/ ||
      !defined $opts{'s'} || $opts{'s'} !~ /^\d+$/ ||
      !defined $opts{'b'} || $opts{'b'} !~ /^\d+$/ ||
      (defined $opts{'i'} && $opts{'i'} !~ /^\d+$/ ) ||
      !defined $opts{'n'} || $opts{'n'} !~ /^\d+$/  );


my $ctgw  = $opts{'w'};
my $highw = $opts{'s'};
my $newW  = $opts{'i'};
my $oriW  = $opts{'n'};
my $valCol= $opts{'b'};
#print join(", ", $ctgw, $highw, $newW, $oriW, $valCol, $opts{'m'}), "\n";

my %wins = ();
open(IN, $opts{'f'}) || die "cannot open $opts{'f'}: $!\n";
while( <IN> )
  {
    chomp;
    my @data = split(/\t/);
    next if( scalar(@data) < $valCol || $data[1] !~ /^\d+$/ );

    push @{ $wins{$data[0]} }, [$data[1], $data[2], $data[$valCol]]      if( defined  $wins{$data[0]} );
    $wins{$data[0]} =        [ [$data[1], $data[2], $data[$valCol]] ]    if( !defined $wins{$data[0]} );
  }
close(IN);

# Algo : start with every win, counting $ctgw windows forward, if >= $highW wins found having value > threshold,
#        store start of 1st win and end of last win with value > threshold as a WSSD region, signal flag,
#        when it goes into a region where <$highW wins of value > threshold in $ctgw countinuous wins
#
my %wssd_hash = ();
foreach my $chrom (sort keys %wins)
{

    my $winlist   = $wins{$chrom};
    my $wssdFlag  = 0;
    my $wssdChrom = '';
    my $wssdStart = 0;
    my $wssdLastStart = 0;                   # variable to store start of last win in WSSD region
    my $wssdEnd   = 0;                       # variable to store end of wssd region
    my $i         = 0;

   @$winlist = sort {$a->[0] <=> $b->[0]} @$winlist;

    while( $i <= scalar( @$winlist ) - $ctgw )
      {

	# for each contiguous window wide
	my $highCt    = 0;    # how many wins in $ctgw with value > threshold
	my $start     = -1;   # to store the start pos of first win with value > threshold
	my $lastStart = -1;   # to store the start pos of last win with value > threshold
	my $end       = -1;   # to store the end pos of last win with value > threshold
	
	for(my $j = $i; $j < $i + $ctgw; $j++)
	  {
	    if( ggcmp($winlist->[$j]->[2], $opts{'c'}, $opts{'m'}) )
	      {
		$highCt++;
		$start     = $winlist->[$j]->[0] if($start == -1);
		$lastStart = $winlist->[$j]->[0];
		$end       = $winlist->[$j]->[1];
	      }
	  }

	if( $highCt >= $highw )
	  {
	    if( !$wssdFlag )
	      {
		($wssdChrom, $wssdStart) = ($chrom, $start);
		$wssdFlag = 1;
	      }
	    $wssdLastStart = $lastStart;
	    $wssdEnd       = $end;
	  }

	# 1st window with windows of high values < highw after start WSSD
	# or last possible countinous $ctgw windows for the chrom, stop WSSD recording.
	# start of last WSSD win (value > threshold) is also stored for finding last refined windows later
	if( $wssdFlag && ( $i == scalar( @$winlist ) - $ctgw || $highCt < $highw ) )
	  {
	    push @{ $wssd_hash{$wssdChrom} }, [$wssdStart, $wssdEnd, $wssdLastStart] if( defined $wssd_hash{$wssdChrom}  );
	    $wssd_hash{$wssdChrom} =      [ [$wssdStart, $wssdEnd, $wssdLastStart] ] if( !defined $wssd_hash{$wssdChrom} );
	    $wssdFlag = 0;
	  }


	$i++;
      }

    # take into consideration when total window # = highw and all of them are high enough
    if( scalar( @$winlist ) == $highw && $highw < $ctgw )
      {
	my $takeWSSD = 1;
	($wssdStart, $wssdEnd) = ( $winlist->[0]->[0],  $winlist->[scalar( @$winlist ) -1]->[1] );
	for(my $j = 0; $j < scalar( @$winlist ); $j++)
	  {
	    if( ! ggcmp($winlist->[$j]->[2], $opts{'c'}, $opts{'m'}) ) {  $takeWSSD = 0; last;	  }
	  }
	if( $takeWSSD ) { $wssd_hash{$chrom} = [ [$wssdStart, $wssdEnd, $wssdStart] ]; }
      }
  }

# recycle the hash for refined window file if given
if( $opts{'k'} )
{
    %wins = ();
    open(IN, $opts{'k'}) || die "cannot open $opts{'k'}: $!\n";

    while( <IN> )
      {
	chomp;
	my @data = split(/\t/);
	next if( $data[1] !~ /^\d+$/ );

	push @{ $wins{$data[0]} }, [$data[1], $data[2], $data[$valCol]]   if(  defined $wins{$data[0]} );
	$wins{$data[0]} =        [ [$data[1], $data[2], $data[$valCol]] ] if( !defined $wins{$data[0]} );
      }
    close(IN);

    open(OUT, ">$opts{'o'}") || die "cannot open $opts{'o'} to write: $!\n";
    print OUT "chrom\twssd_start\twssd_end\n";

    foreach my $key (sort keys %wssd_hash)
      {
	my $refineWin = $wins{$key};
	my $wssdRoA   = $wssd_hash{$key};
	my $refineIdx = 0;
	my @wssdSort  = sort { $a->[0] <=> $b->[0] } @$wssdRoA;
	my @rWinSort  = sort { $a->[0] <=> $b->[0] } @$refineWin;
	my $wssd_prc  = 0;
	my @refWSSD   = ();

	
	foreach my $cdRoA (@wssdSort)
	  {
	    my $j = $refineIdx;
	    #print "@$cdRoA\n";
	    #my $pause = <STDIN>;
	    my @ovalWSSD = ();             # WSSD no matter can be refined or not

	
	    while( $j < scalar(@rWinSort) )
	      {
		
		if( $rWinSort[$j]->[1] <= $cdRoA->[0] )
		  {
		    $j++;
		    next;
		  }
		last if( $rWinSort[$j]->[0] >= $cdRoA->[1] );
		

		#print "wssd=@$cdRoA, refwin=@{$rWinSort[$j]}\n" if ($wssd_prc == 1);
		# the 1st refined window cut in current WSSD region, must have the same start as cur. WSSD
		# but window complete in gap may be removed (original window size larger, so start original
		# window may be partly gaps, but start refined win is completely in gap and removed
		# wssd_prc is used to mark the matching of refined windows in the 1st original window has 

		if( !$wssd_prc && $rWinSort[$j]->[0] >= $cdRoA->[0] )
		  {
		    $refineIdx = $j;           # where to start in refined win for next WSSD record
		    my $foundHigh = 0;

		    while( $j < scalar(@rWinSort) && $j < $refineIdx + $oriW/$newW  )
		      {
			if( $rWinSort[$j]->[2] > $opts{'t'} )
			  {
			    #print OUT $key, "\t", $rWinSort[$j]->[0], "\t";
			    push @ovalWSSD, $rWinSort[$j]->[0];
			    $foundHigh = 1;
			    last;
			  }
			$j++;
		      }

		    #print OUT $key, "\t", $cdRoA->[0], "\t" if( !$foundHigh );         # print the original boundary if no sgl high
		    push @ovalWSSD, $cdRoA->[0] if( !$foundHigh );
		    $wssd_prc = 1;
		    #print "foundHigh=$foundHigh\n";
		    #my $pause = <STDIN>;
		  }


		# the 1st refined window cut in the last WSSD orignal win,
		# (start of last WSSD original win eqals THIS refined win start
		if( $rWinSort[$j]->[0] >= $cdRoA->[2] )
		  {
		    my $curIdx      = $j;
		    my $lastHighEnd = -1;

		    # if count $oriW/$newW refined windows
		    # a window complete in gap may be removed so need to watch the end
		    while( $j < scalar(@rWinSort) && $rWinSort[$j]->[1] <= $cdRoA->[1] )
		      {
			if( ggcmp($rWinSort[$j]->[2], $opts{'t'}, $opts{'m'} ) )
			  {
			    $lastHighEnd = $rWinSort[$j]->[1];
			  }
			$j++;
		      }

		    #print OUT $lastHighEnd, "\n" if( $lastHighEnd > -1 );
		    #print OUT $cdRoA->[1],  "\n" if( $lastHighEnd == -1 );
		    push @ovalWSSD, $lastHighEnd   if( $lastHighEnd > -1 );
		    push @ovalWSSD, $cdRoA->[1]    if( $lastHighEnd == -1 );

		    $wssd_prc = 0;
		    #print "lastHighEnd=$lastHighEnd\n";
		    #my $pause = <STDIN>;
		    last;    # jump out of refine window loop as cur. WSSD is finished
		  }

		$j++;
	      } # loop of refine window

	    # it is possible no 1Kb window in the last 5Kb defining wssd is above 1Kb cutoff so no end is found
	    # in that case, drop this WSSD interval as really all 1Kb windows are below cutoff
	    push @refWSSD, [@ovalWSSD] if ( scalar(@ovalWSSD) == 2);
	  }# loop of raw WSSD
	

	@refWSSD = sort { $a->[0] <=> $b->[0] }  @refWSSD;
	foreach my $wssdRoA (@refWSSD)
	  {
	    print OUT "$key\t", $wssdRoA->[0], "\t", $wssdRoA->[1], "\n";
	  }

      }

    close(OUT);
    exit; # exit if given 2nd file
  }

open(OUT, ">$opts{'o'}") || die "cannot open $opts{'o'} to write: $!\n";
print OUT "chrom\twssd_start\twssd_end\n";
foreach my $key (sort keys %wssd_hash)
  {
    my $wssdRoA = $wssd_hash{$key};
    foreach my $cdRoA (@$wssdRoA)
      {
	print "'$key', '$cdRoA->[0]', '$cdRoA->[1]'\n" if( length($key) == 0 );
	print OUT $key, "\t", $cdRoA->[0], "\t", $cdRoA->[1], "\n";
      }
  }
close(OUT);

sub ggcmp
  {
    my ($val, $cut, $rev) = @_;

    #if($rev ) { print "m is set:     val=$val, cut=$cut\n"; }
    #if(!$rev) { print "m is not set: val=$val, cut=$cut\n"; }
    #my $pause = <STDIN>;

    return 1 if ($rev  && $val < $cut );
    return 1 if (!$rev && $val > $cut );

    return 0;
  }


