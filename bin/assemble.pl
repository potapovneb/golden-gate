#!/usr/bin/perl -w

################################################################################
# Application of a multiplexed ligase end-joining fidelity assay to optimize Golden Gate assembly
# Copyright (C) 2018 New England Biolabs, Inc.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

use strict;
use Getopt::Long;

### command-line options
my $o_np        =  1;
my $o_verbose   = "";
my $o_gap_min   = -1;
my $o_gap_max   = -1;
my $o_epi_max   = -1;
my $o_pro_max   = -1;
my $o_symmetric =  0;
my $o_discard   = "";

my @o_nobreak   = ();
my %o_nobreak   = ();
my @o_keep      = ();
my %o_keep      = ();

GetOptions(
    "np=f"      => \$o_np,
    "verbose"   => \$o_verbose,
    "gap-min=f" => \$o_gap_min,
    "gap-max=f" => \$o_gap_max,
    "pro-max=f" => \$o_pro_max,
    "epi-max=f" => \$o_epi_max,
    "nobreak=s" => \@o_nobreak,
    "keep=s"    => \@o_keep,
    "symmetric" => \$o_symmetric,
    "discard"   => \$o_discard,
    );

for my $x ( split(/,/,join(",",@o_nobreak)) ) { $o_nobreak{$x} = 1; }
for my $x ( split(/,/,join(",",@o_keep)) ) { $o_keep{$x} = 1; }

if( @ARGV == 0 )
{
    print "usage: $0 blast_result.csv subreads_ccs.fasta subreads_ccs.sts.csv inserts.fasta\n\n";
    print "options:\n";
    print " --np\t\tminimum number of passes (default '$o_np')\n";
    print " --verbose\t(default '$o_verbose')\n";
    print " --nobreak\t\tdo not break assemblies into smaller parts at junctions of specified lengths (default '@o_nobreak')\n";
    print " --gap-min\tbreak assembly at any junction where gap length below cutoff (default '$o_gap_min')\n";
    print " --gap-max\tbreak assembly at any junction where gap length above cutoff (default '$o_gap_max')\n";
    print " --pro-max\tbreak assembly at any junction where pro-fragment has unmapped bases (default '$o_pro_max')\n";
    print " --epi-max\tbreak assembly at any junction where epi-fragment has unmapped bases (default '$o_epi_max')\n";
    print " --symmetric\t(default '$o_symmetric')\n";
    print " --discard\tdiscard entire assembly if any error (default '$o_discard')\n";
    print "\n";
    exit;
}

### command-line arguments
my $blast_file   = shift @ARGV;
my $ccs_file     = shift @ARGV;
my $zmw_file     = shift @ARGV;
my $inserts_file = shift @ARGV;

### load consensus reads, zmw stats, assembly data
my $ccs     = load_references($ccs_file);
my $zmws    = load_zmw_data($zmw_file);
my $data    = load_assemblies($blast_file,$zmws);
my $inserts = load_references($inserts_file);

### ----------------------------------------------------------------------------
### define break points, filter PacBio reads, etc.
preprocess_mapped_inserts($data);

### ----------------------------------------------------------------------------
### Naming conventions
###   Inserts: A, B, C, D, ...
###   Overhangs: 1/1', 2/2', 3/3', ...
###
### 0             1              2        9
###   ----------     ----------      ...     ---------- 
###       A       1'     B       2'       9'     J      10' (== 0)

### comvert inserts names to numbers
my %MAP = ();
my @sorted = sort keys %$inserts;

for( my $i = 0; $i < @sorted; $i++ )
{
    $MAP{ $sorted[$i] } = $i + 1;
}

### ----------------------------------------------------------------------------
### blast mapping data
open(T,">","table_00.csv") || die "Can't write 'table_00.csv': $!";
print_mapped_inserts($data,\*T,$zmws);
close(T);

### ----------------------------------------------------------------------------
### insert frequencies
open(T,">","table_01.csv") || die "Can't write 'table_01.csv': $!";
count_insert_frequencies($data,\*T);
close(T);

### ----------------------------------------------------------------------------
### overhang frequencies
open(T,">","table_02.csv") || die "Can't write 'table_02.csv': $!";
count_overhang_frequencies($data,\*T,$inserts);
close(T);

### ----------------------------------------------------------------------------
### assembly frequencies
open(T,">","table_03.csv") || die "Can't write 'table_03.csv': $!";
my $assemblies = count_assembly_frequencies($data,\*T);
close(T);

### ----------------------------------------------------------------------------
### aggreagte by size
open(T,">","table_04.csv") || die "Can't write 'table_04.csv': $!";
aggregate_assemblies_by_size($assemblies,\*T);
close(T);

### ----------------------------------------------------------------------------
### extract A-to-J assemblies
open(T,">","table_05.csv") || die "Can't write 'table_05.csv': $!";
assemblies_A_to_J($assemblies,\*T);
close(T);

################################################################################
#                             SUBROUTINES                                      #
################################################################################

sub load_references {
    my ($file) = @_;
    
    my %refs = ();
    my $name = "";
    
    open(FA,$file) || die "Can't open '$file': $!";
    
    while( my $line = <FA> )
    {
	chomp($line);
	
	if( substr($line,0,1) eq ">" )
	{
	    $name = substr($line,1);
	}
	else
	{
	    $refs{$name} .= $line;
	}
    }
    
    close(FA);
    
    return \%refs;
}

sub load_zmw_data {
    my ($file) = @_;
    
    my %zmws = ();
    my @head = ();

    open(IN,$file) || die "Can't open '$file': $!";
    
    while( my $line = <IN> )
    {
	chomp($line);
	
	my @tokens = split(/,/,$line);
	
	if( @head == 0 )
	{
	    ### extract data fields
	    @head = @tokens;
	}
	else
	{
	    ### store data values
	    my %entry = ();
	    
	    for( my $i = 0; $i < @head; $i++ )
	    {
		$entry{$head[$i]} = $tokens[$i];
	    }

	    $zmws{$entry{"Movie"}}{$entry{"ZMW"}} = \%entry;
	}
    }

    close(IN);

    return \%zmws;
}

sub load_assemblies {
    my ($file,$zmws) = @_;

    my %data = ();
    my @head = ();

    open(CSV,$file) || die "Can't open '$file': $!";

    while( my $line = <CSV> )
    {
	chomp($line);
	
	my @tokens = split(/,/,$line);
	
	if( @head == 0 )
	{
	    @head = @tokens;
	}
	else
	{
	    my %entry = ();
	    
	    for( my $i = 0; $i < @tokens; $i++ )
	    {
		$entry{$head[$i]} = $tokens[$i];
	    }
	    
	    ### keep min and max coordinate for each fragment
	    $entry{"smin"} = $entry{"sstart"} < $entry{"send"} ? $entry{"sstart"} : $entry{"send"};
	    $entry{"smax"} = $entry{"sstart"} > $entry{"send"} ? $entry{"sstart"} : $entry{"send"};

	    ### blast maps "plus" sequences as follows:
	    ###        5' --------------> 3'
	    ### ==========|=============|=======
	    ###         sstart         send

	    ### make sure that blast output is consistent
	    if( $entry{"sstrand"} eq "plus" && $entry{"sstart"} > $entry{"send"} )
	    {
		print STDERR "[ERROR] Unexpected blast output\n\n";
		print STDERR "[ERROR] $line\n";
		exit;
	    }
	    
	    ### blast maps "minus" sequences as follows:
	    ###        3' <-------------- 5'
	    ### ==========|=============|=======
	    ###         send          sstart

	    ### make sure that blast output is consistent
	    if( $entry{"sstrand"} eq "minus" && $entry{"sstart"} < $entry{"send"} )
	    {
		print STDERR "[ERROR] Unexpected blast output\n\n";
		print STDERR "[ERROR] $line\n";
		exit;
	    }

	    my ($movie,$zmw,@tag) = split(/\//,$entry{"sseqid"});
	    
	    ### filter by number of passes
	    next if( $$zmws{$movie}{$zmw}{"NP"} < $o_np );
	    
	    ### store fragment info
	    push @{$data{$entry{"sseqid"}}}, \%entry;
	}
    }

    close(CSV);

    return \%data;
}

sub preprocess_mapped_inserts {
    my ($data) = @_;

    for my $read ( sortccs(keys %$data) )
    {
	### sort fragments by location
	my @fragments = sort { $$a{"smin"} <=> $$b{"smin"} } @{$$data{$read}};
	
	my $correct = 1;

	### calculte distance
	for( my $i = 0; $i < @fragments; $i++ )
	{
	    my $len = "NA";
	    my $pro = "NA";
	    my $epi = "NA";
	    my $brk = "NA";

	    if( $i < @fragments - 1 )
	    {
		my $f1 = $fragments[$i];
		my $f2 = $fragments[$i+1];

		### Schematics
		###
		###     smin         smax  smin         smax
		### =====|============|=====|============|=====
		###      ------------->     ------------->
		###           f1                  f2

		### ligation junction length
		$len = $$f2{"smin"} - $$f1{"smax"} - 1;

		### Schematics
		###
		###               smin              smax
		### ===============|=================|=======================
		###          - - - ------------------- - - - >
		###       clipped  |   mapped part   |  clipped
		###              qstart             qend
		###
		
		### number of bases clipped in f1
		$pro = ($$f1{"sstrand"} eq "plus") ? abs($$f1{"qlen"} - $$f1{"qend"}) : ($$f1{"qstart"} - 1);

		### number of bases clipped in f2
		$epi = ($$f2{"sstrand"} eq "plus") ? ($$f2{"qstart"} - 1) : abs($$f2{"qlen"} - $$f2{"qend"});

		### blast might produce partial (local) alignments and
		### unaligned bases might become part of the ligation junction
		my $gap = $len - $pro - $epi;

		### Schematics
		###                4-base               blunt
		###          ligation junction        ligation
		### ===|===========|====|===========|====**====|============|====|=============|=====
		###    ------------>    ------------>          <-------------    <--------------
		###         f1                f2                    f6                  f5
		###
		### Due to experimental setup, separate assembly products might
		### ligate together (blunt end ligation). Break points define
		### if and where asseblies get broken into smaller assemblies.
		###
		### Additionally, there might be some alignment artifacts. There
		### is a flag ($correct), which indicates whether alignment
		### problems are present for a given PacBio CCS read. Such
		### reads can be optionally discarded from further analysis.

		### define break points
		$brk = 0;

		### unexpected gap length -> flag as a break point
		if( %o_nobreak && ! exists $o_nobreak{$gap} ) { $brk = 1; }
		if( $o_gap_min != -1 && $gap < $o_gap_min   ) { $brk = 1; }
		if( $o_gap_max != -1 && $gap > $o_gap_max   ) { $brk = 1; }

		### flag entire read as problematic based on gap length
		if( %o_keep && ! exists $o_keep{$gap} )
		{
		    $correct = 0;
		}

		### * this is used to control partial fragments (by default,
		### fragments are expected to map end-to-end).
		### * flag entire read as problematic if $pro and $epi beyond expected
		if( $o_pro_max != -1 && $pro > $o_pro_max ) { $brk = 1; $correct = 0; }
		if( $o_epi_max != -1 && $epi > $o_epi_max ) { $brk = 1; $correct = 0; }
	    }

	    $fragments[$i]{"len"} = $len;
	    $fragments[$i]{"pro"} = $pro;
	    $fragments[$i]{"epi"} = $epi;
	    $fragments[$i]{"brk"} = $brk;
	}

	### discard problematic reads
	if( !$correct && $o_discard )
	{
	    delete $$data{$read};
	}
    }
}

sub print_mapped_inserts {
    my ($data,$fh,$zmws) = @_;
    
    for my $read ( sortccs(keys %$data) )
    {
	### sort fragments by location
	my @fragments = sort { $$a{"smin"} <=> $$b{"smin"} } @{$$data{$read}};
	
	### print data
	print_inserts(\@fragments,$fh,$zmws);
    }
}

sub print_inserts {
    my ($fragments,$fh,$zmws) = @_;

    ### print out PacBio CCS read name
    if( ! defined $zmws )
    {
	print $fh $$fragments[0]{"sseqid"}, "\n";
    }
    else
    {
	my ($movie,$zmw,@tag) = split(/\//,$$fragments[0]{"sseqid"});

	### print additional read info (optional)
	print $fh join(",",
		       $$fragments[0]{"sseqid"},
		       $$zmws{$movie}{$zmw}{"NP"},
		       $$zmws{$movie}{$zmw}{"ReadLength"} ), "\n";
    }

    ### print out all mapped fragments
    my @cols = qw/qseqid sstrand pident sstart send slen qstart qend qlen pro len epi brk/;

    print $fh join("\t", @cols), "\n";

    for my $info ( sort { $$a{"smin"} <=> $$b{"smin"} } @$fragments )
    {
	print $fh join( "\t", map { $$info{$_} } @cols), "\n";
    }
}

sub count_insert_frequencies {
    my ($data,$fh) = @_;
    
    my $total = 0;
    my %counts = ();

    for my $read ( sortccs(keys %$data) )
    {
	### sort fragments by location
	my @fragments = sort { $$a{"smin"} <=> $$b{"smin"} } @{$$data{$read}};
	
	for my $f ( @fragments )
	{
	    $total++;
	    $counts{$$f{"qseqid"}}++;
	}
    }

    print $fh join(",",qw/Insert Count Fraction/), "\n";
    
    my $total_fraction = 0;

    for my $insert ( sort keys %counts )
    {
	my $fraction = $counts{$insert} / $total;
	print $fh join(",", $insert, $counts{$insert}, $fraction ), "\n";

	### sanity check
	$total_fraction += $fraction;
    }

    print $fh join(",","Total",$total,$total_fraction);
}

sub count_overhang_frequencies {
    my ($data,$fh) = @_;

    ### collect cross-talk matrix
    my %matrix = ();
    my %list = ();

    ### analyze all pairs of inserts in assembly
    for my $read ( sortccs(keys %$data) )
    {
	### sort fragments by location
	my @fragments = sort { $$a{"smin"} <=> $$b{"smin"} } @{$$data{$read}};

	for( my $i = 0; $i < @fragments - 1; $i++ )
	{
	    my $f1 = $fragments[$i]{"qseqid"};
	    my $f2 = $fragments[$i+1]{"qseqid"};

	    ### insert numbers
	    my $n1 = $MAP{$f1};
	    my $n2 = $MAP{$f2};

	    ### skip break points
	    next if( $fragments[$i]{"brk"} == 1 );

	    ### pull out assembly sequence
	    my $assembly_sequence = $$ccs{$fragments[$i]{"sseqid"}};
	    
	    ### extract overhang sequence
	    my $oseq = substr($assembly_sequence,$fragments[$i]{"smax"},4);
	    $oseq = reverse(complement($oseq));

	    ### deduce idenity of the overhang for [$i] insert
	    if( $fragments[$i]{"sstrand"} eq "plus" )
	    {
	    	###  1  [Insert B == 2]  2
	    	### >>>>>>>>>>>>>>>>>>>>
	    	###     >>>>>>>>>>>>>>>>>>>>
	    	###  1'                  2'

		### overhang name
	    	$n1 = sprintf("%i'", $n1);
	    }
	    else
	    {
	    	###  2  [Insert B == 2]  1
	    	### <<<<<<<<<<<<<<<<<<<<
	    	###     <<<<<<<<<<<<<<<<<<<<
	    	###  2'                  1'

		### overhang name
	    	$n1 = sprintf("%i", $n1-1);
	    }
	    
	    ### overhang names for [$i+1] insert
	    if( $fragments[$i+1]{"sstrand"} eq "plus" )
	    {
	    	$n2 = sprintf("%i", $n2-1);
	    }
	    else
	    {	
	    	$n2 = sprintf("%i'", $n2);
	    }

	    ### rename overhangs in end junctions
	    if( $n1 eq "0" ) { $n1 = "10'"; }
	    if( $n2 eq "0" ) { $n2 = "10'"; }

	    ### cross-talk
	    $matrix{$n1}{$n2}++;

	    ### keep overhang sequences
	    $list{$n1}{$oseq}++;
	}
    }

    ### ordered list of overhang names
    my @names = ();

    for my $insert ( sort keys %MAP )
    {
	push @names, sprintf("%i", $MAP{$insert}), sprintf("%i'", $MAP{$insert});
    }

    ### representative overhang sequences
    my %legend = ();

    for my $o ( @names )
    {
	if( exists $list{$o} )
	{
	    my @sorted = sort {$list{$o}{$b} <=> $list{$o}{$a}} keys %{$list{$o}};

	    $legend{$o} = $sorted[0];
	}
	else
	{
	    $legend{$o} = "NA";
	}
    }

    ### print overhangs cross-talk
    print $fh join(",", "Overhang", map { '"'.$_.'"' } @names, "Sequence"), "\n";
    
    for my $o1 ( @names )
    {
	my @values = ($o1);
	
	for my $o2 ( @names )
	{
	    my $count = 0;
	    
	    $count += $matrix{$o1}{$o2} if( exists $matrix{$o1}{$o2} );
	    $count += $matrix{$o2}{$o1} if( exists $matrix{$o2}{$o1} && $o_symmetric );
	    
	    push @values, $count;
	}
	
	print $fh join(",", @values, $legend{$o1}), "\n";
    }
}

sub count_assembly_frequencies {
    my ($data,$fh) = @_;

    my $total = 0;
    my %assemblies = ();

    for my $read ( sortccs(keys %$data) )
    {
	### sort fragments by location
	my @fragments = sort { $$a{"smin"} <=> $$b{"smin"} } @{$$data{$read}};

	my @pieces = ();
	for( my $i = 0; $i < @fragments; $i++ )
	{
	    ### accumulate inserts comprising an assembly
	    push @pieces, {
		"qseqid" => $fragments[$i]{"qseqid"},
		"sstrand" => $fragments[$i]{"sstrand"}
	    };

	    ### reached the last inserts or break point
	    if( $i+1 == scalar(@fragments) || $fragments[$i]{"brk"} == 1 )
	    {
		### encode assembly as a string of inserts {A;B;C}
		my ($gga,$correct) = assembly_encoding(\@pieces);
		
		$assemblies{$gga}{"Count"}++;
		$assemblies{$gga}{"Correct"} = $correct;
		$assemblies{$gga}{"Size"} = scalar(@pieces);

		$total++;

		### reset
		@pieces = ();
	    }
	}
    }

    ### print out assemblies stats
    print $fh "Assembly,Count,Correct,Size,Fraction\n";
    
    for my $gga ( sort { $assemblies{$b}{"Count"} <=> $assemblies{$a}{"Count"} } keys %assemblies )
    {
	printf( $fh "{%s},%i,%i,%i,%f\n",
		$gga,
		$assemblies{$gga}{"Count"},
		$assemblies{$gga}{"Correct"},
		$assemblies{$gga}{"Size"},
		$assemblies{$gga}{"Count"} / $total );
    }

    return \%assemblies;
}

sub assembly_encoding {
    my ($data) = @_;

    my %counts = ( "minus" => 0, "plus" => 0 );
    my @delta = ();

    for( my $i = 0; $i < @$data; $i++ )
    {
	$counts{ $$data[$i]{"sstrand"} }++;

	if( $i + 1 < @$data )
	{
	    ### get insert serial numbers
	    my $n1 = $MAP{$$data[$i]{"qseqid"}};
	    my $n2 = $MAP{$$data[$i+1]{"qseqid"}};

	    push @delta, ($n1-$n2);
	}
    }

    my $correct = 1;

    ### mixture of differently oriented inserts: this is not a correct assembly
    ### {A;B;C;C-;B-;A-}
    if( $counts{"minus"} != 0 && $counts{"plus"} != 0 )
    {
	$correct = 0;
    }

    ### check insert pairs
    ###
    ### {A;B;C;D;E;F} : correct
    ###  1 2 3 4 5 6
    ###   1 1 1 1 1
    ###
    ### {D-;C-;B-;A-} : correct
    ###  4  3  2  1
    ###   -1 -1 -1
    ###
    ### {A;B;C;C-;B-} : incorrect
    ###  1 2 3 3  2
    ###   1 1 0 -1

    for( my $i = 0; $i < @delta; $i++ )
    {
	if( ( $counts{"minus"} > 0 && $delta[$i] != 1 ) || ( $counts{"plus"} > 0 && $delta[$i] != -1 ) )
	{
	    $correct = 0;
	}
    }

    ### try top keep all assemblies uniform
    ### translate "minus" to "plus" and update orientation
    my @temp = @$data;

    if( $counts{"minus"} > $counts{"plus"} )
    {
	@temp = reverse(@temp);
	
	for( my $i = 0; $i < @temp; $i++ )
	{
	    if( $temp[$i]{"sstrand"} eq "minus" )
	    {
		$temp[$i]{"sstrand"} = "plus";
	    }
	    elsif( $temp[$i]{"sstrand"} eq "plus" )
	    {
		$temp[$i]{"sstrand"} = "minus";
	    }
	}
    }

    my @values = ();
    
    for( my $i = 0; $i < @temp; $i++ )
    {
	my $n = $temp[$i]{"qseqid"};

	if( $temp[$i]{"sstrand"} eq "minus" )
	{
	    push @values, "$n-";
	}
	else
	{
	    push @values, "$n";
	}
    }

    my $encoding = join(":",@values);

    return ($encoding,$correct);
}

sub aggregate_assemblies_by_size {
    my ($assemblies,$fh) = @_;

    my %data = ();

    my @sorted = sort { $$assemblies{$b}{"Size"} <=> $$assemblies{$a}{"Size"} } keys %$assemblies;

    for( my $i = 1; $i <= $$assemblies{$sorted[0]}{"Size"}; $i++ )
    {
	$data{$i}{"Total"} = 0;
	$data{$i}{"Correct"} = 0;
    }

    my $total = 0;

    for my $gga ( keys %$assemblies )
    {
	$total += $$assemblies{$gga}{"Count"};

	$data{$$assemblies{$gga}{"Size"}}{"Total"} += $$assemblies{$gga}{"Count"};
	
	if( $$assemblies{$gga}{"Correct"} )
	{
	    $data{$$assemblies{$gga}{"Size"}}{"Correct"} += $$assemblies{$gga}{"Count"};
	}
    }

    print $fh "AssemblySize,Total,Correct,TotalFraction,CorrectFraction\n";
    
    for my $size ( sort { $a <=> $b } keys %data )
    {
	printf( $fh "%i,%i,%i,%f,%f\n",
		$size,
		$data{$size}{"Total"},
		$data{$size}{"Correct"},
		$data{$size}{"Total"} / $total,
		$data{$size}{"Correct"} / $total,
	    );
    }
}

sub assemblies_A_to_J {
    my ($assemblies,$fh) = @_;

    my %AJ = ();
    my $total = 0;
    my $correct = 0;

    ### iterate through all assemblies
    for my $gga ( keys %$assemblies )
    {
	my @tokens = split(/:/, $gga);
	
	my $b = $tokens[0];  $b =~ s/-//g;
	my $e = $tokens[-1]; $e =~ s/-//g;

	my %set = ( "A" => 1, "J" => 1 );

	### find assemblies that start with A and end with J
	### (irrespective of direction)
	if( exists $set{$b} && exists $set{$e} && $b ne $e )
	{
	    $AJ{$gga} = $$assemblies{$gga};

	    $total += $$assemblies{$gga}{"Count"};

	    if( $$assemblies{$gga}{"Correct"} )
	    {
		$correct += $$assemblies{$gga}{"Count"};
	    }
	}
    }

    print $fh "Assembly,Count,Correct,Size,Fraction\n";

    for my $gga ( sort { $AJ{$b}{"Correct"} <=> $AJ{$a}{"Correct"} || $AJ{$b}{"Count"} <=> $AJ{$a}{"Count"} } keys %AJ )
    {
	printf( $fh "{%s},%i,%i,%i,%f\n",
		$gga,
		$AJ{$gga}{"Count"},
		$AJ{$gga}{"Correct"},
		$AJ{$gga}{"Size"},
		$AJ{$gga}{"Count"} / $total );
    }
}

sub sortccs {
    my (@arr) = @_;

    ### PacBio CCS read naming convention: {movieName}/{holeNumber}/ccs
    ###
    ### Example:
    ### m170525_235416_42200_c101220802550000001823278510171700_s1_p0/163469/ccs

    ### extract name parts
    my @tokens = map { [split(/\//,$_)] } @arr;

    ### sort by movieName and holeNumber
    my @sorted = @arr[ sort { $tokens[$a][0] cmp $tokens[$b][0] || $tokens[$a][1] <=> $tokens[$b][1] } 0..$#arr ];

    return @sorted;
}

sub complement {
    my ($seq) = @_;

    my %bases = (
	"A" => "T",
	"C" => "G",
	"G" => "C",
	"T" => "A",
	"a" => "t",
	"c" => "g",
	"g" => "c",
	"t" => "a",
	"-" => "-",
	);

    my $complement = join("", map {$bases{$_}} split(//, $seq));

    return $complement;
}
