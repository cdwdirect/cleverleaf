#!/usr/bin/env perl

# Build a tab-delineated summary table of timings for plotting time vs nproc
# using plot-samrai-timers.


use ExtractSAMRAITimers;


require Getopt::Long;
my $raw;
Getopt::Long::GetOptions( 'help' => \$help
	  , 'verbose' => \$verbose
	  , 'debug' => \$debug
          , 'raw' => \$raw
	  );

my $script_name = `basename $0`; chop $script_name;

if ( $help ) {
  print <<_EOM_;

Usage: $script_name [options] <log files>

Log files must be of the form *-<nproc>.log*
Each log file name is parsed to find its corresponding nproc.
The same set of TreeCommunication benchmark tests must appear
in all the log files, in the same order.  (They must all be
generated using the same input files.)

Options:

  --raw
    Use raw values instead of normalizing by number of repetitions.

_EOM_
  exit(0);
}

$ExtractSAMRAITimers::verbose = $verbose;
$ExtractSAMRAITimers::debug = $debug;

my @logs = @ARGV;
die "No log files specified" if @logs == 0;

my @nprocs = ExtractSAMRAITimers::get_nprocs_from_log_names(@logs);
$,=' ';
print STDERR ("nprocs:", @nprocs, "\n") if $verbose;
$,='';

$, = "\t";


###################################################################
# Suck in all the timer sections.
###################################################################

$/ = "\n\n";
my @timerparagraphs;
if (defined $timerfile) {
  open( TIMERFILE, "< $timerfile" ) || die "Cannot open file $timerfile";
  @timerparagraphs = <TIMERFILE>;
  close TIMERFILE || die "Cannot close $timerfile";
}
else {
  @timerparagraphs = <DATA>;
}
$/ = "\n";
my ($section_name, $timers_, %timer_label) =
    ExtractSAMRAITimers::parse_timer_selection($timerparagraphs[0]);
my @timers = @$timers_;



###################################################################
# Scan first log file to determine number of tests.
# We have to require that the same tests appear in all logs.
# Store the (nick)names in @test_names.
###################################################################
my @test_names;
open (LOG, $logs[0]) || die "Cannot open $logs[0]";
while ( $_ = <LOG> ) {
  if ( /^Starting test Test([0-9]+) \((.*)\)$/ ) {
    $test_names[$1] = $2;
  }
}
close LOG;
if ( $verbose ) {
  print STDERR "test names:\n";
  for ( 0..$#test_names ) {
    print STDERR "\t$_: $test_names[$_]\n";
  }
}


###################################################################
# Scan log files to collect all the timers.
# Store them in %times, which maps test_name+timer_name to
# the timer value.
###################################################################
my %times;
ExtractSAMRAITimers::extract_timers( \@nprocs, \@logs, \@test_names, \%times );



###################################################################
# Write the output file.
# One paragraph for each test.
###################################################################

# Compose the table headings string.
my $table_heading = "nproc";
for ( 0 .. $#timers ) {
  $table_heading = "${table_heading}\t$timer_label{$timers[$_]}";
}

$,='';

for $itest (0 .. $#test_names) {
  # Each chart is for a named test

  print "$test_names[$itest]\n";
  print "timer*nproc:\n"; # Unused, but needed.
  print "$table_heading\n";

  for $inproc (0 .. $#nprocs) {
    # Each datum is for an nproc

    print $nprocs[$inproc];
    for $itimer (0 .. $#timers) {
      # Each line is for a timer.
      my $value = defined($times{"$itest:$timers[$itimer]"}[$inproc]) ?
          $times{"$itest:$timers[$itimer]"}[$inproc] : 'x';
      print "\t$value";
    }
    print "\n";

  }
  print "\n";

}

# The data section should have paragraphs delinineated by blank lines.
# First line in each paragraph is the header.  If header is commented
# out, skip paragraph.  Each following line in paragraph is a timer
# name.  If timer name is commented out, skip that timer.  If the
# timer name is followed by <tab> and a string, use the timer, but
# substitute the string for timer name.

__DATA__
All timers
#
apps::main::repetitions	Repetition loop
#
apps::main::parent_send	Parent send
apps::main::parent_recv	Parent receive
apps::main::parent_wait	Parent wait
#
apps::main::child_send	Child send
apps::main::child_recv	Child receive
apps::main::child_wait	Child wait
#
apps::main::verify	Verify data
apps::main::processing	Process data
#
apps::main::equiv_MPI	MPI equivalent
