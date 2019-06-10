#! /usr/bin/perl
##
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/trunk/tools/scripts/conversion2.0/renameXd.pl $
## Package:     SAMRAI scripts
## Copyright:   (c) 1997-2018 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 1917 $
## Description: perl script to rename getNumber methods to be getNumberOf 
##

use strict;

use File::Basename;
use File::Find;
use File::Compare;
use Cwd;

my $SAMRAIBinDir = $ARGV[0];
my $replaceDir   = $ARGV[1];

my $SAMRAISourceDir = $SAMRAIBinDir . "/../../..";
my $SAMRAIScriptDir = $SAMRAISourceDir . "/tools/scripts/conversion3.0";

# Flush I/O on write to avoid buffering
# $|=1;

my $debug=4;

my $end_of_line = $/;

#
# Remove duplicated values
#
sub unique {
    foreach my $test (@_){
	my $i = -1;
	my @indexes = map {$i++;$_ eq $test ? $i : ()} @_;
	shift @indexes;
	foreach my $index (@indexes){
	    splice(@_,$index,1);
	}
    }
    return @_;
}

my $pwd = cwd;

my @templatesOnDIMOnly = ();
my $filename=$SAMRAIScriptDir . "/SAMRAI_classes_templated_on_DIM_only.txt";
open FILE, '<', $filename or die "Can't open $filename : $!";
while (<FILE>) {
    my $class=$_;
    chomp($class);
    push @templatesOnDIMOnly,$class;
}
close FILE;


my @templatesOnDIM = ();
my $filename=$SAMRAIScriptDir . "/SAMRAI_classes_templated_on_DIM_and_other.txt";
open FILE, '<', $filename or die "Can't open $filename : $!";
while (<FILE>) {
    my $class=$_;
    chomp($class);
    push @templatesOnDIM,$class;
}
close FILE;

#=============================================================================
# Fixup for all source
#=============================================================================

#
# Excludes files that are in internal source code control directories.
#
my $filePattern = q/(.*\.[ChI]$)|(.*\.CPP$)|(.*\.cpp$)|(.*\.cxx$)|(.*\.CXX$)|(.*\.H$)|(.*\.hxx$)|(.*\.Hxx$)|(.*\.HXX$)|(.*\.txx$)|(.*\.ixx$)/;
my @filesToProcess = ();
sub selectAllSourceFiles {
    if ( $File::Find::name =~ m!/(.svn|CVS|include|scripts|\{arch\})$!o ) {
	$File::Find::prune = 1;
    }
    elsif ( -f && m/$filePattern/ ) {
	push @filesToProcess, $File::Find::name;
	$filesToProcess[$#filesToProcess] =~ s|^\./||o;
    }
}

find( \&selectAllSourceFiles, '.' );

print "@filesToProcess\n" if ($debug);

for my $file (@filesToProcess) {
    print "Working on DIM fixes for source file $file\n";
    my $directory = dirname $file;

    my $filebasename = basename $file;

    my $tempFile = $filebasename . ".samrai.tmp";

    open FILE, "< $file" || die "Cannot open file $file";
    open TEMPFILE, "> $tempFile" || die "Cannot open temporary work file $tempFile";
    while ( my $str = <FILE> ) {

	# Replace variable definitions
	for my $classname (@templatesOnDIMOnly) {
	    $str =~ s/$classname\<DIM\>/$classname/g;
	    $str =~ s/$classname\<NDIM\>/$classname/g;
	    $str =~ s/$classname\<2\>/$classname/g;
	    $str =~ s/$classname\<3\>/$classname/g;
	}
	
	for my $classname (@templatesOnDIM) {
	    $str =~ s/$classname<DIM,\s*(.*)>/$classname<$1>/g;
	    $str =~ s/$classname<DIM,(.*)>/$classname<$1>/g;
	    $str =~ s/$classname<NDIM,\s*(.*)>/$classname<$1>/g;
	    $str =~ s/$classname<2,\s*(.*)>/$classname<$1>/g;
	    $str =~ s/$classname<3,\s*(.*)>/$classname<$1>/g;
	}

	print TEMPFILE $str;
    }

    close FILE || die "Cannot close file $file";
    close TEMPFILE || die "Cannot close file $tempFile";

    # Only replace existing file if a replacement was done.
    if (compare($file,$tempFile) == 0) {
	unlink($tempFile);
    } else {
	unlink($file);
	rename( $tempFile, $file);
    }
}

