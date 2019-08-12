# This file is configured by CMake automatically as DartConfiguration.tcl
# If you choose not to use CMake, this file may be hand configured, by
# filling in the required variables.


# Configuration directories and files
SourceDirectory: /g/g17/wood67/src/cleverleaf/SAMRAI
BuildDirectory: /g/g17/wood67/src/cleverleaf/apollo-test/RelWithDebInfo/build/SAMRAI

# Where to place the cost data store
CostDataFile: 

# Site is something like machine.domain, i.e. pragmatic.crd
Site: quartz1532

# Build name is osname-revision-compiler, i.e. Linux-2.4.2-2smp-c++
BuildName: Linux-g++

# Subprojects
LabelsForSubprojects: 

# Submission information
SubmitURL: http://

# Dashboard start time
NightlyStartTime: 00:00:00 EDT

# Commands for the build/test/submit cycle
ConfigureCommand: "/g/g17/wood67/src/spack/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/cmake-3.14.2-e2vmvqwgjtcnpj6l5vlr3wmenrqhkqle/bin/cmake" "/g/g17/wood67/src/cleverleaf/SAMRAI"
MakeCommand: /g/g17/wood67/src/spack/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/cmake-3.14.2-e2vmvqwgjtcnpj6l5vlr3wmenrqhkqle/bin/cmake --build . --config "${CTEST_CONFIGURATION_TYPE}" -- -i
DefaultCTestConfigurationType: Release

# version control
UpdateVersionOnly: 

# CVS options
# Default is "-d -P -A"
CVSCommand: /usr/bin/cvs
CVSUpdateOptions: -d -A -P

# Subversion options
SVNCommand: /usr/bin/svn
SVNOptions: 
SVNUpdateOptions: 

# Git options
GITCommand: /usr/tce/bin/git
GITInitSubmodules: 
GITUpdateOptions: 
GITUpdateCustom: 

# Perforce options
P4Command: P4COMMAND-NOTFOUND
P4Client: 
P4Options: 
P4UpdateOptions: 
P4UpdateCustom: 

# Generic update command
UpdateCommand: 
UpdateOptions: 
UpdateType: 

# Compiler info
Compiler: /usr/tce/packages/gcc/gcc-4.9.3/bin/g++
CompilerVersion: 4.9.3

# Dynamic analysis (MemCheck)
PurifyCommand: 
ValgrindCommand: 
ValgrindCommandOptions: 
MemoryCheckType: 
MemoryCheckSanitizerOptions: 
MemoryCheckCommand: /usr/tce/bin/valgrind
MemoryCheckCommandOptions: --trace-children=yes --leak-check=full
MemoryCheckSuppressionFile: 

# Coverage
CoverageCommand: /usr/tce/packages/gcc/gcc-4.9.3/bin/gcov
CoverageExtraFlags: -l

# Cluster commands
SlurmBatchCommand: /usr/bin/sbatch
SlurmRunCommand: /usr/bin/srun

# Testing options
# TimeOut is the amount of time in seconds to wait for processes
# to complete during testing.  After TimeOut seconds, the
# process will be summarily terminated.
# Currently set to 25 minutes
TimeOut: 1500

# During parallel testing CTest will not start a new test if doing
# so would cause the system load to exceed this value.
TestLoad: 

UseLaunchers: 
CurlOptions: 
# warning, if you add new options here that have to do with submit,
# you have to update cmCTestSubmitCommand.cxx

# For CTest submissions that timeout, these options
# specify behavior for retrying the submission
CTestSubmitRetryDelay: 5
CTestSubmitRetryCount: 3
