#!/usr/local/bin/perl

# usage: perl bigrun.pl 

$wrkdir = "/star/scratch/starkong/STAR_FlowAnalyzer/";

#system("mkdir tmpAMPT");
#system("mkdir tmpAMPT/log");

# open file as read_only
$filelist ="AuAu54.list";
open(LIST,"< $filelist") or die "Can not open ${filelist} !\n";
@list=<LIST>;

for($i=0; $i<50; $i++) {

    $onefile =$list[$i];

    print $i; print "\n";

    print $onefile; print "\n";
    chomp($onefile);

    system("mkdir $onefile");
    system("mkdir rootfile/$onefile");
    
    chdir("$onefile");
    system("mkdir $onefile");
    system("ln -s $wrkdir/.sl73_gcc485 .");
    system("ln -s $wrkdir/StRoot .");
    system("cp $wrkdir/doEvent.C .");
    system("cp $wrkdir/Analysis.xml .");

# submit job to Condor
   system("star-submit-template -template Analysis.xml -entities run=$onefile");

    chdir("..");
# to sleep in seconds		
   sleep 1;

    
}  # end of while file_list

# to close the file handles
close LIST;

