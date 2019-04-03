#!/usr/local/bin/bash

sylbin=$HOME/sylamer-ca

 
set -e
set -o pipefail

PS4="Running command: "

mkv=4
thedir=syldir
SYLOPTS="-tt 0 -grow 200 --alien-ok -v 0 -m $mkv"
SYLUSEROPTS=""
SYLPLOTOPTS=""
stem=syl
Pvaluecutoff=0
do_clean_up=false
NUMTOPEVD=3

thefastafile=""  # need-fasta-file-here
therankfile=""   # need-rankfile-here
thetable=""
ylim=""
theoutput="out.syl"
thequantile=0.002
thewords=""
thespecies="hsa"
thetitle=""
force=false
mirna_select=
z=--


while getopts :f:w:T:E:p:P:S:t:y:m:r:d:o:q:s:x:z:WFCh opt
do
    case "$opt" in
    m)
      mkv=$OPTARG
      SYLUSEROPTS="$SYLUSEROPTS -m $mkv"
      ;;
    f)
      thefastafile=$OPTARG
      ;;
    x)
      SYLUSEROPTS="$SYLUSEROPTS $OPTARG"
      ;;
    w)
      thewords=$OPTARG
      ;;
    d)
      thedir=$OPTARG
      ;;
    y)
      ylim=$OPTARG
      ;;
    t)
      thetable=$OPTARG
      ;;
    q)
      thequantile=$OPTARG
      ;;
    s)
      stem=$OPTARG
      ;;
    r)
      therankfile=$OPTARG
      ;;
    z)
      z=$OPTARG
      ;;
    p)
      Pvaluecutoff=$OPTARG
      ;;
    T)
      thetitle=$OPTARG
      ;;
    S)
      thespecies=$OPTARG
      ;;
    o)
      theoutput=$OPTARG
      ;;
    W)
      SYLPLOTOPTS="$SYLPLOTOPTS --top=0 --bot=0"
      ;;
    P)
      SYLPLOTOPTS="$SYLPLOTOPTS $OPTARG"
      ;;
    E)
      NUMTOPEVD=$OPTARG
      ;;
    F)
      force=true
      ;;
    C)
      mirna_select=CONS
      ;;
    h)
      cat <<EOU
Options marked with * are required:

-f <fastafile>  *  
-r <rankfile>   *
-d <directory>  *
-t <tablefile>  * (made by syl.mirna-import.pl from Cei files)
-m <num>          Markov correction (default $mkv)
-q <num>          quantile to discard for EVD estimation
-s <stem>         stem for intermediate files)
-S <species>      as used in mirna identifiers; e.g. mmu, hsa
-T <title>        plot title
-w word,word,..   microRNA words to highlight
-p <num>          Pvalue cutoff
-y <ylim>         Y-axis maximum
-P "plot options" passed to syl.plot.R
-E <num>          plot this many top words from EVD analysis. Default $NUMTOPEVD.
-W                turn off automatic highlighting of top words
-F                force sylamer re-runs
-C                select more stringent MIRNA_CONS (conserved) microRNAs as microRNA words only
-x "options"      sylamer options passed through to sylamer
-z tag            skip code (string) matching 'tag'. (statevd|plottop|plotevd).
EOU
      exit
      ;;
    :) echo "Flag $opt needs argument"
        exit 1;;
    ?) echo "Flag $opt unknown"
        exit 1;;
   esac
done



if [[ -z $thetable || -z $thefastafile || -z $therankfile ]]; then
   echo "need -f FASTAFILE -r RANKFILE -t THETABLE"
   false
fi

if [[ ! $thefastafile =~ ^/ ]]; then
   thefastafile=$PWD/$thefastafile
fi
therankfile=$PWD/$therankfile
thetable=$PWD/$thetable


if [[ -z $ylim ]]; then
   ylim=10
fi


function clean_up {
   { set +x; } 2> /dev/null
   cd $OLDPWD
   if [[ $? == 0 ]]; then
      echo "[happy] Script succeeded"
      true
   else
      echo "[grumpy] Script failed"
      true
   fi
   if $do_clean_up; then
      rm -rf $thedir
      echo "I've removed the working directory $thedir"
   else
      echo "I've left the working directory $thedir in place"
   fi
}


trap clean_up SIGTERM  EXIT

mkdir -p $thedir
cd $thedir


echo -e "\n\n========= Run n $thedir ========="
echo "$0 $@" > INVOCATION

set -e

if [[ -z $thetitle ]]; then
   thetitle=${therankfile##/}
fi


##
##    sylamer runs plus plotting of top words.
#
for k in 6 7 8; do
   # echo doing sylamer run $k
   words=""
   if [[ $k == "8" ]]; then
     words="-w .......A"
   fi
   set -x
   if $force || [[ ! -e $stem.$k.gz ]]; then
      sylamer -fasta $thefastafile -universe $therankfile $words -k $k -o - $SYLOPTS $SYLUSEROPTS  | gzip > $stem.$k.gz;
   else
      echo "File $stem.$k.gz exists, skipping"
   fi
   if [[ ! $z =~ plottop ]]; then
      ( R --vanilla --slave --quiet --args --data=$stem.$k.gz --ymin=-$ylim --ymax=$ylim --title="$thetitle" --table="$thetable" --pdf=pdf.$stem.$k.hltop.pdf --add=$thewords --species=$thespecies $SYLPLOTOPTS < $sylbin/syl.plot.R ) >> RLog
      { (( ${PIPESTATUS[0]} != 0 )) && false; set +x; } 2>/dev/null
   fi
done


set -x
if $force || [[ ! -e $stem.SUM4-group ]]; then
   syl.maketable.pl --syl6="$stem".6.gz --syl7="$stem".7.gz --syl8="$stem".8.gz --range=1 --table="$thetable" --mode=SUM4 > $stem.SUM4-group
else
   { set +x; } 2>/dev/null
   echo "File $stem.SUM4-group exists, skipping"
fi


set -x
if [[ ! $z =~ statevd ]]; then
   if ! R --slave --quiet --args --quantile=$thequantile --table="$thetable" --filter="$Pvaluecutoff" --sylsum="$stem".SUM4-group < $sylbin/syl.fitevd.R > "$stem".SUM4-evd; then
      { set +x; } 2>/dev/null
      echo "syl.fitevd.R failed"
      false
   else
      { set +x; } 2>/dev/null
      if $do_clean_up; then
         echo "Output copied to $OLDPWD/$theoutput"
         cp $stem.SUM4-evd $OLDPWD/$theoutput
      fi
   fi
fi


#  TODO: without the +e the process substitution below triggers clean_up, which
#  still claims to be happy.  mmm. is that because the trap function clean_up
#  is triggered by the subprocess EXIT?

set +e
lsowords=$(grep NATIVE_$mirna_select "$stem".SUM4-evd | head -n $NUMTOPEVD | cut -f 1 | tr '\n' ',' | perl -0pe 's/,\Z//')
legend=$(grep -n NATIVE_$mirna_select "$stem".SUM4-evd | head -n $NUMTOPEVD | cut -f 1,8 | perl -ane '$F[0] =~ /(\d+):(\w+)/;printf "%s evd %.1f (rank %d) ## ", $2, -log($F[1])/log(10), $F[0]-1')
set -e

if [[ -z $lsowords ]]; then
   echo "Could not select words from EVD analysis. Please investigate."
else
   echo "Plotting selected words from EVD analysis: $lsowords"
   echo "Using legend $legend"

   for k in 6 7 8; do
      set -x
      if [[ ! $z =~ plotevd ]]; then
         ( R --vanilla --slave --quiet --args --data=$stem.$k.gz --ymin=-$ylim --ymax=$ylim --title="$thetitle" --table="$thetable" --pdf=pdf.$stem.$k.hlevd.pdf --add=$lsowords --species=$thespecies $SYLPLOTOPTS --second-legend="$legend" --top=0 --bottom=0 < $sylbin/syl.plot.R ) >> RLog
         { (( ${PIPESTATUS[0]} != 0 )) && false; set +x; } 2>/dev/null
      fi
   done
fi




