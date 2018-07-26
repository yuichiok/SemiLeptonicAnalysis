#!/bin/bash

#bfaf=$1

#if [ -z $bfaf ] ; then
#  hadd -f root_merge/hadronic_eLpR_kekcc_500k.root rootfile/kekcc/rv01-16-p05_500*I108679*.root rootfile/kekcc/rv01-16-p05_500*I108681*.root
#else
#  hadd -f root_merge/aftervertexrestore_hadronic_eLpR_kekcc_500k.root rootfile/kekcc/aftervertexrestore*I108679*.root rootfile/kekcc/aftervertexrestore*I108681*.root
#fi

hadd -f root_merge/leptonic_yyxylv_eLeR.root rootfile/*I108675*
