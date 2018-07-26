#!/bin/bash

inputxml=semi-leptonic.xml

###======= 500 GeV =======

##=====new 6f_ttbar samples
#6f_ttbar/yyxylv.eL.pR
#dirin=/hsm/ilc/grid/storm/prod/ilc/mc-dbd/ild/dst/500-TDR_ws/6f_ttbar/ILD_o1_v05/v01-16-p05_500/
dirin=/hsm/ilc/grid/storm/prod/ilc/mc-dbd/ild/dst-merged/500-TDR_ws/6f_ttbar/ILD_o1_v05/v01-16-p05_500/
processID=$1
#nfile=-1
nfile=12
#nperjob=50
nperjob=1
./submit2 $dirin $processID 1 $nfile $nperjob $inputxml
