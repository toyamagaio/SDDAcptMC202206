#!/bin/sh

BeDi=55
Nev=100000
output_base="SDDAcptMC_pmuTest_SDD1_"
#./bin/SDDAcptMC -n $Nev -d $BeDi -w root/${output_base}d${BeDi}_n${Nev}_z25.root        -p pdf/${output_base}d${BeDi}_n${Nev}_z25.pdf        -t 0
#./bin/SDDAcptMC -n $Nev -d $BeDi -w root/${output_base}d${BeDi}_n${Nev}_z25_pmuuni.root -p pdf/${output_base}d${BeDi}_n${Nev}_z25_pmuuni.pdf -t 1
#./bin/SDDAcptMC -n $Nev -d $BeDi -w root/${output_base}d${BeDi}_n${Nev}_z25_pmuGau.root -p pdf/${output_base}d${BeDi}_n${Nev}_z25_pmuGau.pdf -t 2


sddx_array=(-20 -15 -10 -5 0 5 10)

for sddx in ${sddx_array[@]}; do
  ./bin/SDDAcptMC -n $Nev -d $BeDi -w root/${output_base}d${BeDi}_n${Nev}_z25_pmuGau_SDDx${sddx}.root -p pdf/${output_base}d${BeDi}_n${Nev}_z25_pmuGau_SDDx${sddx}.pdf -t 2 -x $sddx
done
