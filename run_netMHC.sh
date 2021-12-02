#!/bin/bash

while read line
do

file="HLA/"$line"_HLA_Alleles.txt"
hla=`cat $file`

echo $hla

netMHCpan -p "Peptides/"$line"_Peptides.txt" -a $hla -xls -xlsfile "output/"$line"_netMHC.xls"

done < $1
