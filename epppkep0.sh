#!/bin/bash 
role=0
pairs=$(seq 2 1 5)
size=$(seq 2 1 4)

for i in $pairs; do
	echo "pair $i"
	for j in $size; do
		echo "cycle size: $j"
		./eppkep -r $role -f 0 -x $i -z $size > outfile_r${role}_p${i}_c${j}.txt
	done;
done

for i in $pairs; do
	echo "pair $i"
	for j in $size; do
		echo "cycle size: $j"
		./eppkep -r $role -f 1 -x $i -z $size > outfile_all_r${role}_p${i}_c${j}.txt
	done;
done

echo "success";