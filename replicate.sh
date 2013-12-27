#!/bin/sh

# MGEL
# Surya Saha 05/30/2008 
# Cmd line arg 1 : <input .gff file> 
# Cmd line arg 2 : <# of noise files> 
# Cmd line arg 3 : <length of chromosome>
# to create noise files and mine each one for relationships

if [ $# -ne 3 ]
then
	echo "Usage $0 <input .gff file> <noise file index> <length of chromosome>"
	exit 1
fi

# creating the noise files
for ((i=0;i<$2;i++))
do
	temp=`expr $i + 1`
	/home/ssaha/currwork/rulemining/release/miner.gff2noise.v1.pl $1 $temp $3
	echo "Noise file $temp created.."
done

# finding the relationships
#read noise FILES
FILES="$(echo *.noise.*)"

for FILE in $FILES
do
  /home/ssaha/currwork/rulemining/release/miner.gff2f_itemsets.stage1.v1.pl 15000 $FILE -nocop 1>/dev/null 2> /dev/null
  echo "$FILE mined.."
done
