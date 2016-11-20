#Title:		MirrorFASTA
#Details:	Get the 'mirror' sequences of a FASTA file
#Usage:		MirrorFASTA.sh [INPUT] [OUTPUT]
#Author:		Kevin Blighe
#Date:		9th March 2016

if [ $# -ne 2 ] ;
then
	echo "Illegal number of parameters"
	echo "Usage:   MirrorFASTA.sh [INPUT] [OUTPUT]"
	exit 1
fi

count=0

touch $2

while read line           
do
	((count+=1))

	if ! (($count%2)) ;
	then
		echo $line | rev >> $2
	else
		echo $line >> $2
	fi
done < $1
