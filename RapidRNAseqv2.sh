#!/bin/bash

###########################################################
#Pre-liminary step 1:  Safeguarding the script's operation#
###########################################################

#Check if 0 parameters have been passed
if [ $# -eq 0 ]
then
	echo -e "\n"
	echo "Error - incorrect number of parameters.  Consult the relevant standard operating procedure for correct usage."
	echo "Please type RapidRNAseqv2.sh -h for help."
	echo -e "\n"
	exit 1
fi

#Check if the help parameter was passed
if [ "${1}" == "-h" ]
then
	echo -e "##################################################################################################################"
	echo -e "RapidRNAseq - Kallisto wrapper, developed by Kevin Blighe (k.blighe@qmul.ac.uk) at Queen Mary University of London"
	echo -e "Release date:\t\tThursday 5th May 2016"
	echo -e "Version:\t\tVersion 2"
	echo -e "Mods to version 2:\tNow supports the analysis of sample replicates concurrently"
	echo -e "Mods to version 2:\tRemoved the necessity to supply output file prefix at the command-line"
	echo -e "##################################################################################################################"
	echo -e "\n"


	echo -e "Program syntax:"
	echo -e "Paired-end"
	echo -e "RapidRNAseqv2.sh PAIRED [Rep1MatePair1,Rep1MatePair2,Rep2MatePair1,Rep2MatePair2,...,] [Reference FASTA/FASTA.gz] [Bootstrap quantification value] [Output dir] [User initials]"
	echo -e "\n"

	echo -e "Single-end"
	echo -e "RapidRNAseqv2.sh SINGLE [Fragment/Read length] [Fragment/Read length standard deviation] [Rep1,Rep2,...,] [Reference FASTA/FASTA.gz] [Bootstrap quantification value] [Output dir] [User initials]"
	echo -e "\n"

	exit 0
fi

#Check that no less than 6 parameters have been passed - 6 for paired is minimum possible needed (8 for single)
if [ $# -lt 6 ]
then
	echo -e "\n"
	echo "Error - incorrect number of parameters.  Consult the relevant standard operating procedure for correct usage."
	echo "Please type RapidRNAseqv2.sh -h for help."
	echo -e "\n"
	exit 1
fi

#Check if the first parameter is for single or paired-end sequencing
if [ $(echo "${1}" | awk '{print toupper($0)}') == "PAIRED" ]
then
	#Check that 6 parameters have been passed
	if [ $# -ne 6 ]
	then
		echo -e "\n"
		echo "Error - incorrect number of parameters for paired-end analysis.  Consult the relevant standard operating procedure for correct usage."
		echo "Please type RapidRNAseqv2.sh -h for help."
		echo -e "\n"
		exit 1
	fi

	###

	#Check that the file(s) as $2 exist
	declare -a arraySequences
	declare -i iNumSequences

	#Break up the files into an array
	#sequences=`echo "${2}" | sed 's/,/ /g'`
	IFS=',' read -r -a arraySequences <<< "${2}"

	#Get the number of array elements (i.e. files)
	iNumSequences=${#arraySequences[@]}

	#Loop throgh and test that each exists
	for ((i=0;i<iNumSequences;i++))
	do
		test -e "${arraySequences[i]}"

		if [ $? -ne 0 ]
		then
			echo -e "\n"
			echo "Error - input file "${arraySequences[i]}" does not exist.  Please check the complete file path and re-run."
			echo "Note:  input sample files can be in FASTQ or FASTQ.gz format."
			echo "Please type RapidRNAseqv2.sh -h for help."
			echo -e "\n"
			exit 1
		fi
        done

	###

	#Check that the file as $3 exists
	test -e "${3}"
	if [ $? -ne 0 ]
	then
		echo -e "\n"
		echo "Error - reference file ${3} does not exist.  Please check the complete file path and re-run."
		echo "Note:  reference sequences can be in FASTA or FASTA.gz format"
		echo "Please type RapidRNAseqv2.sh -h for help."
		echo -e "\n"
		exit 1
	fi

	#Check that $4 is numerical
	if ! [ "${4}" -eq "${4}" ] 2>/dev/null
	then
		echo -e "\n"
		echo "Error - bootstrap value ${4} must be an integer value."
		echo "Please type RapidRNAseqv2.sh -h for help."
		echo -e "\n"
		exit 1
	fi

	#Check that $5 is not a file
	if [ -f "${5}" ]
	then
		echo -e "\n"
		echo "Error - output directory ${5} already exists and is a regular file."
		echo "Please choose another output directory"
		echo -e "\n"
		exit 1
	fi

	#Check that $5 is a directory
	if [ ! -d "${5}" ]
	then
		echo -e "\n"
		echo "Output directory ${5} does not exist."
		echo "Creating..."
		mkdir "${5}"
	fi
fi

#Check if the first parameter is for single or paired-end sequencing
if [ $(echo "${1}" | awk '{print toupper($0)}') == "SINGLE" ]
then
	#Check that 8 parameters have been passed
	if [ $# -ne 8 ]
	then
		echo -e "\n"
		echo "Error - incorrect number of parameters for single-end analysis.  Consult the relevant standard operating procedure for correct usage."
		echo "Please type RapidRNAseqv2.sh -h for help."
		echo -e "\n"
		exit 1
	fi

	#Check that $2 is numerical
	if ! [ "${2}" -eq "${2}" ] 2>/dev/null
	then
		echo -e "\n"
		echo "Error - estimated fragment/read length ${2} must be an integer value"
		echo "Please type RapidRNAseqv2.sh -h for help."
		echo -e "\n"
		exit 1
	fi

	#Check that $3 is numerical
	if ! [ "${3}" -eq "${3}" ] 2>/dev/null
	then
		echo -e "\n"
		echo "Error - fragment/read length standard deviation ${3} must be an integer value"
		echo "Please type RapidRNAseqv2.sh -h for help."
		echo -e "\n"
		exit 1
	fi

	###

	#Check that the file(s) as $4 exist
	declare -a arraySequences
	declare -i iNumSequences

	#Break up the files into an array
	#sequences=`echo "${4}" | sed 's/,/ /g'`
	IFS=',' read -r -a arraySequences <<< "${4}"

	#Get the number of array elements (i.e. files)
	iNumSequences=${#arraySequences[@]}

	#Loop throgh and test that each exists
	for ((i=0;i<iNumSequences;i++))
	do
		test -e "${arraySequences[i]}"

		if [ $? -ne 0 ]
		then
			echo -e "\n"
			echo "Error - input file "${arraySequences[i]}" does not exist.  Please check the complete file path and re-run."
			echo "Note:  input sample files can be in FASTQ or FASTQ.gz format."
			echo "Please type RapidRNAseqv2.sh -h for help."
			echo -e "\n"
			exit 1
		fi
        done

	###

	#Check that the file as $5 exists
	test -e "${5}"
	if [ $? -ne 0 ]
	then
		echo -e "\n"
		echo "Error - reference file ${5} does not exist.  Please check the complete file path and re-run."
		echo "Note:  reference sequences can be in FASTA or FASTA.gz format"
		echo "Please type RapidRNAseqv2.sh -h for help."
		echo -e "\n"
		exit 1
	fi

	#Check that $6 is numerical
	if ! [ "${6}" -eq "${6}" ] 2>/dev/null
	then
		echo -e "\n"
		echo "Error - bootstrap value ${6} must be an integer value"
		echo "Please type RapidRNAseqv2.sh -h for help."
		echo -e "\n"
		exit 1
	fi

	#Check if $7 is a regular file
	if [ -f "${7}" ]
	then
		echo -e "\n"
		echo "Output directory ${7} already exists and is a regular file."
		echo "Please choose another output directory"
		echo -e "\n"
		exit 1
	fi

	#Check that $7 is a directory
	if [ ! -d "${7}" ]
	then
		echo -e "\n"
		echo "Output directory ${7} does not exist."
		echo "Creating..."
		mkdir "${7}"
	fi

fi

#If neither single nor paired is selected
if [ $(echo "${1}" | awk '{print toupper($0)}') != "SINGLE" ]
then
	if [ $(echo "${1}" | awk '{print toupper($0)}') != "PAIRED" ]
	then
		echo -e "\n"
		echo "Error - unknown sequencing method selected."
		echo "Sequencing method must be either SINGLE or PAIRED."
		echo "Please type RapidRNAseqv2.sh -h for help."
		echo -e "\n"
		exit 1
	fi
fi



#######################################################################
#Pre-liminary step 2:  Directory-structure creation and begin log file#
#######################################################################
#Check if the first parameter is for single or paired-end sequencing
if [ $(echo "${1}" | awk '{print toupper($0)}') == "PAIRED" ]
then
	#Begin log file
	echo "Beginning script on `date`, run by "${6}" with the following parameters:" > "${5}"/AnalysisLog.txt
	echo -e "\t1\t"${1}"" >> "${5}"/AnalysisLog.txt
	echo -e "\t2\t"${2}"" >> "${5}"/AnalysisLog.txt
	echo -e "\t3\t"${3}"" >> "${5}"/AnalysisLog.txt
	echo -e "\t4\t"${4}"" >> "${5}"/AnalysisLog.txt
	echo -e "\t5\t"${5}"" >> "${5}"/AnalysisLog.txt
	echo -e "\t6\t"${6}"" >> "${5}"/AnalysisLog.txt
	echo -e "\n" >> "${5}"/AnalysisLog.txt

	#1, Single- or paired-end
	#2, Reads separated by commas
	#3, Reference FASTA
	#4, Bootstrap value
	#5, Output directory
	#6, User initials

	echo -e "\n"
	echo -e "##################################################################################################################"
	echo -e "RapidRNAseq - Kallisto wrapper, developed by Kevin Blighe (k.blighe@qmul.ac.uk) at Queen Mary University of London"
	echo -e "Release date:\t\tThursday 5th May 2016"
	echo -e "Version:\t\tVersion 2"
	echo -e "Mods to version 2:\tNow supports the analysis of sample replicates concurrently"
	echo -e "Mods to version 2:\tRemoved the necessity to supply output file prefix at the command-line"
	echo -e "##################################################################################################################"
	echo -e "\n"

	echo -e "Analysis log being written to "${5}"/AnalysisLog.txt"
fi

#Check if the first parameter is for single or paired-end sequencing
if [ $(echo "${1}" | awk '{print toupper($0)}') == "SINGLE" ]
then
	#Begin log file
	echo "Beginning script on `date`, run by "${9}" with the following parameters:" > "${7}"/AnalysisLog.txt
	echo -e "\t1\t"${1}"" >> "${7}"/AnalysisLog.txt
	echo -e "\t2\t"${2}"" >> "${7}"/AnalysisLog.txt
	echo -e "\t3\t"${3}"" >> "${7}"/AnalysisLog.txt
	echo -e "\t4\t"${4}"" >> "${7}"/AnalysisLog.txt
	echo -e "\t5\t"${5}"" >> "${7}"/AnalysisLog.txt
	echo -e "\t6\t"${6}"" >> "${7}"/AnalysisLog.txt
	echo -e "\t7\t"${7}"" >> "${7}"/AnalysisLog.txt
	echo -e "\t7\t"${8}"" >> "${7}"/AnalysisLog.txt
	echo -e "\n" >> "${7}"/AnalysisLog.txt

	#1, Single- or paired-end
	#2, Estimated fragment (read) length
	#3, Stadard deviation of estimated fragment (read) length
	#4, Read
	#5, Reference FASTA
	#6, Bootstrap value
	#7, Output directory
	#9, User initials

	echo -e "\n"
	echo -e "##################################################################################################################"
	echo -e "RapidRNAseq - Kallisto wrapper, developed by Kevin Blighe (k.blighe@qmul.ac.uk) at Queen Mary University of London"
	echo -e "Release date:\t\tThursday 5th May 2016"
	echo -e "Version:\t\tVersion 2"
	echo -e "Mods to version 2:\tNow supports the analysis of sample replicates concurrently"
	echo -e "Mods to version 2:\tRemoved the necessity to supply output file prefix at the command-line"
	echo -e "##################################################################################################################"
	echo -e "\n"

	echo -e "Analysis log being written to "${7}"/AnalysisLog.txt"
fi



#############################################
#Analysis step 1:  Index the reference FASTA#
#############################################
#Check if the first parameter is for single or paired-end sequencing
if [ $(echo "${1}" | awk '{print toupper($0)}') == "PAIRED" ]
then
	if [ ! -f "${3}".idx ]
	then
		echo "Beginning analysis step 1 (indexing reference FASTA) on `date`" >> "${5}"/AnalysisLog.txt
		echo `kallisto index -i "${3}".idx "${3}"`
		echo "Finished on `date`" >> "${5}"/AnalysisLog.txt
		echo -e "\n" >> "${5}"/AnalysisLog.txt
	else
		echo "Index file already exists - skipping analysis step 1 (indexing reference FASTA) on `date`" >> "${5}"/AnalysisLog.txt
		echo -e "\n" >> "${5}"/AnalysisLog.txt
	fi
fi

#Check if the first parameter is for single or paired-end sequencing
if [ $(echo "${1}" | awk '{print toupper($0)}') == "SINGLE" ]
then
	if [ ! -f "${5}".idx ]
	then
		echo "Beginning analysis step 1 (indexing reference FASTA) on `date`" >> "${7}"/AnalysisLog.txt
		echo `kallisto index -i "${5}".idx "${5}"`
		echo "Finished on `date`" >> "${7}"/AnalysisLog.txt
		echo -e "\n" >> "${7}"/AnalysisLog.txt
	else
		echo "Index file already exists - skipping analysis step 1 (indexing reference FASTA) on `date`" >> "${7}"/AnalysisLog.txt
		echo -e "\n" >> "${7}"/AnalysisLog.txt
	fi
fi



###########################################
#Analysis step 2:  Quantify read abundance#
###########################################

#Check if the first parameter is for single or paired-end sequencing
if [ $(echo "${1}" | awk '{print toupper($0)}') == "PAIRED" ]
then
	echo "Beginning analysis step 2 (quantifying read abundance) on `date`" >> "${5}"/AnalysisLog.txt

	command="kallisto quant --bias -i "${3}".idx -o "${5}" -b "${4}""

	for ((i=0;i<iNumSequences;i++))
	do
		command="${command}"" ""${arraySequences[i]}"
	done

	echo `$command`

	echo "Finished on `date`" >> "${5}"/AnalysisLog.txt
fi

#Check if the first parameter is for single or paired-end sequencing
if [ $(echo "${1}" | awk '{print toupper($0)}') == "SINGLE" ]
then
	echo "Beginning analysis step 2 (quantifying read abundance) on `date`" >> "${7}"/AnalysisLog.txt

	command="kallisto quant --bias -i "${5}".idx -o "${7}" --single -l "${2}" -s "${3}" -b "${6}""

	for ((i=0;i<iNumSequences;i++))
	do
		command="${command}"" ""${arraySequences[i]}"
	done

	echo `$command`

	echo "Finished on `date`" >> "${7}"/AnalysisLog.txt
fi

exit 0

