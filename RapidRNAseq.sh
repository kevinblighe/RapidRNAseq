#!/bin/bash

###########################################################
#Pre-liminary step 1:  Safeguarding the script's operation#
###########################################################

#Check if 0 parameters have been passed
if [ $# -eq 0 ]
then
	echo -e "\n"
	echo "Error - incorrect number of parameters.  Consult the relevant standard operating procedure for correct usage."
	echo "Please type RapidRNAseq.sh -h for help."
	echo -e "\n"
	exit 1
fi

#Check if the help parameter was passed
if [ "${1}" == "-h" ]
then
	echo -e "##################################################################################################################"
	echo -e "RapidRNAseq - Kallisto wrapper, developed by Kevin Blighe (k.blighe@qmul.ac.uk) at Queen Mary University of London"
	echo -e "Release date:\tThursday 5th May 2016"
	echo -e "Version:\tVersion 1"
	echo -e "##################################################################################################################"
	echo -e "\n"

	echo -e "Program syntax:"
	echo -e "Paired-end"
	echo -e "RapidRNAseq.sh PAIRED [FASTQ/FASTQ.gz mate-pair 1] [FASTQ/FASTQ.gz mate-pair 2] [Reference FASTA/FASTA.gz] [Bootstrap quantification value] [Output dir] [Output prefix] [User initials]"
	echo -e "\n"

	echo -e "Single-end"
	echo -e "RapidRNAseq.sh SINGLE [Fragment (read) length] [Fragment (read) length standard deviation] [FASTQ/FASTQ.gz] [Reference FASTA/FASTA.gz] [Bootstrap quantification value] [Output dir] [Output prefix] [User initials]"
	echo -e "\n"

	exit 0
fi

#Check that no less than 8 parameters have been passed - 8 for paired is minimum possible needed (9 for single)
if [ $# -lt 8 ]
then
	echo -e "\n"
	echo "Error - incorrect number of parameters.  Consult the relevant standard operating procedure for correct usage."
	echo "Please type RapidRNAseq.sh -h for help."
	echo -e "\n"
	exit 1
fi

#Check if the first parameter is for single or paired-end sequencing
if [ $(echo "${1}" | awk '{print toupper($0)}') == "PAIRED" ]
then
	#Check that 8 parameters have been passed
	if [ $# -ne 8 ]
	then
		echo -e "\n"
		echo "Error - incorrect number of parameters for paired-end analysis.  Consult the relevant standard operating procedure for correct usage."
		echo "Please type RapidRNAseq.sh -h for help."
		echo -e "\n"
		exit 1
	fi

	#Check that the file as $2 exists
	test -e "${2}"
	if [ $? -ne 0 ]
	then
		echo -e "\n"
		echo "Error - input file ${2} does not exist.  Please check the complete file path and re-run."
		echo "Note:  input sample files can be in FASTQ or FASTQ.gz format."
		echo "Please type RapidRNAseq.sh -h for help."
		echo -e "\n"
		exit 1
	fi

	#Check that the file as $3 exists
	test -e "${3}"
	if [ $? -ne 0 ]
	then
		echo -e "\n"
		echo "Error - input file ${3} does not exist.  Please check the complete file path and re-run."
		echo "Note:  input sample files can be in FASTQ or FASTQ.gz format."
		echo "Please type RapidRNAseq.sh -h for help."
		echo -e "\n"
		exit 1
	fi

	#Check that the file as $4 exists
	test -e "${4}"
	if [ $? -ne 0 ]
	then
		echo -e "\n"
		echo "Error - reference file ${4} does not exist.  Please check the complete file path and re-run."
		echo "Note:  reference sequences can be in FASTA or FASTA.gz format"
		echo "Please type RapidRNAseq.sh -h for help."
		echo -e "\n"
		exit 1
	fi

	#Check that $5 is numerical
	if ! [ "${5}" -eq "${5}" ] 2>/dev/null
	then
		echo -e "\n"
		echo "Error - bootstrap value ${5} must be an integer value."
		echo "Please type RapidRNAseq.sh -h for help."
		echo -e "\n"
		exit 1
	fi

	#Check that $6 is not a file
	if [ -f "${6}" ]
	then
		echo -e "\n"
		echo "Error - output directory ${6} already exists and is a regular file."
		echo "Please choose another output directory"
		echo -e "\n"
		exit 1
	fi

	#Check that $6 is a directory
	if [ ! -d "${6}" ]
	then
		echo -e "\n"
		echo "Output directory ${6} does not exist."
		echo "Creating..."
		mkdir "${6}"
	fi
fi

#Check if the first parameter is for single or paired-end sequencing
if [ $(echo "${1}" | awk '{print toupper($0)}') == "SINGLE" ]
then
	#Check that 9 parameters have been passed
	if [ $# -ne 9 ]
	then
		echo -e "\n"
		echo "Error - incorrect number of parameters for single-end analysis.  Consult the relevant standard operating procedure for correct usage."
		echo "Please type RapidRNAseq.sh -h for help."
		echo -e "\n"
		exit 1
	fi

	#Check that $2 is numerical
	if ! [ "${2}" -eq "${2}" ] 2>/dev/null
	then
		echo -e "\n"
		echo "Error - estimated fragment (read) length ${2} must be an integer value"
		echo "Please type RapidRNAseq.sh -h for help."
		echo -e "\n"
		exit 1
	fi

	#Check that $3 is numerical
	if ! [ "${3}" -eq "${3}" ] 2>/dev/null
	then
		echo -e "\n"
		echo "Error - fragment length standard deviation ${3} must be an integer value"
		echo "Please type RapidRNAseq.sh -h for help."
		echo -e "\n"
		exit 1
	fi

	#Check that the file as $4 exists
	test -e "${4}"
	if [ $? -ne 0 ]
	then
		echo -e "\n"
		echo "Error - input file ${4} does not exist.  Please check the complete file path and re-run."
		echo "Note:  input sample files can be in FASTQ or FASTQ.gz format"
		echo "Please type RapidRNAseq.sh -h for help."
		echo -e "\n"
		exit 1
	fi

	#Check that the file as $5 exists
	test -e "${5}"
	if [ $? -ne 0 ]
	then
		echo -e "\n"
		echo "Error - reference file ${5} does not exist.  Please check the complete file path and re-run."
		echo "Note:  reference sequences can be in FASTA or FASTA.gz format"
		echo "Please type RapidRNAseq.sh -h for help."
		echo -e "\n"
		exit 1
	fi

	#Check that $6 is numerical
	if ! [ "${6}" -eq "${6}" ] 2>/dev/null
	then
		echo -e "\n"
		echo "Error - bootstrap value ${6} must be an integer value"
		echo "Please type RapidRNAseq.sh -h for help."
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
		echo "Please type RapidRNAseq.sh -h for help."
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
	echo "Beginning script on `date`, run by "${8}" with the following parameters:" > "${6}"/"${7}"_AnalysisLog.txt
	echo -e "\t1\t"${1}"" >> "${6}"/"${7}"_AnalysisLog.txt
	echo -e "\t2\t"${2}"" >> "${6}"/"${7}"_AnalysisLog.txt
	echo -e "\t3\t"${3}"" >> "${6}"/"${7}"_AnalysisLog.txt
	echo -e "\t4\t"${4}"" >> "${6}"/"${7}"_AnalysisLog.txt
	echo -e "\t5\t"${5}"" >> "${6}"/"${7}"_AnalysisLog.txt
	echo -e "\t6\t"${6}"" >> "${6}"/"${7}"_AnalysisLog.txt
	echo -e "\t7\t"${7}"" >> "${6}"/"${7}"_AnalysisLog.txt
	echo -e "\t8\t"${8}"" >> "${6}"/"${7}"_AnalysisLog.txt
	echo -e "\n" >> "${6}"/"${7}"_AnalysisLog.txt

	#1, Single- or paired-end
	#2, Read1
	#3, Read2
	#4, Reference FASTA
	#5, Bootstrap value
	#6, Output directory
	#7, Prefix for output
	#8, User initials

	echo -e "\n"
	echo -e "##################################################################################################################"
	echo -e "RapidRNAseq - Kallisto wrapper, developed by Kevin Blighe (k.blighe@qmul.ac.uk) at Queen Mary University of London"
	echo -e "Release date:\tThursday 5th May 2016"
	echo -e "Version:\tVersion 1"
	echo -e "##################################################################################################################"
	echo -e "\n"

	echo -e "Analysis log being written to "${6}"/"${7}"_AnalysisLog.txt"
fi

#Check if the first parameter is for single or paired-end sequencing
if [ $(echo "${1}" | awk '{print toupper($0)}') == "SINGLE" ]
then
	#Begin log file
	echo "Beginning script on `date`, run by "${9}" with the following parameters:" > "${7}"/"${8}"_AnalysisLog.txt
	echo -e "\t1\t"${1}"" >> "${7}"/"${8}"_AnalysisLog.txt
	echo -e "\t2\t"${2}"" >> "${7}"/"${8}"_AnalysisLog.txt
	echo -e "\t3\t"${3}"" >> "${7}"/"${8}"_AnalysisLog.txt
	echo -e "\t4\t"${4}"" >> "${7}"/"${8}"_AnalysisLog.txt
	echo -e "\t5\t"${5}"" >> "${7}"/"${8}"_AnalysisLog.txt
	echo -e "\t6\t"${6}"" >> "${7}"/"${8}"_AnalysisLog.txt
	echo -e "\t7\t"${7}"" >> "${7}"/"${8}"_AnalysisLog.txt
	echo -e "\t8\t"${8}"" >> "${7}"/"${8}"_AnalysisLog.txt
	echo -e "\t9\t"${9}"" >> "${7}"/"${8}"_AnalysisLog.txt
	echo -e "\n" >> "${7}"/"${8}"_AnalysisLog.txt

	#1, Single- or paired-end
	#2, Estimated fragment (read) length
	#3, Stadard deviation of estimated fragment (read) length
	#4, Read
	#5, Reference FASTA
	#6, Bootstrap value
	#7, Output directory
	#8, Prefix for output
	#9, User initials

	echo -e "\n"
	echo -e "##################################################################################################################"
	echo -e "RapidRNAseq - Kallisto wrapper, developed by Kevin Blighe (k.blighe@qmul.ac.uk) at Queen Mary University of London"
	echo -e "Release date:\tThursday 5th May 2016"
	echo -e "Version:\tVersion 1"
	echo -e "##################################################################################################################"
	echo -e "\n"

	echo -e "Analysis log being written to "${7}"/"${8}"_AnalysisLog.txt"
fi



#############################################
#Analysis step 1:  Index the reference FASTA#
#############################################
#Check if the first parameter is for single or paired-end sequencing
if [ $(echo "${1}" | awk '{print toupper($0)}') == "PAIRED" ]
then
	if [ ! -f "${4}".idx ]
	then
		echo "Beginning analysis step 1 (indexing reference FASTA) on `date`" >> "${6}"/"${7}"_AnalysisLog.txt
		echo `kallisto index -i "${4}".idx "${4}"`
		echo "Finished on `date`" >> "${6}"/"${7}"_AnalysisLog.txt
		echo -e "\n" >> "${6}"/"${7}"_AnalysisLog.txt
	else
		echo "Index file already exists - skipping analysis step 1 (indexing reference FASTA) on `date`" >> "${6}"/"${7}"_AnalysisLog.txt
		echo -e "\n" >> "${6}"/"${7}"_AnalysisLog.txt
	fi
fi

#Check if the first parameter is for single or paired-end sequencing
if [ $(echo "${1}" | awk '{print toupper($0)}') == "SINGLE" ]
then
	if [ ! -f "${5}".idx ]
	then
		echo "Beginning analysis step 1 (indexing reference FASTA) on `date`" >> "${7}"/"${8}"_AnalysisLog.txt
		echo `kallisto index -i "${5}".idx "${5}"`
		echo "Finished on `date`" >> "${7}"/"${8}"_AnalysisLog.txt
		echo -e "\n" >> "${7}"/"${8}"_AnalysisLog.txt
	else
		echo "Index file already exists - skipping analysis step 1 (indexing reference FASTA) on `date`" >> "${7}"/"${8}"_AnalysisLog.txt
		echo -e "\n" >> "${7}"/"${8}"_AnalysisLog.txt
	fi
fi



###########################################
#Analysis step 2:  Quantify read abundance#
###########################################

#Check if the first parameter is for single or paired-end sequencing
if [ $(echo "${1}" | awk '{print toupper($0)}') == "PAIRED" ]
then
	echo "Beginning analysis step 2 (quantifying read abundance) on `date`" >> "${6}"/"${7}"_AnalysisLog.txt

	echo `kallisto quant --bias -i "${4}".idx -o "${6}" -b "${5}" "${2}" "${3}"`

	echo "Finished on `date`" >> "${6}"/"${7}"_AnalysisLog.txt
fi

#Check if the first parameter is for single or paired-end sequencing
if [ $(echo "${1}" | awk '{print toupper($0)}') == "SINGLE" ]
then
	echo "Beginning analysis step 2 (quantifying read abundance) on `date`" >> "${7}"/"${8}"_AnalysisLog.txt

	echo `kallisto quant --bias -i "${5}".idx -o "${7}" --single -l "${2}" -s "${3}" -b "${6}" "${4}"`

	echo "Finished on `date`" >> "${7}"/"${8}"_AnalysisLog.txt
fi

exit 0
