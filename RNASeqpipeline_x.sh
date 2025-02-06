#!/bin/bash

#time the script
start_time=`date +%s`

# Display copyright information with ASCII art
echo "########################################"
echo "#                                      #"
echo "#  C O P Y R I G H T   N O T I C E    #"
echo "#                                      #"
echo "########################################"
echo ""
echo "  Copyright © 2025 Jaidev.              "
echo "  Licensed under CC BY-NC 4.0.          "
echo ""
echo "########################################"


echo '         ___'
echo '       _(((,|			How do I analyze RNA-seq data?'
echo '      /  _-\\'
echo '     / C o\o \'
echo '   _/_    __\ \     __ __     __ __     __ __     __'
echo '  /   \ \___/  )   /--X--\   /--X--\   /--X--\   /--/'
echo '  |    |\_|\  /   /--/ \--\ /--/ \--\ /--/ \--\ /--/'
echo '  |    |#  #|/          \__X__/   \__X__/   \__X__/'
echo '  (   /     |'
echo '   |  |#  # |'
echo '   |  |    #|'
echo '   |  | #___n_,_'
echo ',-/   7-'"'"' .     `\'
echo '`-\...\\-_   -  o /'
echo '   |#  # `---U--'
echo '   `-v-^-'"'"'/'
echo '     \  |_|_ Wny'
echo '     (___mnnm'


printf "
Welcome to the RNA-seq analysis script!\n
This script will help you process and analyze your RNA-seq data.\n
Please make sure that you have the necessary software and dependencies installed before running this script.\n
1. fastp: (conda install -c bioconda fastp)\n
2. HISAT2: (conda install -c bioconda hisat2)\n
3. samtools: (conda install -c bioconda samtools)\n
5. trimmomatic : (conda install -c bioconda trimmomatic)\n
6. fastqc: (conda install -c bioconda fastqc)\n
7. featureCounts: (conda install -c bioconda subread)\n\n\n
8. Anaconda or Miniconda (script install)
**** Stable internet connection to fetch required files and sufficient space and memory on your system
***** Fastq files are assumed to end with '_R1_001.fastq.gz' or '_R2_001.fastq.gz' to describe paired-end reads\n\n\n
"

declare -a dependencies=("fastp" "trimmomatic" "hisat2" "samtools" "trimmomatic" "fastqc" "featureCounts" "multiqc")
missing_dependencies=()

# Function to install Conda
install_conda() {
    echo "Installing Miniconda..."
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    chmod +x miniconda.sh
    ./miniconda.sh -b -p $HOME/miniconda
    export PATH="$HOME/miniconda/bin:$PATH"
    rm miniconda.sh
    echo "Miniconda installed successfully."
}

# Check for Conda
if ! command -v conda >/dev/null 2>&1; then
    missing_dependencies+=("conda")
fi

# Prompt the user to install Conda
if [[ " ${missing_dependencies[@]} " =~ " conda " ]]; then
    read -p "Conda is not installed. Do you want to install Conda? [y/n] " answer
    if [[ $answer = y ]]; then
        install_conda
    else
        echo "Please install Conda manually before proceeding."
    fi
fi

# Check for other dependencies
for dependency in "${dependencies[@]}"; do
    if ! command -v "$dependency" >/dev/null 2>&1; then
        missing_dependencies+=("$dependency")
    fi
done

# Prompt the user to install missing dependencies
if [[ ${#missing_dependencies[@]} -gt 0 ]]; then
    echo "The following dependencies are missing: ${missing_dependencies[*]}"
    read -p "Do you want to install the missing dependencies? [y/n] " answer

    if [[ $answer = y ]]; then
        for dependency in "${missing_dependencies[@]}"; do
            case $dependency in
            "fastp")
                conda install -c bioconda fastp
                ;;
            "hisat2")
                conda install -c bioconda hisat2
                ;;
            "samtools")
                conda install -c bioconda samtools
                ;;
            "trimmomatic")
                conda install -c bioconda trimmomatic
                ;;
            "fastqc")
                conda install -c bioconda fastqc
                ;;
            "featureCounts")
                conda install -c bioconda subread
                ;;
            "conda")
                # This case should not be reached as Conda was checked earlier
                echo "Conda is not installed. Please install Conda manually."
                ;;
            esac
        done

        echo "Dependencies installed successfully!"
    else
        echo "Please install the missing dependencies manually before running the script."
    fi
else
    echo "All dependencies are installed."
fi


# Store the current working directory path
working_dir=$(realpath "$(pwd)")
echo "Working Directory: $working_dir"

# Prompt the user for the number of threads
echo "Enter the number of threads:"
read num_threads

# Check if the input is a valid positive integer
if ! [[ "$num_threads" =~ ^[1-9][0-9]*$ ]]; then
  echo "Invalid input. Please enter a positive integer."
  exit 1
fi

# Use the input in your script
echo "Running with $num_threads threads..."

# Ask user to select organism
printf "Select the organism for which to download the genome indices:\n 1. Rattus norvegicus\n 2. Homo sapiens (default)\n 3. Mus musculus\n 4. Drosophila melanogaster\n 5. C. elegans\n 6. S. cerevisiae\n"
read -p "Enter your choice (1-6, default is 2): " user_input
user_input=${user_input:-2}  # Default to 2 if no input is provided

# Predefined default paths for each organism
default_human_index="your_path_to/index/hisat2/hg38/grch38"
default_rattus_index="your_path_to/index/hisat2/rn6"
default_mus_index="your_path_to/index/hisat2/grcm38"
default_drosophila_index="your_path_to/index/hisat2/dm6"
default_celegans_index="your_path_to/index/hisat2/celegans"
default_scer_index="your_path_to/index/hisat2/scer"

# Define default paths for each organism
default_human_gtf="your_path_to/human/Homo_sapiens.GRCh38.108.gtf"
default_rattus_gtf="your_path_to/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.108.chr.gtf"
default_mus_gtf="your_path_to/mus_musculus/Mus_musculus.GRCm39.109.gtf"
default_drosophila_gtf="your_path_to/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.32.109.gtf"
default_celegans_gtf="your_path_to/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.109.gtf"
default_scerevisiae_gtf="your_path_to/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.109.gtf"

# Set default paths based on organism selection
case $user_input in
    1)
        default_index="$default_rattus_index"
        default_gtf="$default_rattus_gtf"
        organism_name="Rattus norvegicus"
        ;;
    2)
        default_index="$default_human_index"
        default_gtf="$default_human_gtf"
        organism_name="Homo sapiens"
        ;;
    3)
        default_index="$default_mus_index"
        default_gtf="$default_mus_gtf"
        organism_name="Mus musculus"
        ;;
    4)
        default_index="$default_drosophila_index"
        default_gtf="$default_drosophila_gtf"
        organism_name="Drosophila melanogaster"
        ;;
    5)
        default_index="$default_celegans_index"
        default_gtf="$default_celegans_gtf"
        organism_name="C. elegans"
        ;;
    6)
        default_index="$default_scer_index"
        default_gtf="$default_scerevisiae_gtf"
        organism_name="S. cerevisiae"
        ;;
    *)
        echo "Invalid choice"
        exit 1
        ;;
esac

echo "Selected organism: $organism_name"

# Prompt user to provide a path for the index or choose to download the index files
echo "Would you like to provide a pre-downloaded index path or download the index files?"
echo "1. Provide path to existing HISAT2 index"
echo "2. Download HISAT2 index files"
read -p "Enter your choice (default is 1): " choice
choice=${choice:-1}  # Default to 1 if no input is provided

if [ "$choice" == "1" ]; then
    # Prompt the user for the HISAT2 index directory path
    read -p "Enter the absolute path to your HISAT2 index directory (default is $default_index): " user_provided_path
    
    if [ -z "$user_provided_path" ]; then
        # Use the default path if the user doesn't provide one
        path_index="$default_index"
        echo "Using default HISAT2 index path: $path_index"
    elif [ -d "$user_provided_path" ]; then
        # Validate and use the user-provided path
        path_index="$user_provided_path"
        echo "Using provided HISAT2 index path: $path_index"
    else
        # Handle invalid directory
        echo "Error: Path does not exist or is not a directory. Exiting."
        exit 1
    fi
elif [ "$choice" == "2" ]; then
    # Handle the preset HISAT2 index case
    echo "Using preset HISAT2 index (will download index if necessary)."
    path_index="$default_index"
    # Add your download logic here
else
    # Handle invalid choices
    echo "Error: Invalid choice. Please select a valid option."
    exit 1
fi

# Ask the user whether they want to provide a path or download the GTF file
echo "Do you want to provide a custom path to the GTF file or download it?"
echo "1. Provide path to existing GTF file"
echo "2. Download GTF file"
read -p "Enter your choice (default is 2): " path_choice
path_choice=${path_choice:-2}  # Default to download if no input is provided

if [ "$path_choice" -eq 1 ]; then 
    read -p "Please provide the absolute path to your GTF file (default is $default_gtf): " user_provided_path
    if [ -z "$user_provided_path" ]; then
        # Use default path if no input provided
        path_gtf="$default_gtf"
        echo "Using default GTF file: $path_gtf"
    elif [ -f "$user_provided_path" ]; then
        path_gtf="$user_provided_path"
        echo "Using provided GTF file: $path_gtf"
    else
        echo "Error: GTF file does not exist"
        exit 1
    fi
else
    path_gtf="$default_gtf"
    echo "Using default GTF file: $path_gtf"
    # Add your download logic here if the file doesn't exist
fi

# Print final configuration
echo -e "\nFinal Configuration:"
echo "Organism: $organism_name"
echo "Number of threads: $num_threads"
echo "Index path: $path_index"
echo "GTF file: $path_gtf"

#print the working directory
echo "Current directory: $(pwd)"
echo "fastq files in the working directory"
ls *.fastq.gz
#promt for trimming
# STEP 1: QC and trimming
echo "*************************************************************************************************************************************************************** "
printf "STEP 1: QC and trimming \n\n perform QC and trimming on fastq files\n\n\n\n\n"
# Prompt the user to perform QC and trimming on fastq files
read -p "Do you want to perform QC & trimming of the fastq files? [y/n]" answer

# If the user selects "y", ask them to choose between two methods: fastp or fastqc
	if [[ $answer = y ]]; then
		read -p "Please select the method for analysis: 1) fastp or 2) fastqc? " answer
		if [[ $answer = 1 ]]; then
			echo "Performing QC & trimming of the fastq files using fastp..."

			# Use fastp to trim the reads and output trimmed files and QC reports
			# Input files should be in the current directory and end with "_R1_001.fastq.gz" and "_R2_001.fastq.gz"
			# Output files will be named "trimmed_[filename]" and placed in a subdirectory named "processed_files"
			# QC reports will be named "[filename]-fastp.html" and "[filename]-fastp.json" and placed in a subdirectory named "fastp_results"
			read_file1=($(ls -d *_R1_001.fastq.gz))
			read_file2=($(ls -d *_R2_001.fastq.gz))
			
			echo "trimmed files will be saved to processed_files and QC results to fastp_results"
			mkdir processed_files
			mkdir fastp_results
			
			for f in "${read_file1[@]}"; do
				# Construct the name of the corresponding file in the read_file2 array
				f2=${f/R1_001.fastq.gz/R2_001.fastq.gz}
				f3=${f/R1_001.fastq.gz/R1_R2.fastq.gz}
				echo "Processing $f and $f2..."

				# Run fastp to trim the reads and generate QC reports
				fastp -i $f -I $f2 -o processed_files/trimmed_${f} -O processed_files/trimmed_${f2} -h fastp_results/${f3}-fastp.html -j fastp_results/${f3}-fastp.json
				#fastp -i $f -I $f2 -o trimmed_$f -O trimmed_$f2 -h $f3-fastp.html -j $f3-fastp.json
			done

			echo "QC and trimming completed using fastp!"
							
			echo "Current directory: $(pwd)"
			cd processed_files
			ls
				
				
		elif [[ $answer = 2 ]]; then
			
			echo "trimmed files will be saved to processed_files and QC results to fastqc_results"
			mkdir processed_files
			mkdir fastqc_results
			
			echo "Performing QC & trimming of the fastq files using fastqc and Trimmomatic..."

			# Use fastqc to generate QC reports
			# Input files should be in the current directory and end with ".fastq.gz"
			# QC reports will be named "[filename]_fastqc.html" and placed in the current directory
			
			fastqc *.fastq.gz -o . -t $num_threads
			#mv *.fastqc.html *.fastqc.zip ./fastqc_results &&
				mv *_fastqc.html ./fastqc_results
				mv *_fastqc.zip ./fastqc_results
		
			# Use Trimmomatic to trim the reads with poor quality
			# Input files should be in the current directory and end with ".fastq.gz"
			# Output files will be named "[filename]_trimmed.fastq" and placed in the current directory
			# Trim reads with a trailing quality score of less than 10 (TRAILING:10) and use Phred33 scoring (phred33)
			read -p "Do you want to perform trimming of the fastq files? [y/n]" ans_usr
			if [[ $ans_usr = y ]]; then
				wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
				unzip Trimmomatic-0.39.zip
				java -jar Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 8 *.fastq.gz processed_files/*.fastq_trimmed TRAILING:10 -phred33
				mv *.fastq_trimmed processed_files
				cd processed_files
				ls
    echo "Trimming completed using Trimmomatic!"
			fi	
			echo "QC and trimming aborted!"
			
			
		else
			echo "Invalid input. Please select 1 for fastp or 2 for fastqc."
		fi
	else
		echo "QC and trimming aborted!"
	fi


# STEP 2: Run HISAT2
pwd
echo "***************************************************************************************************************************************************************"
printf "STEP 2: Run HISAT2 \n\n Alignment using HISAT2 \n\n\n\n\n"

# Assign organism and download URL based on input
case $user_input in
    1) 
        organism="Rattus norvegicus" 
        download_url="https://genome-idx.s3.amazonaws.com/hisat/rn6_genome.tar.gz"
        default_index="$default_rattus_index"
        ;;
    2|"")
        organism="Homo sapiens (default)" 
        download_url="https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz"
        default_index="$default_human_index"
        ;;
    3) 
        organism="Mus musculus" 
        download_url="https://cloud.biohpc.swmed.edu/index.php/s/grcm38/download"
        default_index="$default_mus_index"
        ;;
    4) 
        organism="Drosophila melanogaster" 
        download_url="https://genome-idx.s3.amazonaws.com/hisat/dm6.tar.gz"
        default_index="$default_drosophila_index"
        ;;
    5) 
        organism="C. elegans" 
        download_url="https://cloud.biohpc.swmed.edu/index.php/s/bbynxoY2TPpRNQb/download"
        default_index="$default_celegans_index"
        ;;
    6) 
        organism="S. cerevisiae" 
        download_url="https://cloud.biohpc.swmed.edu/index.php/s/Gsq4goLW4TDAz4E/download"
        default_index="$default_scer_index"
        ;;
    *)
        echo "Invalid input! Defaulting to Homo sapiens."
        organism="Homo sapiens (default)"
        download_url="https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz"
        default_index="$default_human_index"
        ;;
esac

echo "Selected organism: $organism"

# Decide whether to use the default path or ask the user for a custom path
if [ "$choice" == "1" ]; then

	echo "Using provided HISAT2 index path: $path_index"
	
elif [ "$choice" == "2" ]; then
    # If the user chose to download the index
    echo "Creating HISAT2 index directory..."
    mkdir -p hisat2-index
    cd hisat2-index || exit
    echo "Current directory: $(pwd)"
    
    echo "Downloading index files for $organism..."
    wget "$download_url" -O genome_index.tar.gz
    tar -xvf genome_index.tar.gz
    path_index="$(pwd)/$(ls -d */ | head -n 1)" # Get the first directory created
    echo "Genome index downloaded and extracted to: $path_index"
    cd ../..
else
    echo "Invalid choice. Exiting."
    exit 1
fi

echo "Genome index is set to: $path_index"


# run alignment
echo "running alignment in :$(pwd)"
read_set1=($(ls -d *R1_001.fastq.gz))
read_set2=($(ls -d *R2_001.fastq.gz))

# a = fwd_read , a2 = bkw_read, a3 = output_file
for a in "${read_set1[@]}"; do
    a2=${a/R1_001.fastq.gz/R2_001.fastq.gz}
    a3=${a/R1_001.fastq.gz/R1_R2}
    echo "Now processing: $a3"
    # hisat2 [options]* -x <ht2-idx> {-1 <m1> -2 <m2> | -U <r>} [-S <sam>]
    hisat2 -p "$num_threads" --dta -x $path_index/genome -1 $a -2 $a2 -S $a3.sam --summary-file $a3.summary.txt
    samtools sort $a3.sam -o $a3.bam
    echo "removing $a3.sam to save space"
    rm $a3.sam
    echo "HISAT2 finished running for $a3!" 
	printf "________________________________________________________________________________________________\n"
done

# STEP 3: Run featureCounts - Quantification
echo "***************************************************************************************************************************************************************"
printf "STEP 3: Run featureCounts - Quantification \n\n Get Counts Matrix \n\n\n\n\n"
## Late Script: Download or process the GTF file based on user input

# Check if the user wants to provide a custom path or download the GTF file
if [ "$path_choice" == "1" ]; then
    # If the user chooses to provide a custom path, ask for it
    
    if [ -f "$user_provided_path" ]; then
        echo "Using provided GTF file: $user_provided_path"
        gtf_file="$user_provided_path"
        path_gtf=$(dirname "$gtf_file")  # Set path to the directory of the GTF file
    else
        echo "Error: The provided path is invalid or the file does not exist. Using default GTF file."
        # Use the default GTF file for the chosen organism
        case $user_input in
            1) gtf_file=${default_rattus_gtf} ;;
            2) gtf_file=${default_human_gtf} ;;
            3) gtf_file=${default_mus_gtf} ;;
            4) gtf_file=${default_drosophila_gtf} ;;
            5) gtf_file=${default_celegans_gtf} ;;
            6) gtf_file=${default_scerevisiae_gtf} ;;
            *)
                echo "Invalid input. Exiting."
                exit 1
                ;;
        esac
        echo "Using default GTF file: $gtf_file"
    fi
elif [ "$path_choice" == "2" ]; then
    # If the user chooses to download the GTF file, set up the directory for genome annotation
    mkdir genome_annotation
    cd genome_annotation
    path_gtf=($(pwd))  # Store current path

    # Download the appropriate GTF file based on the user's organism choice
    case $user_input in
        1)
            # Rattus norvegicus
            wget https://ftp.ensembl.org/pub/release-108/gtf/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.108.chr.gtf.gz
            gunzip Rattus_norvegicus.mRatBN7.2.108.chr.gtf.gz
            gtf_file="Rattus_norvegicus.mRatBN7.2.108.chr.gtf"
            ;;
        2)
            # Homo sapiens
            wget https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.gtf.gz
            gunzip Homo_sapiens.GRCh38.108.gtf.gz
            gtf_file="Homo_sapiens.GRCh38.108.gtf"
            ;;
        3)
            # Mus musculus
            wget https://ftp.ensembl.org/pub/release-109/gtf/mus_musculus/Mus_musculus.GRCm39.109.gtf.gz
            gunzip Mus_musculus.GRCm39.109.gtf.gz
            gtf_file="Mus_musculus.GRCm39.109.gtf"
            ;;
        4)
            # Drosophila melanogaster
            wget https://ftp.ensembl.org/pub/release-109/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.32.109.gtf.gz
            gunzip Drosophila_melanogaster.BDGP6.32.109.gtf.gz
            gtf_file="Drosophila_melanogaster.BDGP6.32.109.gtf"
            ;;
        5)
            # C. elegans
            wget https://ftp.ensembl.org/pub/release-109/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.109.gtf.gz
            gunzip Caenorhabditis_elegans.WBcel235.109.gtf.gz
            gtf_file="Caenorhabditis_elegans.WBcel235.109.gtf"
            ;;
        6)
            # S. cerevisiae
            wget https://ftp.ensembl.org/pub/release-109/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.109.gtf.gz
            gunzip Saccharomyces_cerevisiae.R64-1-1.109.gtf.gz
            gtf_file="Saccharomyces_cerevisiae.R64-1-1.109.gtf"
            ;;
        *)
            echo "Invalid choice. Exiting."
            exit 1
            ;;
    esac

    # Now set the path to the downloaded GTF file
    path_gtf=($(pwd))/$(basename $gtf_file)  # Append the GTF file name to the current directory path

    # Extracted GTF file path
    echo "Your Genome index: $gtf_file ... downloaded and extracted to: $path_gtf"
else
    echo "Invalid choice. Exiting."
    exit 1
fi

# Move to the parent directory to create the 'quants' folder
cd "$working_dir/.." || exit 1
mkdir -p quants

# 🔍 Detect folder containing BAM files
if compgen -G "$working_dir/processed_files/*.bam" > /dev/null; then
    bam_dir="$working_dir/processed_files"
    echo "BAM files found in: $bam_dir"
elif compgen -G "$working_dir/*.bam" > /dev/null; then
    bam_dir="$working_dir"
    echo "BAM files found in: $bam_dir"
else
    echo " No BAM files found in 'processed_files' or the parent directory."
    exit 1
fi

# Path to the GTF file (assuming it's already defined in your previous steps)
path_gtf=$(realpath "$gtf_file")
echo "Using GTF file: $path_gtf"

# Run featureCounts
featureCounts -p -T 60 -a "$path_gtf" -o quants/all-featurecounts.txt "$bam_dir"/*.bam

# Run MultiQC (optional)
echo "Generating MultiQC report..."
# multiqc .

# Calculate runtime
end_time=$(date +%s)
runtime=$((end_time - start_time))
printf "This script took %s seconds\n" "$runtime"

#run multiqc
echo "genrating html report"
multiqc .

end_time=`date +%s`
runtime=$((end_time-start_time))
printf "this script took %s seconds\n" "$runtime"