Pour installer linux en local : 
wsl 
https://learn.microsoft.com/fr-fr/windows/wsl/setup/environment
Dans le powershell de Windows : wsl --install

#Installer conda 
https://docs.conda.io/projects/miniconda/en/latest/   
chmod +x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh 
conda --version

#installer cutadapt
pip install cutadapt


#1/ #Utiliser Fastqc pour connaitre la qualité des données 
Téléchargement : https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
unzip fastqc_vX.XX.zip
#Remplacez fastqc_vX.XX.zip par le nom du fichier que vous avez téléchargé.

#Accès au répertoire FastQC : Naviguez vers le répertoire FastQC nouvellement créé :

cd FastQC/

#Autorisation d'exécution : Avant d'exécuter FastQC, vous devrez peut-être accorder les autorisations d'exécution au script fastqc :

chmod +x fastqc

#Exécution de FastQC :

#Une fois l'installation terminée, vous pouvez exécuter FastQC sur un fichier FASTQ en utilisant la commande suivante dans le terminal :



																		#2/Pour compter le nombre de read en paired-end à chaque étape
zcat file_R1.fastq.gz | echo $((`wc -l`/4))
zcat file_R2.fastq.gz | echo $((`wc -l`/4))
																		 
																		#5/ Pour trimmer en enelvant les adaptateurs, les séquences de mauvaises qualités et trop courtes 
cutadapt -a ADAPTER_SEQUENCE_R1 -A ADAPTER_SEQUENCE_R2 -q QUALITY_THRESHOLD --minimum-length MIN_LENGTH -o trimmed_sample1_R1.fastq.gz -p trimmed_sample1_R2.fastq.gz file_R1.fastq.gz file_R2.fastq.gz > trimmed_sample1_cutadapt.log
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -q 30 -o trimmed_PG.R1.fastq -p trimmed_PG.R2.fastq PG.R1.fastq.gz PG.R2.fastq.gz > trimmed_sample1_cutadapt.log

																		#6/Pour compter le nombre de read après démultiplexage par génotype 
																		#7Utiliser Fastqc pour connaitre la qualité des données 
																		





############################ BWA MEM #############################
bwa index ref_betae.fasta #On indexe la reference
bwa mem ref_betae.fasta trimmed_PG.R1.fastq trimmed_PG.R2.fastq > map_PG_betae.sam #on map tout les read R1 et R2 sur la reference
samtools view -bS map_PG_betae.sam > map_PG_betae.bam #On les transforme en fichier bam
samtools sort -o map_PG_betae.sorted.bam map_PG_betae.bam #On trie les fichiers bam
samtools rmdup map_PG_betae.sorted.bam final.PG_betae.bam #o enleve les duplicats de PCR
samtools index final.PG_betae.bam #On les index

############################ BWA MEM POUR MAPER LE rDNA WHEAT #############################
bwa index Human_rDNA.fasta #On indexe la reference
bwa mem  Human_rDNA.fasta trimmed_PG.R1.fastq trimmed_PG.R2.fastq > map_PG_human.sam #on map tout les read R1 et R2 sur la reference
samtools view -bS map_PG_human.sam > map_PG_human.bam #On les transforme en fichier bam
samtools sort -o map_PG_human.sorted.bam map_PG_human.bam #On trie les fichiers bam
samtools rmdup map_PG_human.sorted.bam final.PG_human.bam #o enleve les duplicats de PCR
samtools index final.PG_human.bam #On les index

############################ BWA MEM POUR MAPER LE rDNA PG et rDNA WHEAT de maniere competitive #############################
cat ref_polymyxa.fasta sequence_rDNA_wheat.fasta rDNA_whaeat_2.fasta rDNA18s_83pb.fasta > combined_references.fasta  #on concatene les ref dans un fichier fasta

bwa index combined_references.fasta

bwa mem combined_references.fasta trimmed_PG.R1.fastq trimmed_PG.R2.fastq > map_PG_no_wheat.sam #on map tout les read R1 et R2 sur la reference
samtools view -bS map_PG_no_wheat.sam > map_PG_no_wheat.bam #On les transforme en fichier bam
samtools sort -o map_PG_no_wheat.sorted.bam map_PG_no_wheat.bam #On trie les fichiers bam
samtools rmdup map_PG_no_wheat.sorted.bam final.map_PG_no_wheat.bam #o enleve les duplicats de PCR
samtools index final.map_PG_no_wheat.bam #On les index

#Trinity 
sudo apt update
sudo apt install trinity
#Blast
# Téléchargement de BLAST (remplacez le lien par la version souhaitée)
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.12.0+-x64-linux.tar.gz

# Décompression de l'archive téléchargée
tar -xzvf ncbi-blast-2.12.0+-x64-linux.tar.gz
# Modifier le fichier ~/.bashrc pour ajouter BLAST à votre PATH
echo 'export PATH=$PATH:/chemin/vers/le/dossier/blast/bin' >> ~/.bashrc
source ~/.bashrc
blastn -version
