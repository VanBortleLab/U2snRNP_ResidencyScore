# Directory of where I put the data in cluster 
fastq_gz_directory = "/home/labs/kvbortle_lab/lab-fastq/unpublished/mRNA/snar_kd_ym155_mrna/fastq/"
sh_file_directory = "/home/labs/kvbortle_lab/sihangz2/RNAseq/snar_kd_ym155_mRNA/sh_files/"
trim_directory = "/home/labs/kvbortle_lab/sihangz2/RNAseq/snar_kd_ym155_mRNA/trim/"
rscript_directory = "/home/labs/kvbortle_lab/sihangz2/RNAseq/snar_kd_ym155_mRNA/R_script/"
star_directory = "/home/labs/kvbortle_lab/sihangz2/RNAseq/snar_kd_ym155_mRNA/star/"


#create an object of the metadata file 


m <- read.table("/home/labs/kvbortle_lab/sihangz2/RNAseq/snar_kd_ym155_mRNA/metadata.txt",header = T, sep = "\t")
#create a loop to write sh_fil for every Run
for(i in 1:nrow(m)){
  info=m$sample[i]
  sh_fil = file(paste(sh_file_directory, info ,"_star_pipeline.sh",sep="")) 
  writeLines(
    c(
      "#!/bin/bash",
      "# ----------------SLURM Parameters----------------",
      "#SBATCH -p normal",
      "#SBATCH -n 1",
      "#SBATCH --mem=200g",
      "#SBATCH -N 1",
      paste("#SBATCH -J", info),
      paste("#SBATCH -D", sh_file_directory),
      "#SBATCH -o pipeline_output.txt",
      "#SBATCH -e pipeline_error_paired.txt",
      "# ----------------Load Modules--------------------",
   
      "# ----------------Commands------------------------",
      "module purge",
      "module load Trim_Galore/0.6.5-IGB-gcc-4.9.4",
      
      paste("trim_galore --dont_gzip -o ",trim_directory, " --paired ", 
            fastq_gz_directory, m$trim_name[i], "_R1_001.fastq.gz ", fastq_gz_directory, m$trim_name[i],"_R2_001.fastq.gz ",
            " --fastqc",sep=""),
    "module purge",
    "module load STAR/2.7.6a-IGB-gcc-8.2.0",
      
      paste("STAR --genomeDir /home/labs/kvbortle_lab/lab-tools/STAR_Index/GenomeDir --readFilesIn ",
            trim_directory, m$trim_name[i],"_R1_001_val_1.fq ",trim_directory, m$trim_name[i], "_R2_001_val_2.fq",
            " --outTmpDir /scratch/$SLURM_JOB_ID",
            " --outFileNamePrefix ", star_directory, info,
            " --outSAMtype BAM SortedByCoordinate", sep = "")
      
      
      
      
      
    ),
    sh_fil)
  
  close(sh_fil)
  
  system(paste("sbatch", " ", sh_file_directory, info,"_star_pipeline.sh",sep=""))
  
}
