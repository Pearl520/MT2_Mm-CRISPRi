config {

  executor="lsf"
  
  commands {   
    Fastqc {
      lsf_request_options="-W 01:00 -M 2000 -n 6"
    }
	
    Trimgalore {
      lsf_request_options="-W 04:00 -M 1000 -n 2"  	
    }	

    Bowtie2 {
	   lsf_request_options="-W 06:00 -M 15000 -n 6"
    }

    Sortbam {
      lsf_request_options="-W 01:00 -M 10000 -n 2"
    }			

    FilterbyMappingQual {
      lsf_request_options="-W 01:00 -M 10000 -n 3"
    }

    Dedup {
      lsf_request_options="-W 01:00 -M 30000 -n 1"
    }

    Macs2 {
      lsf_request_options="-W 01:00 -M 5000 -n 3"
    }

    BamCoverage {
      lsf_request_options="-W 01:00 -M 10000 -n 5"
    }
  }
}

DataDir = "/data/ZYChenlab/Zhiyuan/projects/TE_project/230810_MT2_Mmi_ESCs_Cas9ChIPCUTRUN/raw_data"

def branches = [

	CTRi_ESCs_Cas9_ChIP_rep1        : ["${DataDir}/CTR_Cas9ChIP_1_CKDL230024360-1A_HCNJNDSX7_L3_1.fq.gz",
                                     "${DataDir}/CTR_Cas9ChIP_1_CKDL230024360-1A_HCNJNDSX7_L3_2.fq.gz"],

  CTRi_ESCs_HA_ChIP_rep1          : ["${DataDir}/CTR_HAChIP_1_CKDL230024360-1A_HCNJNDSX7_L3_1.fq.gz",
                                     "${DataDir}/CTR_HAChIP_1_CKDL230024360-1A_HCNJNDSX7_L3_2.fq.gz"],

  MT2Mmi_ESCs_Cas9_ChIP_rep2      : ["${DataDir}/MT2Mmi_Cas9ChIP_2_CKDL230024360-1A_HCNJNDSX7_L3_1.fq.gz",
                                     "${DataDir}/MT2Mmi_Cas9ChIP_2_CKDL230024360-1A_HCNJNDSX7_L3_2.fq.gz"],

  MT2Mmi_ESCs_HA_ChIP_rep2        : ["${DataDir}/MT2Mmi_HAChIP_2_CKDL230024360-1A_HCNJNDSX7_L3_1.fq.gz",
                                     "${DataDir}/MT2Mmi_HAChIP_2_CKDL230024360-1A_HCNJNDSX7_L3_2.fq.gz"],
]

Fastqc = {
	doc title: "QC using fastQC",
	desc: "Using fastQC to generate quality control reports",
	author: "zchen0423@gmail.com"

	var SOUT : "./FastQC"
	var threads: 6

	def basename1 = new File(input1).getName().prefix.replaceAll(".fq","");
	def basename2 = new File(input2).getName().prefix.replaceAll(".fq","");
 	produce("${SOUT}/${branch.name}/${basename1}_fastqc.html","${SOUT}/${branch.name}/${basename2}_fastqc.html"){
       	exec """
        	module load fastqc/0.11.7 &&
         	mkdir -p ${SOUT} && 
         	mkdir -p ${SOUT}/${branch.name} &&
         	fastqc --threads ${threads} --extract  --outdir "${SOUT}/${branch.name}" $inputs
	""","fastqc"
    }
	forward inputs
}

Trimgalore = {
	def basename1 = new File(input1).getName().prefix.replaceAll(".fq","");
	def basename2 = new File(input2).getName().prefix.replaceAll(".fq","");

	produce("./Trimgalore/${branch.name}_trimgalore/${basename1}_val_1.fq.gz","./Trimgalore/${branch.name}_trimgalore/${basename2}_val_2.fq.gz"){
	exec """
      		module load trimgalore/0.6.6 &&
      		mkdir -p Trimgalore && 
	    	mkdir -p Trimgalore/${branch.name}_trimgalore &&
      		trim_galore --paired  --fastqc --output_dir Trimgalore/${branch.name}_trimgalore $inputs
	""","Trimgalore"    
   }
}

bw2_genome = "/data/ZYChenlab/Zhiyuan/genomes_annotations/mm10/genomes/mm10"

Bowtie2 = {
	var genome : bw2_genome
   	var organism : "mm10"
   	var additional_flags : ""

   	def minisize=0
   	def maxisize=1000
   	def fqDir = "./Trimgalore/${branch.name}_trimgalore"

   	def r1 = new File(fqDir).list().find{it=~/_val_1.fq.gz/}
   	def r2 = new File(fqDir).list().find{it=~/_val_2.fq.gz/}

   	if(organism != "mm10"){
       		branch.spikein = "./bowtie2_mapping_${organism}/${branch.name}_${organism}.bam"
       		fqDir = "./unmapped_mm10"

       		r1 = "${branch.name}_1.fq" 
       		r2 = "${branch.name}_2.fq" 

       		maxisize = 700
       		minisize = 0
   	}

   	produce("./bowtie2_mapping_${organism}/${branch.name}_${organism}.bam"){
       		exec """
         		module load gcc/9.3.0 boost/1.79.0 samtools/1.9.0 bowtie2/2.4.2 &&
         		mkdir -p bowtie2_mapping_${organism} &&
         		(bowtie2 -x ${genome}
         		-1 ${fqDir}/$r1
         		-2 ${fqDir}/$r2
         		-p 6 --no-unal --no-mixed --no-discordant
         		-I $minisize
         		-X $maxisize
         		$additional_flags | samtools view -bS - > ./bowtie2_mapping_${organism}/${branch.name}_${organism}.bam) >& "./bowtie2_mapping_${organism}/${branch.name}_${organism}_bowtie2_log.txt"
       		""","Bowtie2"
     }
}

FilterbyMappingQual = {
	def basename = new File(input).getName().prefix
    	def qual = 30
	produce("./filter_qual_Q${qual}/${basename}.Q${qual}.bam"){
		exec """
			module load samtools/1.9.0 sambamba/0.6.8 &&
           		mkdir -p ./filter_qual_Q${qual} &&
			sambamba view -p -t 3 -f bam -F "mapping_quality >= ${qual} and not (unmapped or mate_is_unmapped)" $input > ./filter_qual_Q${qual}/${basename}.Q${qual}.bam
		""","FilterbyMappingQual"
	}
}

Sortbam = {    
	produce("${input.prefix}.sorted.multi.bam"){
	exec """
		module load samtools/1.9.0 &&
		samtools sort $input -o "${input.prefix}.sorted.multi.bam"
	""","Sortbam"
	}    
}

PICARD="/usr/local/picard/2.18.22/picard.jar"
Dedup = {
	def basename = new File(input).getName().prefix
 	branch.basename = basename
	produce("./deduplicated/${basename}.dedup.bam"){
		exec """
			module load java/1.7.0u40 picard/2.18.22 &&
			mkdir -p deduplicated &&
			java -jar $PICARD MarkDuplicates REMOVE_DUPLICATES=true 
			input=$input 
			output=./deduplicated/${basename}.dedup.bam
			METRICS_FILE=./deduplicated/${branch.name}_metrics.txt
			VALIDATION_STRINGENCY=SILENT
			CREATE_INDEX=true
		""","Dedup"
	}
}

Macs2 = {
	produce("./macs2/${branch.name}_multi_peaks.narrowPeak"){
	exec """
		module load MACS/2.1.4 &&
		mkdir -p macs2 &&
      		macs2 callpeak -t $input -f BAMPE -g mm --outdir ./macs2 -n ${branch.name}_multi -B -q 0.0001 --nolambda --nomodel
	""","Macs2"
	forward input 
  }
}

BamCoverage = {
	doc title : "Generate bigwig file"
        author: "zchen0423@gmail.com"
   	def basename = new File(input).getName().prefix
   	def scale_unit=10000
   	def effGenomeSize=2467481108 //https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html (2nd table)
   	produce("./bigwigs/${basename}.bw"){
     		exec """
       			module load deeptools/2.0.0 &&
       			mkdir -p ./bigwigs &&              
       			bamCoverage --bam $input.bam  
                   	--binSize 25 
                   	--extendReads 150
                   	--normalizeUsing RPKM 
                   	--outFileFormat bigwig                   
                   	--scaleFactor 1
                   	--numberOfProcessors 5
                   	-o ./bigwigs/${basename}.bw
    		 ""","BamCoverage"
	}
  	forward input
}

run {
  branches * [    Fastqc + 
	        	  Trimgalore +
                  Bowtie2.using(genome: bw2_genome, organism: "mm10") +                
                  Sortbam +             
                  Dedup + Macs2 + BamCoverage 
             ]
}
