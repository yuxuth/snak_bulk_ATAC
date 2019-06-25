shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")
shell.suffix("; exitstat=$?; echo END at $(date); echo exit status was $exitstat; exit $exitstat")

configfile: "config.yaml"
 
FILES = json.load(open(config['SAMPLES_JSON']))

# CLUSTER = json.load(open(config['CLUSTER_JSON']))

SAMPLES = sorted(FILES.keys())

gtf = config["MYGTF"]
STARINDEX = config['STARINDEX']
hisat = config["hisat"]
stringtie = config['stringtie']
TARGETS = []
bbbuk = config['bbbuk']
adaptors = config['adaptors']
splicesite_index = config['splicesite_index']
## constructe the target if the inputs are fastqs
ALL_FASTQ   = expand("01_merged_seq/{sample}_{read}.fastq.gz", sample = SAMPLES, read = ["R1"]) ## set the SE and PE
ALL_TRIMMED_FASTQ = expand("02_trim_seq/{sample}_{read}.trimmed.fastq.gz", sample = SAMPLES, read = ["R1"])
ALL_FASTQC  = expand("03_fqc/{sample}_{read}.trimmed_fastqc.zip", sample = SAMPLES, read = ["R1"])
ALL_BAM = expand("04_bam/{sample}_Aligned.out.sam", sample = SAMPLES)
ALL_SORTED_BAM = expand("05_sortBam/{sample}.sorted.bam", sample = SAMPLES)
ALL_stringtie_gtf = expand("06_ballgown/{sample}/{sample}.stringtie.gtf", sample = SAMPLES)
ALL_bw = expand("07_bigwig/{sample}_forward.bw", sample = SAMPLES)
TARGETS.extend(ALL_BAM) ##append all list to 
TARGETS.extend(ALL_SORTED_BAM)
TARGETS.extend(ALL_stringtie_gtf)
TARGETS.extend(ALL_FASTQC) ## check later
TARGETS.extend(ALL_FASTQ)
TARGETS.extend(ALL_bw)


## construct the target if the inputs are bams


localrules: all
# localrules will let the rule run locally rather than submitting to cluster
# computing nodes, this is for very small jobs

rule all:
	input: TARGETS


rule merge_fastqs: ## merge fastq
	input:
		# r1 = lambda wildcards: FILES[wildcards.sample]['R1']
		r1 = lambda wildcards: FILES[wildcards.sample]['R1']
		# r2 = lambda wildcards: FILES[wildcards.sample]['R2']
	output:
		"01_merged_seq/{sample}_R1.fastq.gz" #, "01_merged_seq/{sample}_R2.fastq.gz"
	log: "00_log/{sample}_merge_fastq.log"
	params:
		jobname = "{sample}"
	group: "mygroup"
	threads: 1
	message: "merging fastqs {input}: {threads} threads"
	shell:
		"""
		gunzip -c {input.r1} | gzip > {output[0]} 2> {log}
		# echo {wildcards.sample}
		
		"""
		# gunzip -c {input.r2} | gzip > {output[1]} 2>> {log} ## the # need to locate outside of 

rule trim_adapter:
 	input: "01_merged_seq/{sample}_R1.fastq.gz" #, "01seq/{sample}_R2.fastq.gz"
 	output: "02_trim_seq/{sample}_R1.trimmed.fastq.gz" #, "03trim/{sample}_R2.trimmed.fastq.gz"
 	log: "00log/{sample}_trim_adaptor.log"
 	threads: 1
 	group: "mygroup"
 	params : 
 		jobname = "{sample}"
 	message: "trim_adaptor {input}: {threads}"
 	shell:
 		"""
		{bbbuk} in={input} out={output} ref={adaptors} ktrim=r k=23 mink=11 hdist=1 &> {log}
 		"""

rule fastqc:
	input:  "02_trim_seq/{sample}_R1.trimmed.fastq.gz" #, "01seq/{sample}_R2.fastq.gz"
	output: "03_fqc/{sample}_R1.trimmed_fastqc.zip" #, "02fqc/{sample}_R2_fastqc.zip"
	log:    "00log/{sample}_fastqc"
	threads: 1
	group: "mygroup"
	params : jobname = "{sample}"
	message: "fastqc {input}: {threads}"
	shell:
	    """
	    module load fastqc
	    fastqc -o 03_fqc -f fastq --noextract {input}  2> {log}
	    """


rule hisat_mapping:
	input: 
		"02_trim_seq/{sample}_R1.trimmed.fastq.gz"
	output: "04_bam/{sample}_Aligned.out.sam"
	log: "00log/{sample}_hisat_align.log"
	params: 
		jobname = "{sample}",
		# outprefix = "01bam_fq/{sample}"
	threads: 24
	group: "mygroup"
	message: "aligning {input} using hisat: {threads} threads"
	shell:
		"""
		{hisat} -p {threads} \
		--dta\
		--rna-strandness R \
		-x {STARINDEX} \
		-U {input} \
		--known-splicesite-infile {splicesite_index} \
		-S {output} \
		2> {log}
		"""

rule sortBam:
	input: "04_bam/{sample}_Aligned.out.sam"
	output: "05_sortBam/{sample}.sorted.bam"
	log: "00log/{sample}_sortbam.log"
	params:
		jobname = "{sample}"
	threads: 12
	group: "mygroup"
	message: "sorting {input} : {threads} threads"
	shell:
		"""
		module load samtools
		samtools sort  -@ {threads} -T {output}.tmp -o {output} {input} 2> {log}
		"""

rule index_bam:
    input:  "05_sortBam/{sample}.sorted.bam"
    output: "05_sortBam/{sample}.sorted.bam.bai"
    log:    "00log/{sample}.index_bam"
    threads: 1
    params: jobname = "{sample}"
    message: "index_bam {input}: {threads} threads"
    shell:
        """
        samtools index {input} 2> {log}
        """

rule stringtie_FPKM_caculation:
	input: "05_sortBam/{sample}.sorted.bam"
	output: "06_ballgown/{sample}/{sample}.stringtie.gtf"
	log: "00log/{sample}_stringtie.log"
	params:
		jobname = "{sample}"
	threads: 24
	group: "mygroup"
	message: "stringtie {input} : {threads} threads"
	shell:
		"""
		# -p for paried-end, counting fragments rather reads
		{stringtie} -e -B -p 8 -G {gtf} -o {output} {input}
		"""
## add the bowgrow later

rule make_bigwigs:
    input : "05_sortBam/{sample}.sorted.bam", "05_sortBam/{sample}.sorted.bam.bai"
    # input : expand("05_sortBam/{sample}.sorted.bam", sample = SAMPLES)
    output: "07_bigwig/{sample}_forward.bw", "07_bigwig/{sample}_reverse.bw"
    log: "00log/{sample}.makebw"
    threads: 64
    params: jobname = "{sample}"
    message: "making bigwig for {input} : {threads} threads"
    shell:
        """
    	# no window smoothing is done, for paired-end, bamCoverage will extend the length to the fragement length of the paired reads
        bamCoverage -b {input[0]}  --skipNonCoveredRegions --normalizeUsing RPKM --samFlagExclude 16  -p {threads}  -o {output[0]} 2> {log}
        bamCoverage -b {input[0]}  --skipNonCoveredRegions --normalizeUsing RPKM --samFlagInclude 16  -p {threads}  -o {output[1]} 2>> {log}
        """

rule update_bigwigs:
	input : ALL_bw
