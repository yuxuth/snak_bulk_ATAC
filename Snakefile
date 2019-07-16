shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")
shell.suffix("; exitstat=$?; echo END at $(date); echo exit status was $exitstat; exit $exitstat")

configfile: "config.yaml"
 
FILES = json.load(open(config['SAMPLES_JSON']))

# CLUSTER = json.load(open(config['CLUSTER_JSON']))

SAMPLES = sorted(FILES.keys())

bowtie2_index = config['bowtie2_index'] ## update as the bowtie index 
Genrich = config['Genrich']
bbbuk = config['bbbuk']
adaptors = config['adaptors']

## constructe the target if the inputs are fastqs
ALL_FASTQ   = expand("01_merged_seq/{sample}_{read}.fastq.gz", sample = SAMPLES, read = ["R1", "R2"]) ## set the SE and PE
ALL_TRIMMED_FASTQ = expand("02_trim_seq/{sample}_{read}.trimmed.fastq.gz", sample = SAMPLES, read = ["R1", "R2"])
ALL_FASTQC  = expand("03_fqc/{sample}_{read}.trimmed_fastqc.zip", sample = SAMPLES, read = ["R1", "R2"])
##
ALL_BAM = expand("04_bam/{sample}_sorted.out.bam", sample = SAMPLES)

ALL_SORTED_BAM = expand("05_sortBam/{sample}.sorted.bam", sample = SAMPLES)

ALL_BIGWIG = expand("07bigwig/{sample}.bw", sample = SAMPLES)

ALL_PEAKS = expand("06_peak/{sample}.peak", sample = SAMPLES)


# ALL_NUCLEO = expand("09nucleoATAC/{sample}_nucleoATAC.occpeaks.bed.gz", sample = ALL_SAMPLES)
# ALL_QC = ["10multiQC/multiQC_log.html"]
# ALL_ATAQV = expand("04aln/{sample}.sorted.bam.ataqv.json", sample = ALL_SAMPLES)
# ALL_ATAC_QC = ["11ATAC_qc_html"]

TARGETS = []
TARGETS.extend(ALL_FASTQ)
TARGETS.extend(ALL_FASTQC) 
TARGETS.extend(ALL_BAM) 
# # TARGETS.extend(ALL_SORTED_BAM)
# # TARGETS.extend(ALL_BIGWIG)
TARGETS.extend(ALL_PEAKS)
# TARGETS.extend(ALL_NUCLEO)
# TARGETS.extend(ALL_QC)
# TARGETS.extend(ALL_ATAQV)
# TARGETS.extend(ALL_ATAC_QC)

localrules: all
# localrules will let the rule run locally rather than submitting to cluster
# computing nodes, this is for very small jobs

rule all:
	input: TARGETS


rule merge_fastqs: ## merge fastq
	input:
		r1 = lambda wildcards: FILES[wildcards.sample]['R1'],
		r2 = lambda wildcards: FILES[wildcards.sample]['R2']
	output:
		"01_merged_seq/{sample}_R1.fastq.gz" , "01_merged_seq/{sample}_R2.fastq.gz"
	log: "00_log/{sample}_merge_fastq.log"
	params:
		jobname = "{sample}"
	# group: "mygroup"
	threads: 1
	message: "merging fastqs {input}: {threads} threads"
	run:
		if len(input.r1) == 1:
			shell("ln -sr {input.r1} {output[0]} && ln -sr {input.r2} {output[1]}")
			print("generate link")
		else:
			shell("bgunzip -c {input.r1} | gzip > {output[0]} 2> {log} && gunzip -c {input.r2} | gzip > {output[1]} 2>> {log}")
			print("merge_fastq")

rule trim_adapter:
 	input: "01_merged_seq/{sample}_R1.fastq.gz" , "01_merged_seq/{sample}_R2.fastq.gz"
 	output: "02_trim_seq/{sample}_R1.trimmed.fastq.gz" , "02_trim_seq/{sample}_R2.trimmed.fastq.gz"
 	log: "00_log/{sample}_trim_adaptor.log"
 	threads: 1
 	# group: "mygroup"
 	params : 
 		jobname = "{sample}"
 	message: "trim_adaptor {input}: {threads}"
 	shell:
 		"""
		{bbbuk} in1={input[0]} in2={input[1]} out1={output[0]} out2={output[1]} ref={adaptors} ktrim=r k=21 mink=11 hdist=1 &> {log}
 		"""

rule fastqc:
	input:  "02_trim_seq/{sample}_R1.trimmed.fastq.gz" , "02_trim_seq/{sample}_R2.trimmed.fastq.gz"
	output: "03_fqc/{sample}_R1.trimmed_fastqc.zip" , "03_fqc/{sample}_R2.trimmed_fastqc.zip"
	log:    "00_log/{sample}_fastqc"
	threads: 1
	# group: "mygroup"
	params : jobname = "{sample}"
	message: "fastqc {input}: {threads}"
	shell:
	    """
	    module load fastqc
	    fastqc -o 03_fqc -f fastq --noextract {input}  2> {log}
	    """

rule bowtie_mapping:
	input: 
		"02_trim_seq/{sample}_R1.trimmed.fastq.gz", "02_trim_seq/{sample}_R2.trimmed.fastq.gz"
	output: "04_bam/{sample}_sorted.out.bam"
	log: "00log/{sample}_bowtie_align.log"
	params: 
		jobname = "{sample}",
		# outprefix = "01bam_fq/{sample}"
	threads: 12
	# group: "mygroup"
	message: "aligning {input} using bowtie2: {threads} threads"
	shell:
		"""
		module load bowtie samtools
		bowtie2 --threads {threads} --very-sensitive  -k 10 -x {bowtie2_index} \
		 -1 {input[0]} -2 {input[1]} 2> {log} \
		 | samtools view  -u  - \
		 | samtools sort -m 2G -@ 5 -T {input[1]}.tmp -n  -o {output}  -
		"""
## change to wi

# rule index_bam:
#     input:  "05_sortBam/{sample}.sorted.bam"
#     output: "05_sortBam/{sample}.sorted.bam.bai"
#     log:    "00log/{sample}.index_bam"
#     threads: 1
#     params: jobname = "{sample}"
#     message: "index_bam {input}: {threads} threads"
#     shell:
#       """
#         samtools index {input} 2> {log}
# 		"""

rule Genrich_calling:
    input:  "04_bam/{sample}_sorted.out.bam"
    output: "06_peak/{sample}.peak"
    log:    "00log/{sample}.peak_calling"
    threads: 1
    params: jobname = "{sample}"
    message: "index_bam {input}: {threads} threads"
    shell:
        """
        {Genrich}  -t {input}  -o {output}  -j  -y  -r  -e chrM  -v 2> {log}
        """

# ## add the bowgrow later

# rule make_bigwigs:
#     input : "05_sortBam/{sample}.sorted.bam", "05_sortBam/{sample}.sorted.bam.bai"
#     # input : expand("05_sortBam/{sample}.sorted.bam", sample = SAMPLES)
#     output: "07_bigwig/{sample}_forward.bw", "07_bigwig/{sample}_reverse.bw"
#     log: "00log/{sample}.makebw"
#     threads: 64
#     params: jobname = "{sample}"
#     message: "making bigwig for {input} : {threads} threads"
#     shell:
#         """
#     	# no window smoothing is done, for paired-end, bamCoverage will extend the length to the fragement length of the paired reads
#         bamCoverage -b {input[0]}  --skipNonCoveredRegions --normalizeUsing RPKM --samFlagExclude 16  -p {threads}  -o {output[0]} 2> {log}
#         bamCoverage -b {input[0]}  --skipNonCoveredRegions --normalizeUsing RPKM --samFlagInclude 16  -p {threads}  -o {output[1]} 2>> {log}
#         """

# rule update_bigwigs:
# 	input : ALL_bw
