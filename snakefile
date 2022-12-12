try:
	shell("mkdir txt_files")
except:
	#print("le dossier txt_files existe déjà :) ")
	pass
try: 
	shell("mkdir sra_files")
except:
	#print("le dossier sra_files existe déjà")
	pass
try:
	shell(" mkdir fastq_files")
except: 
	#print("le dossier fastq_files existe déjà :) ")
	pass

try:
	shell(" mkdir ref")
except: 
	#print("le dossier ref existe déjà :) ")
	pass


shell("chmod +x txt_files") # on s'assure qu'il n'y aura pas de problème de droit d'accès
shell("chmod +x fastq_files") 
shell("chmod +x sra_files")
shell("chmod +x ref") 


LIST_SRA = ['SRR628589', 'SRR628588', 'SRR628587', 'SRR628586', 'SRR628585', 'SRR628584', 'SRR628583', 'SRR628582']

# on crée un fichier txt pour chaque numero d'accession
for item in LIST_SRA :
	with open("txt_files/{item}.txt".format(item = item), "w") as f:
		f.write(str(item))
	
rule all:
	input : "condition_cancer_vs_control_dge.csv"

rule import_rsa:
	input : "txt_files/{list_sra}.txt"
	output: "sra_files/{list_sra}.sra"
	container : "docker://pegi3s/sratoolkit:2.10.0"
	shell: 	
		"""
		prefetch -v {wildcards.list_sra} > sra_files/{wildcards.list_sra}.sra
		rm -rf {wildcards.list_sra}
		"""
	
rule fastq_dump:
	input:"sra_files/{list_sra}.sra"
	output: "fastq_files/{list_sra}_1.fastq", "fastq_files/{list_sra}_2.fastq"
	container : "docker://pegi3s/sratoolkit:2.10.0"
	shell: 
		"""
		fastq-dump --split-files {wildcards.list_sra} --outdir fastq_files/
		"""

rule fastqc:
	input: expand("fastq_files/{list_sra}_1.fastq", list_sra = LIST_SRA), expand("fastq_files/{list_sra}_2.fastq", list_sra = LIST_SRA)
	output: directory('fastqc_files') 
	container : "docker://biocontainers/fastqc:v0.11.9_cv8"
	shell:
    		"""
		mkdir fastqc_files
   		fastqc "fastq_files/{list_sra}_1.fastq", "fastq_files/{list_sra}_2.fastq" -O fastqc_files/
   		"""

# Maintenant, passons à STAR. 
# Dans un premier temps, on télécharge le génome de référence. 
# On télécharge tous les chromosomes et on les concatene en un unique fichier de référence.

print("Creation de l'index")
rule get_genome_ref:
	output: fa = "ref.fa",
		gtf = "ref.gtf"
	shell: 
		"""
		source get_genome_ref
		""" 

rule index:
	input: 
		fa = "ref.fa",
		gtf = "ref.gtf",
	output: directory('index_done')
	threads: 16
	container: "docker://antoinedegenne/star:latest"
	shell: 
		"""
		STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir ref --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf} --sjdbOverhang 100
		mkdir index_done   # on crée ce dossier pour que cette rule se lance effectiveùent durant le snkfl 
		"""
 
print('avant mapping')

rule mapping:
	input:
		index_done = directory('index_done'),
		fastq1 = "fastq_files/{list_sra}_1.fastq",
		fastq2 = "fastq_files/{list_sra}_2.fastq"
		

	output: 
		bam1 = "{list_sra}Aligned.sortedByCoord.out.bam"
	
	threads : 16 
	container : "docker://antoinedegenne/star:latest"
	shell: 
		'STAR --runThreadN {threads} '
		'--genomeDir ref ' 
		'--runMode alignReads '		
		'--readFilesIn {input.fastq1} {input.fastq2} '
		'--outFileNamePrefix {wildcards.list_sra} '
		'--outSAMtype BAM SortedByCoordinate '
		
rule feature_counts:
	input :
		gtf = "ref.gtf",
		bam = "{list_sra}Aligned.sortedByCoord.out.bam"
	
	output : 
		counts = "{list_sra}Counts"
	threads : 16
	container : "docker://drylse/featurecounts:featurecounts"
	shell : 
		'featureCounts -T {threads} -t gene -g gene_id -p -s 0 -a {input.gtf} -o {output.counts} {input.bam}'

rule deseq2:
	input: expand("{list_sra}Counts", list_sra = LIST_SRA)
	output: "condition_cancer_vs_control_dge.csv"
	container : "docker://nanozoo/deseq2:1.34.0--c670fa0"
	script: "deseq2.R"

						
