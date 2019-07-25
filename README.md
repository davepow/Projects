# RNAseq-NF ENCODE pipeline 

NOTE: this is a clone of this repository: 

https://github.com/nextflow-io/rnaseq-encode-nf


with a very slight modification to main.nf:

Switch from ftp URLs to http URLs:

params.transcriptome = "http://ftp.ensembl.org/pub/release-90/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"


url = "http://ftp.sra.ebi.ac.uk/vol1/fastq/{0}/".format(dbxref[0:6])


