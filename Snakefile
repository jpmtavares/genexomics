configfile: "snakeconfig.yml"

READS = ["R1", "R2"]
LANES = ["L001", "L002"]

SAMPLES = config["SAMPLES"]
SOURCE_FASTQ_DIR = config["SOURCE_FASTQ_DIR"]
TRANSCRIPTOME = config["TRANSCRIPTOME_DIR"]

rule all:
    input:
        expand("{sample_name}_{sample_number}/fastq/{sample_name}_{sample_number}_{lane}_{read}_001.fastq.gz",
               sample_name=SAMPLES.values(), sample_number=SAMPLES.keys(), lane=LANES, read=READS),
        expand("{sample_name}_{sample_number}/cellranger/count",
               sample_name=SAMPLES.values(), sample_number=SAMPLES.keys()),
        expand("{sample_name}_{sample_number}/analysis/plots",
               sample_name=SAMPLES.values(), sample_number=SAMPLES.keys())

rule copy_raw_data:
    input:
        original_fastq=lambda wildcards: config["SOURCE_FASTQ_DIR"] + f"{wildcards.sample_name}_{wildcards.sample_number}_{wildcards.lane}_{wildcards.read}_001.fastq.gz"
    output:
        "{sample_name}_{sample_number}/fastq/{sample_name}_{sample_number}_{lane}_{read}_001.fastq.gz"
    shell:
        """
        cp {input.original_fastq} {output}
        """

rule cellranger_count:
    input:
        fastq_files=expand("{sample_name}_{sample_number}/fastq/{sample_name}_{sample_number}_{lane}_{read}_001.fastq.gz",
                           sample_name=SAMPLES.values(), sample_number=SAMPLES.keys(), lane=LANES, read=READS),
        transcriptome=TRANSCRIPTOME
    output:
        directory("{sample_name}_{sample_number}/cellranger/count")
    params:
        id="{sample_name}_{sample_number}",
        sample="{sample_name}",
        fastq_dir="{sample_name}_{sample_number}/fastq/"
    shell:
        """
        cellranger count --id={params.id} --transcriptome={input.transcriptome} \
        --fastqs={params.fastq_dir} --sample={params.sample} --create-bam=true \
        --output-dir={output} \
        --localcores 4
        """

rule plot_matrix:
    input:
        matrix_dir=expand("{sample_name}_{sample_number}/cellranger/count",
                          sample_name=SAMPLES.values(), sample_number=SAMPLES.keys())
    output:
        directory("{sample_name}_{sample_number}/analysis/plots")
    shell:
        """
        python matrix_plotting.py --input_dir {input.matrix_dir}/outs/filtered_feature_bc_matrix/
        """
