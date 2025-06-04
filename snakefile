from pathlib import Path

INDIR = Path("inputs")
OUTDIR = Path("results")
SAMPLESHEET = INDIR / "samples.csv"
include: "rules/common.snk"

dictReads = parseSamples(SAMPLESHEET)
BARCODES = list(dictReads.keys())
NUM_REPEAT = 2

rule all:
    input:
        expand(
            OUTDIR / "subsampled" / "{barcode}.replicate{repeat_counter}.txt", 
            barcode = BARCODES,
            repeat_counter = range(0, NUM_REPEAT)
        ),
        OUTDIR / "vircov" / "all.combined.tsv"

rule rasusa:
    input:
        R1 = get_input_r1, 
        R2 = get_input_r2,
    output:
        status = OUTDIR / "subsampled" / "{barcode}.replicate{repeat_counter}.txt",
        R1 = OUTDIR / "subsampled" / "{barcode}" / "{barcode}_replicate{repeat_counter}_R1.fastq.gz",
        R2 = OUTDIR / "subsampled" / "{barcode}" / "{barcode}_replicate{repeat_counter}_R2.fastq.gz",
        log = OUTDIR / "logs" / "{barcode}.replicate{repeat_counter}.txt",
    params:
        limit = get_input_limit,
    conda: "config/rasusa.yaml"
    threads: 10
    shell:"""
    rasusa reads \
        {input.R1} \
        {input.R2} \
        --num {params.limit} \
        -o {output.R1} \
        -o {output.R2} 2> {output.log}
    touch {output.status}
    """

rule vircov:
    input:
        status = rules.rasusa.output.status,
        R1 = rules.rasusa.output.R1,
        R2 = rules.rasusa.output.R2,
    output:
        results = OUTDIR / "vircov" / "{barcode}" / "{barcode}-replicate{repeat_counter}.tsv",
        log = OUTDIR / "logs" / "{barcode}" / "{barcode}.replicate{repeat_counter}.vircov.txt",
    params:
        wkdir = temp(OUTDIR / "vircov" / "{barcode}" / "wkdir" / "{repeat_counter}"),
        vircovArgs = "--consensus-min-depth 3 --min-depth-coverage 1 --annotation-preset default --remap-exclude-bins 'Homo sapiens'",
        index = "/raid/VIDRL-USERS/RESOURCES/DATABASES/ICTV/ictv",
        reference = "/raid/VIDRL-USERS/RESOURCES/DATABASES/ICTV/ictv.fasta",
    conda: "config/vircov.yaml"
    threads: 10
    shell:"""
    vircov run \
        -i {input.R1} \
        -i {input.R2} \
        -o {output.results} \
        --aligner bowtie2 \
        --index {params.index} \
        --reference {params.reference} \
        --scan-threads {threads} \
        --remap-threads {threads} \
        --remap-parallel {threads} \
        --workdir {params.wkdir} \
        --keep {params.vircovArgs} 2> {output.log}
    """

rule concat_vircov:
    input:
        expand(OUTDIR / "vircov" / "{barcode}" / "{barcode}-replicate{repeat_counter}.tsv",
               barcode = BARCODES,
               repeat_counter = range(0, NUM_REPEAT)
        )
    output:
        OUTDIR / "vircov" / "all.combined.tsv"
    params:
        min_completeness = 30
    conda: "config/vircov.yaml"
    shell:
        """
        vircov tools concat-output \
            --input {input} \
            --output {output} \
            --file-id
        """
