from pathlib import Path

INDIR = Path("inputs")
OUTDIR = Path("results")
SAMPLESHEET = INDIR / "samples.csv"
include: "rules/common.snk"

dictReads = parseSamples(SAMPLESHEET)
BARCODES = list(dictReads.keys())
SUBREADS = [subdict["limit"] for subdict in dictReads.values()]

rule all:
    input:
        expand(
            OUTDIR / "subsampled" / "{barcode}.{limit}.{repeat_counter}.txt", 
            barcode = BARCODES,
            limit = SUBREADS,
            repeat_counter = range(0, 5)
            )

rule rasusa:
    input:
        R1 = get_input_r1, 
        R2 = get_input_r2,
    output:
        status = OUTDIR / "subsampled" / "{barcode}.{limit}.replicate{repeat_counter}.txt",
        R1 = OUTDIR / "subsampled" / "{barcode}" / "{barcode}_{limit}_replicate{repeat_counter}_R1.fastq.gz",
        R2 = OUTDIR / "subsampled" / "{barcode}" / "{barcode}_{limit}_replicate{repeat_counter}_R2.fastq.gz",
    params:
        limit = "{limit}",

    conda: "config/rasusa.yaml"
    threads: 1
    shell:"""
    rasusa reads \
        {input.R1} \
        {input.R2} \
        --num {params.limit} \
        -o {output.R1} \
        -o {output.R2}
    touch {output.status}
    """

# rule vircov:
#     input:
#         status = rules.rasusa.output.status,
#         R1 = rules.rasusa.output.params.R1,
#         R2 = rules.rasusa.output.params.R2,
#     output:
#         results = OUTDIR / "vircov" / "{barcode}" / "{barcode}_{limit}_{repeat_counter}.vircov.tsv",
#     params:
#         wkdir = temp(OUTDIR / "vircov" / "{barcode}" / "wkdir" / "{repeat_counter}_{limit}"),
#         vircovArgs = "--consensus-min-depth 3 --min-depth-coverage 1 --annotation-preset default --remap-exclude-bins 'Homo sapiens'",
#         index = "/raid/VIDRL-USERS/RESOURCES/DATABASES/ICTV/",
#         reference = "/raid/VIDRL-USERS/RESOURCES/DATABASES/ICTV/ictv.fasta",
#     conda: "config/vircov.yaml"
#     shell:"""
#     vircov run \
#         -i {input.R1} \
#         -i {input.R2} \
#         -o {output.results} \
#         --aligner bowtie2 \
#         --index {params.index} \
#         --reference {params.reference} \
#         --scan-threads {threads} \
#         --remap-threads {threads} \
#         --remap-parallel {threads} \
#         --workdir {params.wkdir} \
#         --keep {params.vircovArgs}
#     """