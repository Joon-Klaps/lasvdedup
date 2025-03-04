configfile: "config.yaml"

# Accept command line parameters with defaults from config
CONTIGS_TABLE = config.get("contigs_table", "")
SEQ_DATA_DIR = config.get("seq_data_dir", "")

# Override with command line parameters if provided
if workflow.overwrite_config:
    if "contigs_table" in workflow.overwrite_config:
        CONTIGS_TABLE = workflow.overwrite_config["contigs_table"]
    if "seq_data_dir" in workflow.overwrite_config:
        SEQ_DATA_DIR = workflow.overwrite_config["seq_data_dir"]

SEGMENTS = config["segments"]

# Validate required parameters
if not CONTIGS_TABLE or not SEQ_DATA_DIR:
    raise ValueError("Missing required parameters: contigs_table and/or seq_data_dir")

rule all:
    input:
        expand("iqtree-out/lasv-{segment}.denovo.trimal.renamed.treefile", segment=SEGMENTS)

rule prepare_directories:
    output:
        directory("raw"),
        directory("iqtree-out")
    shell:
        """
        mkdir -p {output}
        """

rule prepare_base_alignments:
    output:
        "raw/lasv-base-{segment}.aln"
    input:
        "raw"
    shell:
        """
        touch {output}
        # Note: Copy your existing base alignments here or modify this rule
        # cp /path/to/your/base/alignment/{wildcards.segment}.aln {output}
        """

rule extract_sequences:
    output:
        "raw/lasv-{segment}.denovo.fasta"
    input:
        contigs=CONTIGS_TABLE,
        dir="raw"
    shell:
        """
        csvtk -t filter2 -f '${{(annotation) acronym}}=="LASV" && ${{(annotation) segment}}=="{wildcards.segment}"' {input.contigs} | \
        cut -f1 | grep -v "_consensus" | \
        xargs -I @ find {SEQ_DATA_DIR} -name "*@*" -type f -exec cat {{}} + > {output}
        """

rule align_sequences:
    output:
        "lasv-{segment}.denovo.aln"
    input:
        sequences="raw/lasv-{segment}.denovo.fasta",
        base_aln="raw/lasv-base-{segment}.aln"
    params:
        threads=config["mafft"]["threads"]
    shell:
        """
        mafft --add {input.sequences} --adjustdirection --thread {params.threads} {input.base_aln} > {output}
        """

rule trim_alignment:
    output:
        "lasv-{segment}.denovo.trimal.aln"
    input:
        "lasv-{segment}.denovo.aln"
    shell:
        """
        trimal -in {input} -out {output} -automated1
        """

rule build_tree:
    output:
        tree="iqtree-out/lasv-{segment}.denovo.trimal.treefile",
        log="iqtree-out/lasv-{segment}.denovo.trimal.log"
    input:
        alignment="lasv-{segment}.denovo.trimal.aln",
        outdir="iqtree-out"
    params:
        model=config["iqtree"]["model"],
        bootstrap=config["iqtree"]["bootstraps"],
        prefix="iqtree-out/lasv-{segment}.denovo.trimal"
    threads: workflow.cores
    shell:
        """
        iqtree -m {params.model} -B {params.bootstrap} -T {threads} -redo -s {input.alignment} \
        -pre {params.prefix} > {output.log} 2>&1
        """

rule rename_tree:
    output:
        "iqtree-out/lasv-{segment}.denovo.trimal.renamed.treefile"
    input:
        "iqtree-out/lasv-{segment}.denovo.trimal.treefile"
    shell:
        """
        sed -e 's/\.consensus_ivar//g' -e 's/_R_//g' -e 's/|/\//g' {input} > {output}
        """