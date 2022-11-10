from glob import glob

def get_species_list(cow):
    all_file_1_timepoint = glob(f"/Users/lnmerk/midas_db_v1.2/rep_genomes/{cow}_*")
    species_list = []

    for file in all_file_1_timepoint:
        species_list.append(file.split("/")[-1])

    return ",".join(species_list)

rule MIDAS_snps:
    input:
        "workflow/out/complete.txt",
        r1 = lambda wildcards: glob(join(config["fastqdir"],"{cow}_{timepoint}_*R1*.fastq.gz".format(cow=wildcards.cow, timepoint=wildcards.timepoint))),
        r2 = lambda wildcards: glob(join(config["fastqdir"],"{cow}_{timepoint}_*R2*.fastq.gz".format(cow=wildcards.cow, timepoint=wildcards.timepoint)))

    output:
        "workflow/out/midasOutput/{cow}_{timepoint}/snps/summary.txt"
    params:
        lambda wildcards: get_species_list(wildcards.cow)
    threads: config['maxCPUs']
    shell:
        "/Users/lnmerk/MIDAS/scripts/run_midas.py snps workflow/out/midasOutput/{wildcards.cow}_{wildcards.timepoint} \
         -1 {input.r1} -2 {input.r2} -t {threads} --species_id {params}"

def get_samples_list(cow):
    all_file_1_cow = glob(f"workflow/out/midasOutput/{cow}_*")

    return ",".join(all_file_1_cow)

rule mergeSNPs:
    input:
        lambda wildcards: expand("workflow/out/midasOutput/{cow}_{timepoint}/snps/summary.txt",cow = wildcards.cow, timepoint = cow_timepoints[wildcards.cow])
    output:
        "workflow/out/{cow}/complete_merge.txt"
    params:
        sampledirs=lambda wildcards: get_samples_list(wildcards.cow),
        outdir="workflow/out/midasOutput/snps/"
    threads: config['maxCPUs']
    shell:
        """
        mkdir {params.outdir}
        /Users/lnmerk/MIDAS/scripts/merge_midas.py snps {params.outdir}/{wildcards.cow} -t list -i {params.sampledirs}  \
        --sample_depth 1 --site_depth 1 --min_samples 1 --max_species 150 --site_prev 0.0 --threads {threads} \
        --fract_cov 0 --all_sites --all_samples
        echo > workflow/out/{wildcards.cow}/complete_merge.txt
        """
