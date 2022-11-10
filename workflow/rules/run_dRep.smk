from glob import glob
from pathlib import Path
from shutil import copyfile, rmtree, copyfileobj
import gzip

rule drep:
    input:
        # From https://bioinformatics.stackexchange.com/questions/7184/mix-globbing-and-wildcards-when-specifying-rule-input
        glob(join(config["bindir"],"*.fa"))
    output:
        "workflow/out/dRep_output/data_tables/genomeInfo.csv"
    params:
        outdir = "workflow/out/dRep_output/"
    threads:
        config['threads']
    conda:
        "../../workflow/envs/dreppython37_yay.yml"
    shell:
        "dRep compare {params.outdir} -g {input} -p {threads}"
   