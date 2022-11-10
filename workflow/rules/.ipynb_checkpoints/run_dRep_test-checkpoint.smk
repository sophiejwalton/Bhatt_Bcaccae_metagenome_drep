from glob import glob
from pathlib import Path
from shutil import copyfile, rmtree, copyfileobj
import gzip

rule drep:
    input:
        # From https://bioinformatics.stackexchange.com/questions/7184/mix-globbing-and-wildcards-when-specifying-rule-input
        lambda wildcards: glob(join(config["bindir"],"{cow}_*.fa".format(cow=wildcards.cow)))
    output:
        "workflow/out/dRep_output/{cow}/data_tables/genomeInfo.csv"
    params:
        outdir = "workflow/out/dRep_output/{cow}/",
        genome_info = "config/genomeInformation.csv"
    conda:
        "../../workflow/envs/dreppython37_yay.yml"
    shell:
        """
        dRep dereplicate {params.outdir} -g {input} -pa 0.5 -sa 0.95 -p 7 
        """
        # dereplicate {params.outdir} -g {input} -pa 0.5 -sa 0.95 -p 7"
        #config['dRepdir'] + " dereplicate {output} -g {input} --completeness 0.8 -pa 0.95 \
        #--run_tax --S_algorithm fastANI -p 7 --cent_index /Users/lnmerk/centrifuge-1.0.3-beta/indices/p+h+v"

