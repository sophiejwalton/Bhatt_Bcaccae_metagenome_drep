from os.path import join
import pandas as pd
from glob import glob

configfile: "config/config.yaml"

rule all:
    input:
        "workflow/out/dRep_output/data_tables/genomeInfo.csv",
        
# Check the dag
#  snakemake --forceall --dag | dot -Tpng >dag.png
include: "workflow/rules/run_dRep.smk"



