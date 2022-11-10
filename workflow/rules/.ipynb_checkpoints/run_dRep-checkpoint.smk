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
        outdir = "workflow/out/dRep_output/{cow}/"
    shell:
        config['dRepdir'] + " dereplicate {params.outdir} -g {input} -pa 0.5 -sa 0.95 -p 7"
        #config['dRepdir'] + " dereplicate {output} -g {input} --completeness 0.8 -pa 0.95 \
        #--run_tax --S_algorithm fastANI -p 7 --cent_index /Users/lnmerk/centrifuge-1.0.3-beta/indices/p+h+v"

def bin_pairs():
    # We need to find all the bins, so we search through the dereplicated rep_genomes
    # We create a path: folder name directory, which we need to copy files from the folder
    # directory, and zip them into the new directory.
    bin_dir_list = Path('workflow/out/dRep_output/').glob('**/dereplicated_genomes/*.fa')

    pair_dict = {}

    for f in bin_dir_list:
        pair_dict[f.name.replace('.fa', '')] = str(f)

    return pair_dict

pair_dict = bin_pairs()

rule add_MIDAS_db:
# Add representative genomes (order 60-100) to MIDAS database.
# This finds all dereplicated bins, zips the fasta file, and places it in a directory
# named after the bin inside the MIDAS DB. It also adds an empty features file
# Lastly, it updates the speciesinfo and genomeinfo, which are needed in merge snps.
    input:
        # The first input tells us when dRep is finished running
        # The second will be all the file names of the bins, which we will
        # make individual directories for.
        expand("workflow/out/dRep_output/{cow}/data_tables/genomeInfo.csv", cow =cows),
        lambda wildcards: pair_dict.values()
    output:
        # This will tell midas snps when we've updated the database completely
        "workflow/out/complete.txt"
    run:
        # Make 2 lists for the txt files we will need to add
        genome_info_to_add = []
        species_info_to_add = []

        for f in input:
            bin = f.split('/')[-1].split('.fa')[0]
            new_dir = config['midasDBdir'] + "/" + bin
            new_path = new_dir +  "/genome.fna"

            # If genome is already in the rep genome database, delete it and replace.
            if os.path.isdir(new_dir):
                rmtree(new_dir)

            # Make the directory
            os.mkdir(new_dir)

            # This is the information we will add to the text files
            # The numbers are random, might need to come back and change
            genome_info_to_add.append(f"{bin}	{bin}	1	500000	200	{bin}")
            species_info_to_add.append(f"{bin}	{bin}	1")

            # gzip it
            with open(f, 'rb') as f_in:
                with gzip.open(new_path + ".gz", 'wb') as f_out:
                    copyfileobj(f_in, f_out)

            # For now, put the empty features file
            output_file = gzip.GzipFile(new_dir +  "/genome.features.gz", "w")

        # Now, we need to add the new bins to genome info and species info_fields
        # First, let's read in the headers since taking the set later on will scramble
        header_genome = open(config['midasDB'] + 'genome_info.txt').readline().strip()
        header_species = open(config['midasDB'] + 'species_info.txt').readline().strip()

        # Now, get every other line of the file and save as a set.
        # This will make prevent us from adding our new genomes multiple times
        curr_genome_info = set(line.strip() for i, line in enumerate(open(config['midasDB'] + 'genome_info.txt')) if i > 0)
        curr_species_info = set(line.strip() for i, line in enumerate(open(config['midasDB'] + 'species_info.txt')) if i > 0)

        # Add all the new lines to the set to avoid duplicates
        curr_genome_info.update(genome_info_to_add)
        curr_species_info.update(species_info_to_add)

        # Rewrite the two info files
        with open(config['midasDB'] + 'species_info.txt', 'w') as fp:
            # First, add the header
            fp.writelines(header_species + "\n")
            # Then, add the rest of the files.
            fp.writelines([l + "\n" for l in curr_species_info])

        with open(config['midasDB'] + 'genome_info.txt', 'w') as fp:
            fp.writelines(header_genome + "\n")
            fp.writelines([l + "\n" for l in curr_genome_info])

        with open('workflow/out/complete.txt', 'w') as fp:
            pass
