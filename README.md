# Supp.Material.OA.eDNA
Data and scripts needed to process the output of the demultiplexer for dada2 for the 7 MiSeq runs used in the manuscript by Gallego et al, 2019.

The workflow starts with the output from 7 independent runs of the pipeline (https://github.com/ramongallego/demultiplexer_for_dada2) on 7 MiSeq runs. To replicate that step, the user needs to download the raw fastq files and install the pipeline and dependencies. Each run of the pipeline produces a folder that includes an abundance table (ASV_table.csv), a Hash key to keep record of the unique sequences found (Hash_key.csv) and the metadata associated with that run. 

So the input for this workflow consists on 7 folders with the output of the 7 Miseq runs. They are all inside the folder /1.Sequencing_runs_output.


## Dependencies

In order to make this work, you need to have the following up and running:
- R (v3.6.0 or above)
- RStudio (Version 1.2.1335 or above)
- The following packages:

  -- tidyverse
  -- MASS
  -- vegan
  -- proxy
  -- reshape2
  -- sequinr
  
If you want to run the occupancy modelling, then you also need
  
  -- rjags
 - A JAGS installation in your system, instructions here http://mcmc-jags.sourceforge.net/

# Step 1: Clean up

## Running the script

The first script to be run is the Markdown file: Denoising.all.runs.Rmd from the folder Scripts. It will run smoothly until step *Cleaning Process 4* where it needs to import the file "Occ.fate.csv". You will find one file with that name now - run by me - but if you want, you can run the occupancy yourself once you have run the script to the point it generates the file "Cleaning.before.Occ.model.rds"
  
## Occupancy modelling

The script uses Occupancy modeling to determine if the presence of an Sequence in a sample is real or not. To do so, there is a second script (Rjags.tunning.Rmd). It takes a lot of computer power / time. So be patient. Or skip it - and trust the iteration we used.

## Output (s)

All your hard work will be paid with:

- An .html file with your progress, and some graphs
- The output files with the clean dataset (an Abundance table - ASV_table_all_together.csv -  and a hash key -Hash_Key_all_together.csv)
- The summary file of all the clean-up operations, with the number of samples, ASVs and reads that passed each of the cleaning steps.


- A summary file of the trimmings

# Step 2: Taxonomic Annotation

The taxonomic annotation of DNA sequences is likely to change in the future: "From its inception, GenBank has grown exponentially, and continues to do so with the number of sequence records doubling approximately every 35 months". I see this not as a problem, but a potential to revisit this dataset in the future. Hopefully, these scripts will be able to cope with a new annotation.


  
