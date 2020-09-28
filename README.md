# Supp.Material.OA.eDNA
Data and scripts needed to process the output of the demultiplexer for dada2 for the 7 MiSeq runs used in the manuscript by Gallego et al, 2019.

The workflow starts with the output from 7 independent runs of the pipeline (https://github.com/ramongallego/demultiplexer_for_dada2) on 7 MiSeq runs. To replicate that step, the user needs to download the raw fastq [files](https://drive.google.com/drive/folders/1XtuDMdlnk9V2acNiPJo5T6O1Tp0KjaCz?usp=sharing) and install the pipeline and dependencies. Each run of the pipeline produces a folder that includes an abundance table (ASV_table.csv), a Hash key to keep record of the unique sequences found (Hash_key.csv) and the metadata associated with that run. 

So the input for this workflow consists on 7 folders with the output of the 7 Miseq runs. They are all inside the folder `/1.Sequencing_runs_output`.


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
  
  -- insect
  
  -- RColorBrewer
  
  -- rstanarm
  
  -- fitdistrplus
  
  -- ggridges
  
  -- ggforce
  
  -- cowplot
  
  -- concaveman
  
  -- lubridate
  
  -- ggrepel
  
  -- factoextra
  
  -- nnet
  
If you want to run the occupancy modelling, then you also need
  
 - rjags
 - A JAGS installation in your system, instructions here http://mcmc-jags.sourceforge.net/
 
 # Installation
 
 I think the easiest way of using this project is by downloading the repository directly into Rstudio. On the right top corner of Rstudio click on `New project` -> From GIT -> and the URL is https://github.com/ramongallego/eDNA.and.Ocean.Acidification.Gallego.et.al.2020.
 
 Alternatively, click download on the github page and all the files you need will be downloaded on a zip file

# Step 1: From MiSeq outputs to a clean, decontaminated dataset

## Running the script

The first script to be run is the Markdown file:`Denoising.all.runs.Rmd` from the folder Scripts. It will run smoothly until step *Cleaning Process 4* where it needs to import the file "Occ.fate.csv". You will find one file with that name now - run by me - but if you want, you can run the occupancy yourself once you have run the script to the point it generates the file "Cleaning.before.Occ.model.rds"
  
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

Here we used two different methods to obtain a taxonomic annotation:

- Phylogenetic tree placement: Using the package `insect`, an ID is returned by placing the query sequence in a curated, reference phylogenetic tree. 
- CRUX & Bowtie2: Using a reference database made through *in-silico* PCR we used Anacapa's classifier to retrieve the most probable taxonomical annotation. See DOI: 10.1111/2041-210X.13214. This was performed 

CRUX database had a higher coverage and thus retrieved more positive IDs. 

The CRUX-COI database is still embargoed and won't be released until May 2020.

## Taxonomic metadata

For each taxa with a taxonomic annotation to at least the family level, we proceeded to annotate their mineralization stage (both in larvae and adult forms), their lifestyle (benthos, Plankton, None), trophic level. 

The file with the Annotation is "all.taxonomy.20190130.csv". The process to obtain that file is in the Rmarkdown `taxonomy.biom.food.env.Rmd`

## Higher taxonomy

Using NCBI's taxonomy, the file with the annotation is "higher_taxonomy.csv"


# Environmental data

The environmental dataset is on the file "env.data.updated.csv"

## All taxa annotated

Combining previous files: environment, taxonomy and ASV - the output is "Combined_Biol_Env_Plankton.csv". You can see the process and replicate it using the rmarkdown script `From.asv.to.data.plankton.Rmd`

# Final Step: Manuscript Main Analysis

This script relies on the previous outputs: the abundance table (`Input/Combined_Biol_Env_Plankton.csv`), the taxonomy and meta_annotation file (`Input/higher_taxonomy.csv`) and the models.

The script uses the Stan models for the multinomial regression - the models are already calculated on standarized values of pH and temperature. They are loaded within the object models at the beginning of the script. The code that leads to the generation of these models is in the file `furtherModeling.R`, which needs the file `Prepare.data.for.stan.Rmd` for precisely that.

The manuscript also relies in a simulated data of pH and temperature for the year 2095 for the Hood Canal and San Juan Island. The simulation occurs within the main script, although it requires functions and estimations made in the file `Future_climate_conditions.r`

The Script that runs the analysis is `Manuscript.Main.Analysis.Rmd`

The output of this script includes all the Figures that go in the Manuscript, the calculations behind other stats reported in the Main Manuscript and the Supplementary Figures
