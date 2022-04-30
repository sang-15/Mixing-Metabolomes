# Overview
This is the repository created for COMP483 Project by group 2: <br />
A Python wrapper to provide a convienacne way to generate K identifiers required for [KEGG Mapper](https://www.genome.jp/kegg/mapper/reconstruct.html), by giving assembled genome sequences. 

## Perequisites
- **[Python3](https://www.python.org/)** <br />
  - **[Biopython](https://biopython.org/)** <br />

- **[SRA-Toolkit](https://www.ncbi.nlm.nih.gov/sra)** <br />
For downloading FASTA from NCBI <br /> 

- **[Prokka](https://github.com/tseemann/prokka)** <br />
For rapid prokaryotic genome annotation <br />

## Included within wrapper
This wrapper uses **[SilentGene - prokka2kegg_batch](https://github.com/SilentGene/Bio-py/tree/master/prokka2kegg)** to automatically convert annotated gene from Prokka's (.gtk file) to K ids required. <br />
The python script and database required for prokka2kegg are included in the wrapper, user does not need to download or install anything for this part.

## Running wrapper
**-i**: (required) Specifiey the  FASTA source and Taxonomy information for each sample in a json format<br />
**-e**: (oprional) Email used for Biopython <br />
**-o**: (optional) An optional flag for user to name the output folder for each run. Defalt will give output in a folder named 'formatkresults'.

### GenBank assembly accession

The following code will run the wrapper by giving GenBank assembly accession <br />

- Runing 1 sample: [GCA_002861225.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_002861225.1)
```
python3 formatk.py -i '{"GCA_002861225.1": "Escherichia coli"}' 
```
- Runing multiple samples: [GCA_002861225.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_002861225.1) and [GCA_002861815.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_002861815.1)
```
python3 formatk.py \
-i '{"GCA_002861225.1": "Escherichia coli", "GCA_002861815.1": "Lactobacillus crispatus"}'
```

### User specified FASTA

The following code will run the wrapper with user specified FASTA<br />

- Runing 1 sample: Escherichia_coli.FASTA
```
python3 formatk.py -i '{"DIR/Escherichia_coli.FASTA": "Escherichia coli"}' 
```
- Runing multiple samples: Escherichia_coli.FASTA and Lactobacillus_crispatus.FASTA
```
python3 formatk.py \
-i '{"DIR/Escherichia_coli.FASTA": "Escherichia coli", "DIR/Lactobacillus_crispatus.FASTA": "Lactobacillus crispatus"}' 
```

### Mixing both GenBank assembly accession and user specified FASTA
The following code will run the wrapper with mixed GenBank assembly accession and user specified FASTA<br />
- Runing multiple samples: [GCA_002861225.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_002861225.1) and Lactobacillus_crispatus.FASTA
```
python3 formatk.py \
-i '{"GCA_002861225.1": "Escherichia coli", "DIR/Lactobacillus_crispatus.FASTA": "Lactobacillus crispatus"}' 
```
### With optional flags -o and -e
The following code will run the wrapper by giving GenBank assembly accession ([GCA_002861225.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_002861225.1)) and give the results in a folder named **'TESTRUN'** under '$HOME/' directory with your own email
```
python3 formatk.py -i '{"GCA_002861225.1": "Escherichia coli"}' -e useremail -o TESTRUN 
```

## Output
We rename each sample by the end of the process based on sample taxonomy information and it's input file name or GenBank assembly accession as 'genus_species_filename' or 'genus_species_GenBank_ID'. <br />

The wrapper will generate a 'formatkresults' folder under '$HOME/' directory (or user specified folder if '-o' argument is supplied), and the folder contains the following: <br />

- **formatk_out.txt** <br />
Final gene list required by KEGG MAPPER
- **formatk_out_order.txt** <br />
The order of each sample in formatk_out.txt and KEGG 
- Others
  - Downloaded assemble genome data from NCBI (if applicable)
  - Prokka folder <br /> 
    Prokka generated results for each entry
  - GBK folder <br /> 
    Prokka generated gbk files output for each entry
   - 2kegg <br />
    SlientGene generated results for each entry


## Test data
- E. coli: [GCA_002861225.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_002861225.1) 
- L. crispatus: [GCA_002861815.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_002861815.1) 
- P. mirabilis: [GCA_012030515.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_012030515.1) 

## [Color cods for KEGG MAPPER](https://www.genome.jp/kegg/kegg1c.html#mapping)

- **Pathway mapping of two organisms (global/overview maps)** <br />
Organism 1	      	![#00cc33](https://via.placeholder.com/15/00cc33/000000?text=+) `#00cc33` <br />
Organism 2	      	![#ff3366](https://via.placeholder.com/15/ff3366/000000?text=+) `#ff3366`<br />
Organisms 1 and 2	      	![#3366ff](https://via.placeholder.com/15/3366ff/000000?text=+) `#3366ff`<br />

- **Pathway mapping of multiple organisms (regular maps)** <br />
Organism 1	      	![#bfffbf](https://via.placeholder.com/15/bfffbf/000000?text=+) `#bfffbf`<br />
Organism 2	      	![#ffbbcc](https://via.placeholder.com/15/ffbbcc/000000?text=+) `#ffbbcc`<br />
Organism 3	      	![#bbccff](https://via.placeholder.com/15/bbccff/000000?text=+) `#bbccff`<br />
Organism 4	      	![#cfffcf](https://via.placeholder.com/15/cfffcf/000000?text=+) `#cfffcf`<br />
Organism 5	      	![#ffcfef](https://via.placeholder.com/15/ffcfef/000000?text=+) `#ffcfef`<br />
Organism 6	      	![#cfefff](https://via.placeholder.com/15/cfefff/000000?text=+) `#cfefff`<br />
Organism 7	      	![#dfefcf](https://via.placeholder.com/15/dfefcf/000000?text=+) `#dfefcf`<br />
Organism 8	      	![#ffefcc](https://via.placeholder.com/15/ffefcc/000000?text=+) `#ffefcc`<br />
Organism 9	      	![#dfccff](https://via.placeholder.com/15/dfccff/000000?text=+) `#dfccff`<br />
Organism 10	      	![#dfdfcc](https://via.placeholder.com/15/dfdfcc/000000?text=+) `#dfdfcc`<br />
