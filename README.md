# Overview
This is the repository created for COMP483 group project #2: <br />
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
### Default setting
The following code will run the wrapper with default setting by giving GenBank assembly accession (eg. [GCA_002861225.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_002861225.1))
```
python3 formatek.py -i '{\"Escherichia coli\": \"GCA_002861225.1\"}' -e useremail
```

### User specified FASTA
The following code will run the wrapper with user specified FASTA
```
python3 formatek.py -i '{\"Escherichia coli\": \"LOCATION\.FASTA"}' -e useremail
```


## Output
The wrapper will generate a 'results' folder under '$HOME/' directory, and the folder contains the following: <br />

- (Default) Downloaded assemble genome data from NCBI
- Prokka folder <br /> Prokka generated outputs 
- gbk_results folder <br /> prokka2kegg generated outputs 
- formatek_out.txt <br />
Final gene list required by KEGG MAPPER


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
