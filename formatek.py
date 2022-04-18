import glob
from Bio import Entrez 
import os
import urllib
import os.path
import json
import argparse

#Setup
#Create the output folder
os.system('mkdir $HOME/results/')

#Set path
path = os.path.expanduser('~')
path = path + '/results/'

#### path for linux: $HOME/results/ ####
#### path for python: path ####
#### new files in path for python: path + 'xxx.txt' ####


#Retrieve data 

test ='{"Escherichia coli": "GCA_002861225.1","Lactobacillus crispatus": "GCA_002861815.1"}'


parser = argparse.ArgumentParser(description="Enter in genus species and its corresponding accession number in python dictonary notation with -i (Be sure to include quotes!) and your email for entrez with -e. \n Ex. {\"Escherichia coli\": \"GCA_002861225.1\",\"Lactobacillus crispatus\": \"GCA_002861815.1\"}") 
#above line: create parser object and set description for user to learn input format
parser.add_argument('-i','--input', type=json.loads) #-i or --input set to take a json.load as argument 
parser.add_argument('-e', '--email', help='enter your email so entrez knows who you are')


# test input = python3 Final2.py -e lgonzalez7@luc.edu -i {"Escherichia coli": "GCA_002861225.1","Lactobacillus crispatus": "GCA_002861815.1"}
#args = parser.parse_args(["-e","lgonzalez7@luc.edu", "-i", test]) # read input for py
args = parser.parse_args() #reads input for linux


speciesDict = args.input #create a variable for dictionary so you don't have to call args.input


terms = list(speciesDict.values()) #user input in list form


Entrez.email=args.email #email to use entrez
#Entrez.email="lgonzalez7@luc.edu"

#1= E. coli GCA_002861225.1
#L. crispatus: GCA_002861815.1 
#P. mirabilis: GCA_012030515.1 

files = [] #empty list of file names

for item in terms: #loop through accession inputs

    handle = Entrez.esearch(db="assembly", term=item, retype="text") #search assembly database for accession inputs
    record = Entrez.read(handle) #format handle


    for id in record['IdList']: #use id(s) found in handle 
        
        # Get Assembly Summary
        esummary_handle = Entrez.esummary(db="assembly", id=id, report="full")
        esummary_record = Entrez.read(esummary_handle)
        
        # set ftp url for downloading
        url = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
        

    label = os.path.basename(url) #format ftp url for downloading
    link = os.path.join(url,label+'_genomic.fna.gz') #navigating to folder on website
    link = link.replace(os.sep, '/') #format -> replace \ with /
    print("currently downloading " + label + "...\n" ) #show progress
    
    
    urllib.request.urlretrieve(link, f'{label}.fna.gz') #command to download file
    
    files.append(label + '.fna.gz') #add file name to file list

    handle.close()

        

#Prokka
#Use Prokka to annotate inputted genome
full_species_name = list(speciesDict.keys()) #user input for genus/species name in list

split_name = [] #list for split names
for name in full_species_name: #loop through list of genus/species
    split_name.append(name.split()) #split genus and species into list
    
genus = [] #list of genus
species = [] #list of species
for spec in split_name: #loop through list of split genus/species
    genus.append(spec[0]) #add first element to genus list
    species.append(spec[1]) #add second element to species list

fasta = [] #list of fasta files
for file in files:
    os.system('gunzip ' + file) #unzip fasta files
    file = file[:-3] #remove .gz 
    fasta.append(file) #add file to list of fasta files
    
prokka_results = path + 'Prokka/' #path for Prokka results
gbk_results = path + 'GBK'
os.system('mkdir ' + gbk_results) #make directory for Prokka .gbk files

for i in range(len(fasta)): #loop over each retrieved record
    prokka_prefix = 'prokka_' + genus[i] + '_' + species[i] #prefix for prokka files
    os.system('prokka --outdir ' + prokka_results + genus[i] + '_' + species[i] + ' --prefix prokka_' + genus[i] + '_' + species[i] + ' --genus ' + genus[i] + ' --species ' + species[i] + ' ' + fasta[i]) #Prokka command
    os.system('mv ' + prokka_results + genus[i] + '_' + species[i] + '/' + prokka_prefix + '.gbk ' + gbk_results) #move .gbk to folder for batch script


#SilentGene
#retrieve and use prokka2kegg_batch to convert to K IDs

if not os.path.exists(path + 'prokka2kegg_batch.py'): #check if script and db already exist
    os.system('wget https://raw.githubusercontent.com/SilentGene/Bio-py/master/prokka2kegg/prokka2kegg_batch.py -P ' + path) #download prokka2kegg_batch script
    os.system('wget https://github.com/SilentGene/Bio-py/blob/master/prokka2kegg/idmapping_KO.tab.gz?raw=true -P ' + path) #download prokka2kegg database

os.system('python3 ' + path + 'prokka2kegg_batch.py -i ' + gbk_results + ' -o ' + path + '2kegg/ -d ' + path + 'idmapping_KO.tab.gz?raw=true') #convert to K id's

#Aggrgation

#Query data for KEGG MAPPER
#Two-column dataset with K numbers in the second column, optionally preceded by the user's identifiers in the first column.
#This is consistent with the output files of automatic annotation servers, BlastKOALA, GhostKOALA, KofamKOALA and KAAS.
#The dataset may contain lists of K numbers for multiple organisms, each list preceded by the comment line starting with # and optional color specification.

# Open output file
outfile = open(path+'formatk_out.txt','w') 

#Read each prokka outputs 
for i in glob.glob("$HOME/results/2kegg/*.gbk.ko.out"):

    #Load Silent_gene output
    sg_out = open(i, 'r').read().rstrip().split('\n')

    #Set organism name
    org = i[7:-11]
    
    #write the first organism name
    outfile.write('# ' + org + '\n') #To-do: format the name later

    #Loop through each line, identify the line with K ids and write the line to output file named 'formatk_out.txt'
    for line in sg_out:
        if '\t' in line:
            outfile.write(line + '\n')
            
outfile.close()
