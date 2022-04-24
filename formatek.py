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

from Bio import Entrez 
import json
import os
import urllib
import argparse
import os.path

test ='{"Escherichia coli": "GCA_002861225.1","Lactobacillus crispatus": "GCA_002861815.1"}'
homeDir = os.path.expanduser('~')

newpath = homeDir + '/results/downloads' 
if not os.path.exists(newpath):
    os.system('mkdir ' + newpath)

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
if Entrez.email == None:
	Entrez.email="ylin22@luc.edu"


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
    
        
    urllib.request.urlretrieve(link, homeDir + f'/results/downloads/{label}.fna.gz') #command to download file linux
    

    handle.close()

      

#Prokka
#Use Prokka to annotate inputted genome
full_species_name = list(speciesDict.keys()) #user input for genus/species name in list

split_name = dict() #dictionary for split names
for name in speciesDict.keys(): #loop through list of names
    split_name[speciesDict[name]] = name.split() #split genus and species into list and add to dict
 
fasta = dict() #list of fasta files
for accession in files.keys(): #loop through file paths
    if '.gz' in files[accession]: #if file is zipped
        os.system('gunzip ' + files[accession]) #unzip fasta files
        file = files[accession][:-3] #remove .gz 
        fasta[accession] = file #add file to list of fasta files
    else:
        fasta[accession] = files[accession] #if not zipped, append file path without further action
    
prokka_results = path + 'Prokka/' #path for Prokka results
gbk_results = path + 'GBK'
os.system('mkdir ' + gbk_results) #make directory for Prokka .gbk files

for entry in fasta.keys(): #loop over each accession
    genus = split_name[entry][0] #retrieve genus name
    species = split_name[entry][1] #retrieve species name
    prokka_prefix = 'prokka_' + genus + '_' + species #prefix for prokka files    
    os.system('prokka --outdir ' + prokka_results + genus + '_' + species + ' --prefix prokka_' + genus + '_' + species + ' --genus ' + genus + ' --species ' + species + ' ' + fasta[entry]) #Prokka command
    os.system('mv ' + prokka_results + genus + '_' + species + '/' + prokka_prefix + '.gbk ' + gbk_results) #move .gbk to folder for batch script


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
#The dataset may contain lists of K numbers for multiple organisms, each list preceded by the comment line starting with #.

# Open output file
outfile = open(path+'formatk_out.txt','w') 

#Read each prokka outputs 
for i in glob.glob(path+"2kegg/*.gbk.ko.out"):

    #Load Silent_gene output
    sg_out = open(i, 'r').read().rstrip().split('\n')

    #Set organism name
    org = i.split('/')
    org = org[-1]
    org = org[7:-11]
    
    #write each organism name
    outfile.write('# ' + org + '\n') 

    #Loop through each line, identify the line with K ids and write the line to output file named 'formatk_out.txt'
    for line in sg_out:
        if '\t' in line:
            outfile.write(line + '\n')
            
outfile.close()
