
import glob
from Bio import Entrez 
import os
import urllib
import os.path
import json
import argparse

##Setup
#Set flags for wrapper input and implement each flag to the wrapper 
parser = argparse.ArgumentParser(description='Thank you for using formatk! Please refer to READ.ME or our GitHub page for more detailed user information: https://github.com/sang-15/Mixing-Metabolomes.') 
#above line: create parser object and set description for user
parser.add_argument('-i','--input', type=json.loads, 
                    help='Enter in accession number or data location in python dictonary notation. (Be sure to include quotes outside the brackets as well!) \n Ex. {: \"GCA_002861225.1\":\"Escherichia coli\", \"GCA_002861815.1\": \"Lactobacillus crispatus\"}') 
#-i or --input set to take a json.load as argument 
parser.add_argument('-e', '--email', default = 'ylin22@luc.edu', help='Enter your email so entrez knows who you are')
parser.add_argument('-o', '--output', default = 'formatkresults', help ='Name for output folder')

args = parser.parse_args() #reads input for linux
outputdir = args.output

#Create the output folder
os.system('mkdir ' + outputdir)

#Set path
path = os.path.expanduser('~')
path = path + '/' + outputdir + '/'

newpath = path + 'downloads/' #path to download files
if not os.path.exists(newpath):
    os.system('mkdir ' + newpath) #create path to download files

    
##Retrieve data
speciesDict = args.input #create a variable for dictionary so you don't have to call args.input
terms = list(speciesDict.keys()) #user input in list form
Entrez.email=args.email #email to use entrez
files = dict() #empty dict of file names

for item in terms: #loop through accession inputs

    if '.fasta' in item or '.fna' in item: # if term is user supplied file
        print('User supplied file: ' + item)
        files[item] = item #append path to files list directly
    
    else:

        handle = Entrez.esearch(db="assembly", term=item, retype="text") #search assembly database for accession inputs
        record = Entrez.read(handle) #format handle


        for id in record['IdList']: #use id(s) found in handle 
        
            # Get Assembly Summary
            esummary_handle = Entrez.esummary(db="assembly", id=id, report="full")
            esummary_record = Entrez.read(esummary_handle)
        
            # set ftp url for downloading
            url = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
    
        print("this is the url variable: " + url)

        label = os.path.basename(url) #format ftp url for downloading
        link = os.path.join(url,label+'_genomic.fna.gz') #navigating to folder on website
        link = link.replace(os.sep, '/') #format -> replace \ with /
        print("currently downloading " + label + "...\n" ) #show progress
    
       
        urllib.request.urlretrieve(link, newpath + f'{label}.fna.gz') #command to download file linux       
        
        files[item] = newpath + label + '.fna.gz' #add file name to file dict
    
        handle.close()

      
#Prokka
#Use Prokka to annotate inputted genome
#full_species_name = list(speciesDict.values()) #user input for genus/species name in list

split_name = dict() #dictionary for split names
for i in speciesDict.keys(): #loop through list of names
    split_name[i] = speciesDict[i].split() #split genus and species into list and add to dict
 
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
    suffix = genus + '_' + species + '_' + entry #genus_species_HHMMS
    prokka_prefix = 'prokka_' + suffix #prefix for prokka files    
    os.system('prokka --outdir ' + prokka_results + suffix + ' --prefix prokka_' + suffix + ' --genus ' + genus + ' --species ' + species + ' ' + fasta[entry]) #Prokka command
    os.system('mv ' + prokka_results + suffix + '/' + prokka_prefix + '.gbk ' + gbk_results) #move .gbk to folder for batch script


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
