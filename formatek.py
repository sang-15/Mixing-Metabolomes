import glob
from Bio import Entrez 
#from Bio import SeqIO
import os
import urllib


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
Entrez.email="lgonzalez7@luc.edu" #email to use entrez

#1= E. coli GCA_002861225.1
#L. crispatus: GCA_002861815.1 
#P. mirabilis: GCA_012030515.1 

terms= ["GCA_002861225.1","GCA_002861815.1 ","GCA_012030515.1"] #user input in list form
files = [] #list of files to be produced during downloads

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
    print("currently downloading " + label + "..." ) #show progress
    
    
    urllib.request.urlretrieve(link, f'{label}.fna.gz') #command to download file
    
    files.append(label + '.fna.gz') #add file name to file list

    handle.close()
        

#Prokka
#Use Prokka to annotate inputted genome
fasta = [] #list of fasta files
for file in files:
    os.system('gunzip ' + file) #unzip fasta files
    file = file[:-3] #remove .gz 
    fasta.append(file) #add file to list of fasta files
    
genus = ['Escherichia', 'Proteus', 'Lactobacillus'] #CHANGE: list of genus
species = ['coli', 'mirabilia', 'crispatus'] #CHANGE: list of species
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
for i in glob.glob("$HOME/results/GBK/*.gbk.ko.out"):

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
