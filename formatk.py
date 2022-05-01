import glob
from Bio import Entrez 
import os
import urllib
import os.path
import json
import argparse
import urllib.request
from urllib.error import HTTPError

##Setup
#Set flags for wrapper input and implement each flag to the wrapper 
parser = argparse.ArgumentParser(description='Thank you for using formatk! Please refer to READ.ME or our GitHub page for more detailed user information: https://github.com/sang-15/Mixing-Metabolomes.') 
#above line: create parser object and set description for user
parser.add_argument('-i','--input',
                    help='Enter in accession number or local assembled genome file location in a python dictonary notation. \n Ex. {\"GCA_002861225.1\":\"Escherichia coli\", \"DIR.fasta\": \"Lactobacillus crispatus\"} (Be sure to include quotes outside the brackets as well!)') 
#-i or --input set to take a json.load as argue4xiyment 
parser.add_argument('-e', '--email', default = 'ylin22@luc.edu', help='Enter your email so NCBI Entrez API knows who you are')
parser.add_argument('-o', '--output', default = 'formatkresults', help ='Name for output folder')

#reads input for linux
args = parser.parse_args() 

#Check input format
#End the program if the input format is wrong
try:
    a_json = json.loads(args.input)
except:
    print("ERROR: Something is wrong with the input format. Please try again! Use -h for help!")
    print("ERROR: End program.")
    os._exit(0)

#Check input duplicates
#End the program if there is duplicated local assembled genome file or NCBI accession number 
n = []
t = str(args.input).replace('{','').replace('}','').replace(' ','').replace('"','').split(',')
for i in t:
    n.append(i.split(':')[0])

if len(n) != len(set(n)):
    print("ERROR: Duplicated local assembled genome file or NCBI accession number found in input! Please remove them and try again!")
    print("ERROR: End program.")
    os._exit(0)

#Check input sample number (maximum 10 allowed)
#End the program if there are more than 10 samples
if len(a_json) > 10:
    print("ERROR: Too many samples! We can only take upto 10 samples. Please try again!")
    print("ERROR: End program.")
    os._exit(0)

    
#Set path
outputdir = args.output
path = os.path.expanduser('~')
path = path + '/' + outputdir + '/'

## Check the availability of the output folder
#End the program if the outfolder already exsit and print out info
if os.path.exists(path):
    print("ERROR: Looks like you already have folder " + path +", plese use -o to rename your output folder and try again!")
    print("ERROR: End program.")
    os._exit(0)

#Set output folder
os.system('mkdir ' + path)

#Set up email to use entrez
Entrez.email=args.email 


##Retrieve data
speciesDict = json.loads(args.input) #create a variable for dictionary 
terms = list(speciesDict.keys()) #user input in list form
files = dict() #empty dict of file names

for item in terms: #loop through accession inputs

    #For user specified input, check if the file exsit, if not print error and exit program
    if '.fasta' in item or '.fna' in item: 
        if not os.path.exists(os.path.expanduser(item)):
            print('ERROR: We CANNOT locate your file ' + item + ', plese check your input and try again!')
            print("ERROR: End program.")
            os._exit(0)
            
    else: # if user is not user supplied file
        newpath = path + 'downloads/' #path to download files
        os.system('mkdir ' + newpath) #create path to download files
        break
      
for item in terms: #loop through accession inputs

    if '.fasta' in item or '.fna' in item: # if term is user supplied file
        print('User supplied file: ' + item)
        files[item] = item #append path to files list directly
    
    else:

        handle = Entrez.esearch(db="assembly", term=item, retype="text") #search assembly database for accession inputs
        record = Entrez.read(handle) #format handle

        #Check if accession number input can be found in NCBI
        #End program if the accession number is wrong
        if 'ErrorList' in record:
            print('ERROR: We CANNOT locate your accession number ' + item + ' in NCBI assembly database, plese check your input and try again!')
            print("ERROR: End program")
            os._exit(0)

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
        
       
        urllib.request.urlretrieve(link, newpath + f'{label}.fna.gz') #command to download file linux       
        print("currently downloading " + label + "...\n" ) #show progress
        files[item] = newpath + label + '.fna.gz' #add file name to file dict
    
        handle.close()

if os.path.exists('esummary_assembly.dtd'):
    os.system('mv esummary_assembly.dtd ' + path) #Move the files to output file

      
#Prokka
#Use Prokka to annotate inputted genome

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
gbk_results = path + 'GBK/'
os.system('mkdir ' + gbk_results) #make directory for Prokka .gbk files


for entry, name in fasta.items(): #loop over each file
    newname = name.split('/')[-1]
    newname = newname.split('.')[0] # Add file name
    genus = split_name[entry][0] #retrieve genus name
    species = split_name[entry][1] #retrieve species name
    suffix = genus + '_' + species + '_' + newname #genus_species_ID
    prokka_prefix = 'prokka_' + suffix #prefix for prokka files    
    os.system('prokka --usegenus --outdir ' + prokka_results + suffix + ' --prefix prokka_' + suffix + ' --genus ' + genus + ' --species ' + species + ' ' + fasta[entry]) #Prokka command
    os.system('mv ' + prokka_results + suffix + '/' + prokka_prefix + '.gbk ' + gbk_results) #move .gbk to folder for batch script




#SilentGene
#retrieve and use prokka2kegg_batch to convert to K IDs

if not os.path.exists(path + 'prokka2kegg_batch.py'): #check if script and db already exist
    os.system('wget https://raw.githubusercontent.com/SilentGene/Bio-py/master/prokka2kegg/prokka2kegg_batch.py -P ' + path) #download prokka2kegg_batch script
    os.system('wget https://github.com/SilentGene/Bio-py/blob/master/prokka2kegg/idmapping_KO.tab.gz?raw=true -P ' + path) #download prokka2kegg database

os.system('python3 ' + path + 'prokka2kegg_batch.py -i ' + gbk_results + ' -o ' + path + '2kegg/ -d ' + path + 'idmapping_KO.tab.gz?raw=true') #convert to K id's

#Aggrgation

#Query data for KEGG MAPPER

#formatk_out.txt: 
#Two-column dataset with K numbers in the second column, optionally preceded by the user's identifiers in the first column.
#This is consistent with the output files of automatic annotation servers, BlastKOALA, GhostKOALA, KofamKOALA and KAAS.
#The dataset may contain lists of K numbers for multiple organisms, each list preceded by the comment line starting with #.

#formatk_out_order.txt:
#The final order of the KEGG output which corresponding to the input order

# Open output files
outfile = open(path+'formatk_out.txt','w') 
of = open(path+'formatk_out_order.txt','w')

count = 1

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
    
    #Provide final order information corrsponding to the initial input
    of.write(org + 'is used as input number ' + str(count) + ' to KEGG.' + '\n')
    count += 1

    #Loop through each line, identify the line with K ids and write the line to output file named 'formatk_out.txt'
    for line in sg_out:
        if '\t' in line:
            outfile.write(line + '\n')
                
outfile.close()
