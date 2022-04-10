#Here is the code for the wrapper

#Importing libraries
import glob


#Retrieve data 
from Bio import Entrez 
#from Bio import SeqIO

import os
import urllib

Entrez.email="lgonzalez7@luc.edu" #email to use entrez

#1= E. coli GCA_002861225.1
#L. crispatus: GCA_002861815.1 
#P. mirabilis: GCA_012030515.1 

terms= ["GCA_002861225.1","GCA_002861815.1 ","GCA_012030515.1"] #user input in list form

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
    

    handle.close()
        

#Prokka + Silent Gene
#Use Prokka to annotate inputted genome
#Use Prokka output in Silent Gene script to get K numbers

#results = output dir, genus = genus name from input, fasta = genome from input (update later)
os.system('prokka --outdir $HOME/' + results + ' --prefix prokka --genus ' + genus + ' ' + fasta) #Prokka command

os.system('python3 prokka2kegg.py -i prokka.gbk -d idmapping_KO.tab.gz -o sample.kegg.out.txt') #Silent Gene command

#Aggrgation

#Query data for KEGG MAPPER
#Two-column dataset with K numbers in the second column, optionally preceded by the user's identifiers in the first column.
#This is consistent with the output files of automatic annotation servers, BlastKOALA, GhostKOALA, KofamKOALA and KAAS.
#The dataset may contain lists of K numbers for multiple organisms, each list preceded by the comment line starting with # and optional color specification.

# Open output file
outfile = open('formatk_out.txt','w') 

#Read each prokka outputs 
for i in glob.glob("*.gbk.ko.out"):

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
