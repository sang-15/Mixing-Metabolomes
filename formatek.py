#Here is the code for the wrapper

#Retrieve data 

#TEST DATA:
#1= E. coli GCA_002861225.1
#L. crispatus: GCA_002861815.1 
#P. mirabilis: GCA_012030515.1 

from Bio import Entrez 
from Bio import SeqIO
#import requests

"""
NOTES:
    we are in the assembly database
    tying to use accession number with efetch
    
    *bad request error* due to databse name
    
    https://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.T._entrez_unique_identifiers_ui/?report=objectonly
    
    above link does not list assembly database reference...
    
TODO:
   - https://www.biostars.org/p/141581/ 
"""

url = "https://www.ncbi.nlm.nih.gov/assembly/"
gban = "GCA_002861815.1"

#option 1
#-------------------------------------------------

Entrez.email="lgonzalez7@luc.edu" #email to use entrez

#def fetchGBFile(genBankAN): #future method name and argument
handle = Entrez.efetch(db='assembly', id=gban, rettype="gb") #retrieving the file
"""
link = gban + ".gb" #preparing to save file with appropriate name and extension 
file = open(link, 'w') #open file
file.write(handle.read()) #write data to file
handle.close() #close handle
file.close() #close file

def main(): #pseudo main method 
    fetchGBFile(gban) #calling file
    
if __name__ == '__main__': #actual main method
    main()
    
#TODO:
#convert .gb file to fasta file



#option 2
#-------------------------------------------------

url = url + gban #set up url

res = requests.get(url) #use request to download from url
if res.status_code != 200: #if status code != 200 then there is an error
    print("ERROR COULD NOT GET FILE") #display error
    
with open("GCA_002861225.1", 'w') as fh: #open file
    fh.write(res.text) #write to text document to show it is working
    
for seq_record in SeqIO.parse("GCA_002861225.1", "fasta"): #loop through records
    print(seq_record.id) #print out some data
    print(len(seq_record)) #print out some data
   

"""

#OUTPUT:the output for both optioin 1 and 2 should be species name and fasta


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

#Load Silent_gene output from previous analysis
sg_out = open('sample.kegg.out.txt', 'r').read().rstrip().split('\n')

#open output file and write the first organism name
outfile = open('formatk_out.txt','w') 
outfile.write('# ' + 'organism1'+ '\n') #To-do: format the name later

#Loop through each line, identify the line with K ids and write the line to output file named 'formatk_out.txt'
for line in sg_out:
    if '\t' in line:
        outfile.write(line + '\n')
outfile.close()
