#Here is the code for the wrapper








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
