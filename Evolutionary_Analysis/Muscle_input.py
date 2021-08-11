#this script prepares input for MUSCLE to perform multiple sequence alignment on single copy orthogroups of our 4 species
# it uses the list of single copy orthogroups got from orthofinder
# and the protein fasta files from each of the species
# we had 4949 SCOs
# so it creates 4949 files. Each file as 1 gene from each species with the respective AA sequences


#function to create dictionaries
def parse_diction(file_name):
    diction = {}
    with open(file_name,"r") as file: #opens the file while needed
        file = file.readlines()
        for y in range(len(file)):
            if file[y].startswith(">"):
                key = file[y].strip()
                value = file[y+1].strip()
                diction[key]= value
            else:
                pass
    return diction

#call function to create dictionaries
Dictio_dc = parse_diction("Dictyo_SCOs.fa")
Myrio_dc = parse_diction("Myriotrichia_SCOs.fa")
Desma_dc = parse_diction("Desmarestia_SCOs.fa")  
Saccho_dc = parse_diction("Sacchoriza_SCOs.fa")


#read single copy orthogroups file
scos = open("OrthologsIDS_test.txt","r").readlines() #read SCOs and put into a list


for line in scos:
    geneids = line.split(" ")
    file = open("muscle_input/" + geneids[0] + ".fa","w")


#write files
    file.write(">" + geneids[1] + "\n" + Dictio_dc[">" + geneids[1]] + "\n")
    file.write(">" + geneids[2] + "\n" + Myrio_dc[">" + geneids[2]] + "\n")
    file.write(">" + geneids[3] + "\n" + Desma_dc[">" + geneids[3]] + "\n")
    file.write(">" + geneids[4] + Saccho_dc[">" + geneids[4].strip()] + "\n")
    
    file.close()

    #Tara....
