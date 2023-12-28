import argparse
import os

def seperate_multi_fasta(filename,output_folder) :

    if os.path.exists(output_folder) == False :
        os.mkdir(output_folder)

    f=open(filename,"r")
    opened = False
    for line in f :
        if(line[0] == ">") :
            if(opened) :
                fasta.close()
            opened = True
            fasta=open("%s/%s.fa" % (output_folder,line[1:].rstrip()), "w")
            print(line[1:].rstrip())
        fasta.write(line)
    fasta.close()

def main() :
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",help="path of multiple sequence fasta")
    parser.add_argument("-o", "--output",help="path of multiple sequence fasta")
    args = parser.parse_args()

    seperate_multi_fasta(args.input,args.output)

if __name__ == '__main__' :
    main()