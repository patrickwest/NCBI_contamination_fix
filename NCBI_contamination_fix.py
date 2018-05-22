#Trim sequences based on an NCBI sequence error/contamination report
#Patrick West 2018

from Bio import SeqIO
import argparse

def main():
    contamination_fh = open(args.i)
    genome_fh = open(args.f)
    out_fh = open(args.o,'w')

    seqs_to_trim = parse_contamination(contamination_fh)
    trim_sequences(genome_fh, seqs_to_trim, out_fh)

    out_fh.close()

def parse_contamination(fh):
    to_trim = []
    #Contamination file has a huge header, so have this function to skip to only the relevant lines
    trim = False
    for line in fh:
        if 'Trim:' in line:
            trim = True
            #Has secondary header after 'Trim:' line, so skip it
            next(fh)
        elif trim == True:
            sline = line.strip().split()
            #contamination file has blank line at the end that causes errors, skip it by checking how many fields after spliting the line
            if len(sline) >= 3:
                to_trim.append(sline)
    return to_trim

def trim_sequences(fh, to_trim, out_fh):
    for record in SeqIO.parse(fh, "fasta"):
        for trim in to_trim:
            if trim[0] == record.id:
                trim_start, trim_end = trim[2].split('..')
                #trim
                record.seq = str(record.seq)[:int(trim_start)-1] + str(record.seq)[int(trim_end)+1:]
                #Remove any Ns on the end of sequence
                while record.seq[0] == 'N':
                    record.seq = str(record.seq)[1:]
                while record.seq[-1] == 'N':
                    record.seq = str(record.seq)[:-1]
        #check length. NCBI requires trimmed seqeuences to still be longer than 200bp
        if len(str(record.seq)) > 200:
            out_fh.write('>' + record.description + '\n')
            out_fh.write(str(record.seq) + '\n')

if __name__ == '__main__':
    # argument parser
    parser = argparse.ArgumentParser(description='Modifies a genome file for NCBI submission based upon errors or issues in the "contamination.txt" file.')
    # input/output
    parser.add_argument('-i', required = True, help = 'NCBI contamination text file')
    parser.add_argument('-f', required = True, help = 'genome fasta file')
    parser.add_argument('-o', required = True, help = 'output file name for new genome fasta')
    args = parser.parse_args()

    main()
