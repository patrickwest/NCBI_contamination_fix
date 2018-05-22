# NCBI_contamination_fix

Following a NCBI WGS genome submission, 
Trim sequences based on an NCBI sequence error/contamination report

### Example Usage ###
```
python NCBI_contamination_fix.py -i Contamination.txt -f genome.fa -o trimmed_genome.fa
```
See 'Example_Contamination.txt' as an example of an NCBI Contamination file

### Prerequisites ###
biopython
