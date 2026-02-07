#Retriving sequence from uniprot database in fasta formate.
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
print("All successful")

                #2. Sequence Quality Analysis

record = SeqIO.read("NS1.fasta", "fasta")
print("The Sequence ID:", record.id)
print("The Sequence Description:", record.description)
print("Sequence:", record.seq)
print("Type of the seq record:", type(record))

print("The Length of the Amino acid:", len(record))
 # Amino acid composition
sequence = str(record.seq)
aa_count = {}
for aa in sequence:
    aa_count[aa] = aa_count.get(aa, 0) + 1

print("Amino acid composition:")
for aa, count in sorted(aa_count.items()):
    print(f"{aa}: {count}")

               #3. Sequence Filtering & Validation

if len(record.seq) < 50:
    print("Too short sequence")
else:
    print("Sequence is fine")

                #4. Homology Search (BLAST)

result_handle = NCBIWWW.qblast("blastp", "nr", record.seq)
with open("blast_results.xml", "w") as out_b:
    out_b.write(result_handle.read())
    
print("BLAST search complete. Results saved to blast_results.xml")
with open("blast_results.xml") as b:
    blast_record = NCBIXML.read(b)

top_hit = blast_record.alignments[0]
print("Top BLAST hit:")
print(top_hit.hit_def)

top_hsp = top_hit.hsps[0]

                #closest homolog

print("Closest homolog:")
print("Protein:", top_hit.hit_def)
print("Length:", top_hit.length)
print("E-value:", top_hsp.expect)
print("Identity:", top_hsp.identities, "/", top_hsp.align_length)
print("Percent identity:", (top_hsp.identities / top_hsp.align_length) * 100)

                #Identify Conserved Regions

print("Query sequence:")
print(top_hsp.query)

print("Match:")
print(top_hsp.match)

print("Subject sequence:")
print(top_hsp.sbjct)

                #Step 5: Functional Annotation

from Bio import SeqIO
record = SeqIO.read("NS1.fasta", "fasta")
print(record.id)
print(record.description)
print(len(record.seq))
print(record.seq)

                #Predict: Function, Biological role, Organism relevance

#1)Chemically, the NS1 protein does not act as an enzyme and does not catalyze biochemical reactions. 
# Instead, it binds viral RNA and interacts with host proteins, which allows it to regulate host cellular processes.

#2)Biologically, NS1 plays an important role in interfering with the host immune system. 
# It blocks interferon signaling pathways, which normally help the host detect and eliminate viral infections.

#3)For the virus, NS1 is essential because it helps the virus escape immune detection by the host. 
# This immune evasion allows the virus to replicate efficiently and establish infection. 

#Step 6: Biological Interpretation (Research Outcome)

#What does this sequence likely do?

#This sequence likely encodes an immune-modulating protein that helps the virus interfere with the host’s immune response. 
# Specifically, it suppresses interferon signaling, allowing the virus to evade immune detection and replicate efficiently.

#Why do they think so?

#Researchers think this because the sequence shows similarity to known NS1 proteins in BLAST analysis 
# and contains conserved regions associated with RNA binding and host-protein interaction. 
# These features are characteristic of proteins involved in immune suppression rather than enzymatic activity.

#What evidence supports their claim?

#Evidence comes from experimental studies showing that viruses lacking this protein fail to inhibit interferon responses 
# and are rapidly cleared by the host immune system. 
# Structural and functional studies further demonstrate that NS1 directly interacts with RNA and host immune factors,
# confirming its role in immune evasion.