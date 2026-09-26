import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

input_file = "ctcf_positive_clean.fasta"
output_file = "ctcf_negative_clean.fasta"

records = list(SeqIO.parse(input_file, "fasta"))

random.seed(42)

negative_records = []

for i, record in enumerate(records):
    sequence = str(record.seq)

    shuffled = list(sequence)
    random.shuffle(shuffled)
    negative_sequence = "".join(shuffled)

    negative_records.append(
        SeqRecord(
            Seq(negative_sequence),
            id=f"negative_{i+1}",
            description=""
        )
    )

SeqIO.write(negative_records, output_file, "fasta")

print("Positive sequences:", len(records))
print("Negative sequences:", len(negative_records))
print("Saved:", output_file)
