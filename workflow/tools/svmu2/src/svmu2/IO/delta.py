'''
Docstring for IO.delta
read delta file format into models.alignment.Alignment object
'''

from svmu2.models.alignment import Alignment

def parse_delta_file(file_path):
    alignments = []
 
    with open(file_path) as f:
        next(f)
        next(f)

        aln = None
        index = 0

        for line in f:
            line = line.strip()

            if line.startswith(">"):
                if aln:
                    alignments.append(aln)

                parts = line.strip(">").split()
                aln = Alignment(parts[0], parts[1], parts[2], parts[3])
                index = 0

            elif len(line.split()) == 7:
                indel_map = []

                while True:
                    nxt = next(f).strip()
                    val = int(nxt)
                    indel_map.append(val)
                    if val == 0:
                        break

                aln.add_alignment_block(line.split(), indel_map, index)
                index += 1

            else:
                raise ValueError("Potentially malformed delta file")

    if aln:  # Append the last alignment if exists
        alignments.append(aln)

    return alignments
