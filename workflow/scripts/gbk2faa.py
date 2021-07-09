#/usr/bin/python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import sys

infile = sys.argv[1]

gbk_filename = infile
faa_filename = "{}.faa".format(gbk_filename)
input_handle  = open(gbk_filename, "rt")
output_handle = open(faa_filename, "w")

for seq_record in SeqIO.parse(input_handle, "genbank") :
    print("Dealing with GenBank record {}".format(seq_record.id))
    for seq_feature in seq_record.features :
        if seq_feature.type=="CDS" :
            try:
               assert len(seq_feature.qualifiers['translation'])==1
               output_handle.write(">{} from {}\n{}\n".format(
                   seq_feature.qualifiers['locus_tag'][0],
                   seq_record.name,
                   seq_feature.qualifiers['translation'][0]))
            except:
                if seq_feature.location.strand == 1:
                   seq_feature.qualifiers['translation'] = seq_record.seq[seq_feature.location.start:seq_feature.location.end].translate(to_stop=True)
                else:
                   seq_feature.qualifiers['translation'] = seq_record.seq[seq_feature.location.start:seq_feature.location.end].reverse_complement().translate(to_stop=True)
                output_handle.write(">{} from {}\n{}\n".format(seq_feature.qualifiers['locus_tag'][0],seq_record.name,seq_feature.qualifiers['translation']))

output_handle.close()
input_handle.close()
