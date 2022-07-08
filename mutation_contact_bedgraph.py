import sys
import pysam
import os

dfh1 = open("path/to/mutation/mutation.csv")
dfh2 = pysam.AlignmentFile(sys.argv[1])

window = int(sys.argv[2])

sample_ID = sys.argv[1].split(".markdup")[0]

n = 0

for line in dfh1:
        if n != 0:
                line = line.rstrip().split(",")
                line = [x[1:-1] for x in line]
                variant_ID = line[0]
                if variant_ID == "":
                        variant_ID = "N"
                chrom = line[3]
                pos = int(line[4])
                pnext = pos + 1
                variant_classify = line[6]
                variant_type = line[7]
                ref = line[8]
                alt = line[9]

                is_none = 1
                rfh = open(f"{chrom}.{pos}.{sample_ID}.{variant_ID}_target_mate_depth.txt",'w')
                rfh.write(f"#{sample_ID}\t{variant_ID}\t{chrom}\t{pos}\t{ref}\t{alt}\n")
                for read in dfh2.fetch(chrom, pos, pnext):
                        line2 = str(read).split()
                        rfh.write(dfh2.get_reference_name(read.next_reference_id)+'\t'+str((int(read.next_reference_start)//window)*window)+'\n')

                        is_none = 0
                rfh.close()

                if is_none:
                        os.system(f"rm {chrom}.{pos}.{sample_ID}.{variant_ID}_target_mate_depth.txt")
                else:
                        os.system(f"{HNT_Dir}/pair_depth_bedgraph.sh {chrom}.{pos}.{sample_ID}.{variant_ID}_target_mate_depth.txt {sys.argv[2]}")
                        os.system(f"rm {chrom}.{pos}.{sample_ID}.{variant_ID}_target_mate_depth.txt")

        n += 1
os.system(f"cat chr*.bedgraph > all.merged.{sys.argv[1][:-4]}.bedgraph")
os.system(f"sort -k1,1V -k2,2V all.merged.{sys.argv[1][:-4]}.bedgraph > ../Bedgraph/all.merged.{sys.argv[1][:-4]}.sort.bedgraph")
os.system(f"rm all.merged.{sys.argv[1][:-4]}.bedgraph")
os.system(f'rm chr*bedgraph')
