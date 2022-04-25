
import sys
import pysam
import os

dfh1 = open(sys.argv[2])	# mutation.csv file
dfh2 = pysam.AlignmentFile(sys.argv[1])		# alignment.bam file

Pipeline_Dir =
window = int(sys.argv[3])

for line in dfh1:
	line = line.rstrip().split()
	chrom = line[3]
	pos = int(line[4])
	pnext = pos + 1
	rfh = open(f"{chrom}.{pos}.{sys.argv[1]}_target_mate_depth.txt",'w')
	if_true = 1
	for read in dfh2.fetch(chrom, pos, pnext):
		if read.reference_start <= pos:
			if None not in [read.reference_start, read.reference_end]:
				if dfh2.get_reference_name(read.reference_id) == chrom: 
					line2=str(read).split()
					rfh.write(dfh2.get_reference_name(read.reference_id)+'\t'+str((int(read.reference_start)//window)*window)+'\n')
					rfh.write(dfh2.get_reference_name(read.next_reference_id)+'\t'+str((int(read.next_reference_start)//window)*window)+'\n')
					if_true = 0
	rfh.close()

	if if_true:
		os.system(f"rm {chrom}.{pos}.{sys.argv[1]}_target_mate_depth.txt")
	else:
		os.system(f"{Pipeline_Dir}/pair_depth_bedgraph.sh {chrom}.{pos}.{sys.argv[1]}_target_mate_depth.txt {sys.argv[2]}")
		os.system(f"rm {chrom}.{pos}.{sys.argv[1]}_target_mate_depth.txt")


