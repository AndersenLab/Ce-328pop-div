from __future__ import division
import sys

with open(sys.argv[1], 'r') as alignment:
	alignment_dict = {}
	for line in alignment:
		if line.startswith(">"):
			header = line.rstrip("\n").replace(">","")
			alignment_dict[header] = ''
		else:
			alignment_dict[header] += line.rstrip("\n")

for header,seq in alignment_dict.items():
	outlist = []
	for header2,seq2 in alignment_dict.items():
		identical, total = 0, 0
		for i in range(0, len(seq)):
			if seq[i] == "-" or seq2[i] == "-":
				pass
			elif seq[i] == seq2[i]:
				identical += 1
				total += 1
			else:
				total += 1
		outlist.append(str(round((identical/total)*100, 2)))
	print(header + "\t" + " ".join(outlist))
