fi = "colo829_chr_8_mpileup"
splitDir = "split_8"
splitFile = splitDir+"/colo829_chr_8_mpileup_0"

f = open(fi, 'r')
sf = open(splitFile, 'w')
n_lines=0
index=0

while True:
	line = f.readline()

	if not line: 
		break

	n_lines += 1
	if n_lines < 20000:
		sf.write(line)
	elif n_lines==20000:
		line = line.rstrip("\n")
		# split by tabs
		splitLine = line.split("\t")
		# second field: position
		last_pos = int(splitLine[1])
		sf.write(line+"\n")
	else:
		line = line.rstrip("\n")
		# split by tabs
		splitLine = line.split("\t")
		# second field: position
		line_pos = int(splitLine[1])
		# gap at least 1000 residues - we can split here
		if line_pos>last_pos+1000:
			sf.close()
			print("File " + str(index) + " finished, it has " + str(n_lines) + " lines.")
			index += 1
			n_lines = 1
			sf = open(splitDir+"/colo829_chr_8_mpileup_"+str(index), 'w')
			sf.write(line+"\n")
		else:
			last_pos = line_pos
			sf.write(line+"\n")

print("File " + str(index) + " finished, it has " + str(n_lines) + " lines.")
print("Finished.")
sf.close()
f.close()

