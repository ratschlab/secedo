import sys
import getopt

def varCalling_metrics(ind):
	resFile = "varCalls_chr"+str(ind)
	trueFile = "true_chr"+str(ind)

	rf = open(resFile, 'r')
	tf = open(trueFile, 'r')

	fp = 0
	fn = 0
	tp = 0

	rpos = int(rf.readline().split("\t")[1])
	tpos = int(tf.readline().split("\t")[1])

	while True:
		if rpos == float('inf') and tpos == float('inf'):
			break
		elif rpos<tpos:
			fp += 1
			line = rf.readline()
			if not line:
				rpos = float('inf')
			else:
				rpos = int(line.split("\t")[1])
		elif rpos>tpos:
			fn += 1
			line = tf.readline()
			if not line:
				tpos = float('inf')
			else:
				tpos = int(line.split("\t")[1])
		else:
			tp += 1
			line = rf.readline()
			if not line:
				rpos = float('inf')
			else:
				rpos = int(line.split("\t")[1])
			line = tf.readline()
			if not line:
				tpos = float('inf')
			else:
				tpos = int(line.split("\t")[1])


	rf.close()
	tf.close()
	
	with open("res_all", 'a') as f:
		f.write(','.join([str(ind), str(tp), str(fp), str(fn)])+"\n")


def main(argv):
	# parse the arguments
	try:
		opts, args = getopt.getopt(argv,"i:")
	except getopt.GetoptError:
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-i':
			ind = arg
		
	varCalling_metrics(ind)



if __name__ == "__main__":
    main(sys.argv[1:])