import numpy
import sys


def main():

	if len(sys.argv) >=2:
		infile = sys.argv[1]
	else:
		print 'FAIL!'

	fi = open(infile,'rU')
	outfile = infile[:-4] + '.txt'
	fo = open(outfile,'w')

	for line in fi:
  		cols = line.split('  ')
  		lam = cols[0]
  		trans = cols[1]
		
		lam = float(lam)*1e-10
		trans = float(trans)

		print >> fo,lam,'  ',trans
		#print lam,' ',trans

  
if __name__ == '__main__':
	main()
