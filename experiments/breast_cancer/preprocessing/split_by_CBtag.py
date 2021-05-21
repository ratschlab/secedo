##### Code has not been tested on unsorted bam files, sort on barcode (CB):
##### samtools sort -t CB unsorted.bam > sorted_tags.bam
###
##### INPUT: .bam file to be sorted and output directory to place split BC
##### OUTPUT: .bam file for each unique barcode, best to make a new directory

### Python 3.6.8
import pysam
import time
import getopt
import sys



def main(argv):
  # file to split on
  unsplit_file = ''
  # where to place output files
  out_dir = ''
  # file with allowed CB tags
  allowedTags_file = ''

  try:
    opts, args = getopt.getopt(argv,"f:o:t:")
  except getopt.GetoptError:
    sys.exit(2)

  for opt, arg in opts:
    if opt == '-f':
      unsplit_file = arg
    elif opt == '-o':
      out_dir = arg
    elif opt == '-t':
      allowedTags_file = arg


  #start_time = time.time()
  # read the allowed tags
  with open(allowedTags_file, 'r') as f:
    allowedTags = f.readlines()
  allowedTags = [tag[:-1] for tag in allowedTags]
  #print(allowedTags)
  

  # variable to hold barcode index
  CB_hold = 'unset'
  itr = 0
  allowed = False
  # read in upsplit file and loop reads by line
  samfile = pysam.AlignmentFile( unsplit_file, "rb")
  for read in samfile.fetch( until_eof=True):
      # barcode itr for current read
      CB_itr = read.get_tag( 'CB')
      # if change in barcode or first line; open new file  
      if( CB_itr!=CB_hold or itr==0):
          # close previous split file, only if not first read in file
          if( itr!=0):
              split_file.close()
          CB_hold = CB_itr
          
          # if in allowed tags, open a new file and set allowed to true and increase the counter
          if CB_itr in allowedTags:
            allowed = True
            split_file = pysam.AlignmentFile( out_dir + "CB_{}.bam".format( CB_itr), "wb", template=samfile)
            itr+=1
          else:
            allowed = False

      # if allowed, write read with same barcode to file
      if allowed:
        split_file.write( read)
  split_file.close()
  samfile.close()

  #print("--- %s seconds ---" % (time.time() - start_time))




   

if __name__ == "__main__":
   main(sys.argv[1:])




