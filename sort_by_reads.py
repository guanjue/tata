import numpy as np
def sort_by_reads(tabular_file,output_file):
	data0=open(tabular_file,'r')
	result=open(output_file,'w')
	data01_name=[]
	data01_reads=[]
	result.write(data0.readline())
	for records in data0:
		tmp=[x.strip() for x in records.split('\t')]
		data01_name.append(tmp[0])
		data01_reads.append(tmp[2:])
	### convert to np array
	data01_name=np.array(data01_name)
	data01_reads=np.array(data01_reads,dtype=float)

	### 
	data01_reads_rowsum=np.sum(data01_reads, axis=1)
	sort_id=np.argsort(data01_reads_rowsum)[::-1]

	### sort array
	data01_name_sort=data01_name[sort_id]
	data01_reads_sort=data01_reads[sort_id,]

	for name,reads in zip(data01_name_sort,data01_reads_sort):
		result.write(name+'\t'+name+'\t')
		#print(reads)
		for read in reads:
			result.write(str(read)+'\t')
		result.write('\n')

	result.close()
	data0.close()


############################################################################
import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hf:o:")
	except getopt.GetoptError:
		print 'gff3_expand_sort.py -s <image d2> -t <image d3> -f <image d4> -i <filter1 out> -e <filter2 out> -u <full connect1 out> -l <filter1 d1> -o <filter1 d2> -r <filter2 d1> -v <filter2 d1> -a iter_num -b batch_size -p train_pos_num -n train_neg_num -x max_pooling_1 -y max_pooling_2 -k visual_k-mer_mutation -m visual_Second_k-mer_mutation -d training_speed -j seq_num -q shape_num -w file_list_all -z file_list_random'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'gff3_expand_sort.py -s <image d2> -t <image d3> -f <image d4> -i <filter1 out> -e <filter2 out> -u <full connect1 out> -l <filter1 d1> -o <filter1 d2> -r <filter2 d1> -v <filter2 d1> -a iter_num -b batch_size -p train_pos_num -n train_neg_num -x max_pooling_1 -y max_pooling_2 -k visual_k-mer_mutation -m visual_Second_k-mer_mutation -d training_speed -j seq_num -q shape_num -w file_list_all -z file_list_random'
			sys.exit()
		elif opt=="-f":
			tabular_file=str(arg.strip())
		elif opt=="-o":
			output_file=str(arg.strip())
	sort_by_reads(tabular_file,output_file)

if __name__=="__main__":
	main(sys.argv[1:])
