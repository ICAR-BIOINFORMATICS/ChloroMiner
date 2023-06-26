#liraries for executing system command and performing folder-related operations
import os
import sys

arg_dic={}

for i in range(1,len(sys.argv)):
	if sys.argv[i]=='-ref':
		arg_dic['ref']=i
	if sys.argv[i]=='-adap':
		arg_dic['adap']=i
	if sys.argv[i]=='-threads':
		arg_dic['threads']=i

#print(arg_dic)

#exit()

#first argument - adapter library
adap_link = sys.argv[arg_dic['adap']+1]

#second argument - chloroplast reference genome
ref_path = sys.argv[arg_dic['ref']+1]

trim_thread=0
try:
	trim_thread = int(sys.argv[arg_dic['threads']+1])
	print("Thread count:", trim_thread)
except:
	print("Thread count by default")
	pass


#get the path of folder where the script is situated
fold_path = os.getcwd()

print(fold_path)

#get the list of folders present in the current folder
#each folder carries different NGS file (forward and reverse)
fold_list = os.listdir(fold_path)

print(fold_list)

#for each NGS library - execute everything present in the loop
for fold in fold_list:

	#if it is the code file, ignore it
	if (fold.find('.') != -1):
		print("Ignoring file:",fold)

		continue
	
	print("Processing folder:", fold)
	#continue


	#prepare the trimmomatic command
	file_list = os.listdir(fold_path+"/"+fold+"/")


	f = file_list[0][file_list[0].find(".")-1]
	#print("f value:", f)
	if f == '1':
		for_file = file_list[0]
		rev_file = file_list[1]
	else:
		for_file = file_list[1]
		rev_file = file_list[0]


	for_file_path = fold_path+'/'+fold+'/'+for_file
	rev_file_path = fold_path+'/'+fold+'/'+rev_file
	print("Path of forward NGS file:", for_file_path)
	print("Path of reverse NGS file:", rev_file_path)

	for_file_paired_path = fold_path+'/'+fold+'/'+fold+"_1_trimmed.fq.gz"
	for_file_unpaired_path = fold_path+'/'+fold+'/'+"forward_unpaired.fq.gz"

	rev_file_paired_path = fold_path+'/'+fold+'/'+fold+"_2_trimmed.fq.gz"
	rev_file_unpaired_path = fold_path+'/'+fold+'/'+"reverse_unpaired.fq.gz"

	if trim_thread==0:
		trimmo_cmd = "TrimmomaticPE "+for_file_path+" "+rev_file_path+" "+for_file_paired_path+" "+for_file_unpaired_path+" "+rev_file_paired_path+" "+rev_file_unpaired_path+" -phred33 ILLUMINACLIP:"+adap_link+":2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36 SLIDINGWINDOW:4:30 MINLEN:36"
	else:
		trimmo_cmd = "TrimmomaticPE "+for_file_path+" "+rev_file_path+" "+for_file_paired_path+" "+for_file_unpaired_path+" "+rev_file_paired_path+" "+rev_file_unpaired_path+" -phred33 ILLUMINACLIP:"+adap_link+":2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36 SLIDINGWINDOW:4:30 MINLEN:36 -threads "+str(trim_thread) 

	print("\n\nTrimmomatic command: ",trimmo_cmd,"\n\n")

	#execute the trimmomatic command
	os.system(trimmo_cmd)


	#unzip the output of trimmomatic command
	unzip_cmd = "gunzip "+for_file_paired_path
	
	print(unzip_cmd,"\n")

	os.system(unzip_cmd)

	unzip_cmd = "gunzip "+rev_file_paired_path
	
	print(unzip_cmd,"\n")

	os.system(unzip_cmd)

	#prepare the cminer command
	#ref_path = fold_path+'/'+'OS_chl.fasta'
	out_path = fold_path+'/'+fold+'/'
	out_name = fold+'_chloro'

	c_miner_cmd = 'python3 '+fold_path+'/c_miner.py '+ref_path+' '+for_file_paired_path[:-3]+' '+rev_file_paired_path[:-3]+' '+out_path+' '+out_name

	print(c_miner_cmd,"\n")

	#execute the cminer
	os.system(c_miner_cmd)

	#break
