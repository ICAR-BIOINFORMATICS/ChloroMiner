'''
This script takes raw short-reads file (.fastq) and a reference chloroplast genome file (.fasta) as input.
Script's output is the aligned and merged chlorplast genome (.fasta) file from the input short reads.

Workflow of the script -
1 - The reference genome sequence is stored in dictionary format (for searching the short-reads, during alignment)
2 - Empty align_list is created to store the counts of the nucleotides (of the matched short-reads)
3 - Open the input short-reads file and store the reads-batch in dataframe
4 - Alignment of the short reads based on hash-map
        a - search all the possible locations where reads have possibility of occurance
        b - on each location check the short-read in navie manner with some mismatch-exception
        c - if the short is fit for the particular position, increment the align-list position score of the particular base-pair
5 - repeat step 3-4 until the input file is finished
6 - repeat step 3-5 until all the input files are read - if file has negative reads, then reverse complement reads before processing
7 - merge the align_list - generate the chloroplast genome, select the base-pairs which have the highest scores in the list

Design and Development :
Mr. Samarth Godara
Scientist, ICAR-IASRI, New Delhi, India
Email : samarth.godara@icar.gov.in, Mobile : +918114499630

Dr. Shbana Begam
Scientist, ICAR-NIPB, New Delhi, India
Email : shaba.shb@gmail.com

Date of publishing v1 : 22-06-2023
'''

import sys
#pandas used to manipulate the short-reads in dataframe format
import pandas as pd
#to plot the quality-profiling graphs
import matplotlib.pyplot as plt
#to print the starting and finishing time of the execution
from datetime import datetime

#function to create dictionary for searching the possible occurances of short reads
#against the reference genome. Substrings of size key_len are taken from the reference string and used as key.
#The list of positions of those substrings are stored as their corresponding value in the dictionary.
def create_seq_dict(ref_file_name, overlap_size, key_len):
  #open the reference genome file
  rice_ref = open(ref_file_name,'r')
  #print the first line of the file
  line = rice_ref.readline()
  #string to store the whole reference sequence
  seq = ''
  while True:
    #read the file line by line and append the lines in seq variable
    line = rice_ref.readline()
    #remove extra end-line characters from the line
    line = line.replace("\n","")
    #append the line into sequence
    if len(line):
      seq+=line
    else:
      break

  #remove any next line characters from the sequence  seq = seq.replace("\n","")
  print("Total size of the reference genome : ",len(seq))
  #close the reference genome file
  rice_ref.close()
  #convert all characters into UPPER case - atgc => ATGC
  seq=seq.upper()
  #create dictionary for storing the locations of k-mers
  print("Creating hash map of the reference sequence...")
  seq_hash = {}
  #add some part of the end of the reference genome sequence to its start
  #this is done because the chloroplast genome is circular in shape
  #this helps in alignment
  seq = seq[(len(seq)-overlap_size):]+seq
  #creating the key and values for the dictionary
  for i in range(len(seq)-key_len):
    #size of key is defined by key_len
    key = seq[i:(i+key_len)]
    #if the key is already present, then append the presently found location in the list
    if (seq_hash.get(key) != None):
      x = seq_hash[key]
      x.append(i)
      seq_hash[key] = x
    #otherwise if the key is being added in the dictionary for the first time
    else :
      seq_hash[key] = [i]

  print("Size of the sequence dictionary : ",len(seq_hash))
  #display the key which has the maximum locations associated to it
  max_len=0
  max_key = ''

  #for all the keys
  for i in seq_hash.keys():
    #get the list of the particular key
    x = seq_hash[i]
    #check the size of the list
    if len(seq_hash[i])>max_len:
      max_len = len(seq_hash[i])
      max_key = i
  #print the key (k-mer) with the largest number of possible occurances
  print("Biggest key value : ",max_key)
  print("Largest list : ",max_len)
  return (seq,seq_hash)

def create_align_list(seq,overlap_size):
  #creating lists for storing the number of overlaps per base pair against the reference genome
  #these lists are used for merging the aligned reads
  a_list = [0 for i in range(len(seq)+overlap_size)]
  c_list = [0 for i in range(len(seq)+overlap_size)]
  g_list = [0 for i in range(len(seq)+overlap_size)]
  t_list = [0 for i in range(len(seq)+overlap_size)]
  return [a_list, c_list, g_list, t_list]

#open the input files (short read files)
#later this fuction can be used to check the name of file (only .fastq files should be entered)
def open_input_files(p_read_path):
  return open(p_read_path,'r')

#function to check if the line is sequence (ATCG) or not
def isNotSeq(line, comp):
  #for each character of the line, check if it is present in the dictionary comp
  for i in line:
    #if the dictionary doesnt have the character, then the line is not a sequence
    if comp.get(i)==None:
      return True
  #if the whole line is finished and still all the characters are from comp, then the line is a sequence
  return False

#function to return a dataframe including 2 columns - read, score
#the reads and scores are read from the input file, reading limit is defined (maximum batch size) to keep
#the program from crashing the ram, if there is a large number of reads in the file
#rev_comp is boolean used to reverse and compliment the reads before storing them in the dataframe (used in case of negative reads)
def create_dataframe(reads_file, read_limit, rev_comp):
  #there are reads left in the file to be read
  status=True
  i=0
  #dictionary used to compliment the negative reads
  comp = {'A': 'T','T': 'A','G': 'C','C': 'G','N': 'N', '\n':'\n','a': 't','t': 'a','g': 'c','c': 'g','n': 'n'}
  #creating a new dataframe to store the reads in
  dataset = pd.DataFrame(columns = ['read', 'score'])
  #empty lists for the dataframe
  reads = []
  scores = []

  while True:
    #read the lines from each of the two files, in a group of 4
    fileline = reads_file.readline()
    #if the current line is not a sequence then continue the loop and read the next line
    if  isNotSeq(fileline,comp):
      continue
    #if the line is a sequence, then convert it into UPPER case and skip one line after it
    p_read = fileline.upper()
    p_sign = reads_file.readline()
    #after skipping, the next line corresponds to the scores of the sequence read
    p_score = reads_file.readline()

    #if the file contains negative reads, then reverse compliment the reads
    if rev_comp:
      #reverse the negative read and score string
      p_read = p_read[::-1]
      p_score = p_score[::-1]
      p_read = list(p_read)

      #compliment the negative read nucleotides
      for j in range(len(p_read)):
        p_read[j] = comp[p_read[j]]
      p_read= ''.join(p_read)

    #remove any extra spaces from the lines
    if len(p_read)>0:
      p_read = p_read.replace(" ","")
      p_read = p_read.replace("\n","")
      p_read = p_read.replace("\t","")

      p_score = p_score.replace(" ","")
      p_score = p_score.replace("\n","")
      p_score = p_score.replace("\t","")

      #store the positive and negative reads and their scores in one go
      reads.append(p_read)
      scores.append(p_score)
    #if there is no more lines in the file to read
    else:
      #change the status
      status=False
      #close the file
      reads_file.close()
      print("File_closed...")
      break
    #print message after each 1M sequence read
    i+=1
    if i%1000000==0 :
      print("Reads stored : ",i)
    if i>=read_limit:
      break

  #store the lists into dataframe format
  dataset['read']=reads
  dataset['score']=scores

  return (dataset,status)

#converts the first 'mark_len' nucleotides of the short read into a number of base 5
#numb is the dictionary used to map each nucleotide to a number
def create_ada_mark(read, numb, mark_len):
  #extract the begining nucleotide sequence
  mark = read[:mark_len]
  j=0
  number=0
  #begin from the last character
  for i in range(mark_len-1,-1,-1):
    #convert each char to number, then multiply by the 5^position-1
    number += numb[mark[i]]*(5**j)
    j+=1
  return number

#function to detect adapter in the dataset with reads and their scores in it
def adapter_detection(dataset, mark_len, ada_th_per):
  print("Detecting Adapter...")

  #dictionary used to make adapter-markers
  numb = {'A': 0,'C': 1,'G': 2,'T': 3, 'N':4, 'a': 0,'c': 1,'g': 2,'t': 3, 'n':4}

  #convert the begining section of the read into a number of base 5
  dataset['ada_mark'] = dataset.apply(lambda x: create_ada_mark(x['read'], numb, mark_len), axis=1)

  #check the count of each number present in the dataset
  mark_count = dataset['ada_mark'].value_counts()
  #print(mark_count[:10])

  #if the count of the numbers in the dataset exceed the threshold then adapter is present in the dataset
  if (mark_count[0]/dataset.shape[0])*100 > ada_th_per :
    print("Adapter is Present")
    adapter_stat=True
  else :
    print("Adapter is Absent")
    adapter_stat=False

  return (dataset,adapter_stat)

#function to create quality-profiling lists, these lists are used to save the min, max and avg
#read-scores corresponding to the base pari's position
def create_q_prof(read_score_list_len):
  print("Creating score lists...")
  #creating empty lists to store the min, max and avg values of reads
  dum_list1=[500 for i in range(read_score_list_len)]
  dum_list2=[0 for i in range(read_score_list_len)]
  dum_list3=[0 for i in range(read_score_list_len)]
  #            min       max       average
  return [dum_list1, dum_list2, dum_list3]

#function to take a record as input, and store the min, max and avg quality score of the reads
def get_score_stat(q_prof, score):
  #for all nucleotides of the reads
  j=0
  for i in score:
    #convert the quality score in integer form
    s=ord(i)
    #if the score is the minimum at that particular position, then store it
    if s < q_prof[0][j]:
      q_prof[0][j]=s
    #if the score is the maximum at that particular position, then store it
    if s > q_prof[1][j]:
      q_prof[1][j]=s
    #add up all the scores at that particular position, for calculating average
    q_prof[2][j]+=s
    j+=1

#function to trim the reads based on their quality score
def trim_reads(th_score, read, score):
  final_read=''

  #skip the initial lower-quality reads
  i=0
  while i < len(read):
    if (ord(score[i]))>th_score:
      break
    i+=1

  #capture the central high-quality reads
  while i < len(read):
    if (ord(score[i])<th_score):
      break
    final_read+=read[i]
    i+=1

  #return the trimmed read
  return final_read

def quality_profiling(dataset, qp, min_read_size, read_score_list_len, overlap_size):
  print("Trimming reads based on quality score...")
  q_prof=create_q_prof(read_score_list_len)
  #apply the get_score_stat function to all the records of the dataframe
  #get the information about the min, max and avg score of each base pair position
  dataset.apply(lambda x: get_score_stat(q_prof, x['score']), axis=1)

  #get the number of reads9
  no_read = dataset.shape[0]
  #calculate the average score at each position
  for i in range(len(q_prof[0])):
    q_prof[2][i] = q_prof[2][i]/no_read

  #trim the score lists (till the present scores in the list)
  #get the length of the longest read present in the dataset
  read_len=0
  for i in range(len(q_prof[0])):
    if q_prof[0][i]==500:
      read_len=i
      break

  #calculate the overlap of the reference genome from both ends
  #a little extra than the longest read in the dataset
  if overlap_size < read_len:
    print("Warning : Reference-genome overlap size is smaller than the largest short-read...")

  #trim the score lists
  for i in range(3):
    q_prof[i]=q_prof[i][:read_len]

  #creating list for the graph plots
  x=[i for i in range(1,len(q_prof[0])+1)]

  #plot the graph of qality score profiling of the positive reads
  fig = plt.figure()
  ax = fig.gca()
  ax.plot(x,q_prof[0])
  ax.plot(x,q_prof[1])
  ax.plot(x,q_prof[2])
  plt.xlabel('Base Pair Position')
  plt.ylabel('Quality Score')
  plt.title('Quality Score profiling of the Short Reads')
  plt.legend(['Minimum Quality Score','Maximum Quality Score','Average Quality Score'])
  plt.show()

  #obtain the minimum and maximum scores of the reads
  mini = min(q_prof[0])
  maxi = max(q_prof[1])
  th_score = int (mini+((maxi-mini)*qp))

  print("Minimum Score : ",mini)
  print("Maximum Score : ",maxi)
  print("Threshold Score : ",th_score)

  return th_score


#function to trim and filter the reads of the dataset, based on the threshold-score and minimum-read-size
def quality_trimming_filtering(dataset, min_read_size):
  #apply the trimming to all the records
  data = pd.DataFrame()
  data['t_reads'] = dataset['read']

  return data

'''
#function to trim and filter the reads of the dataset, based on the threshold-score and minimum-read-size
def quality_trimming_filtering(dataset, th_score, min_read_size):
  #apply the trimming to all the records
  print("Trimming the reads...\nThis may take a few minutes...")
  dataset['t_reads'] = dataset.apply(lambda x: trim_reads(th_score, x['read'], x['score']), axis=1)

  print("Filtering the reads based on their length...")
  #get the length of every read after trimming
  dataset['read_l'] = dataset.apply(lambda x: len(x['t_reads']), axis=1)

  #filter the reads based on their size : min_read_size
  dataset = dataset[dataset['read_l'] >= min_read_size]['t_reads']
  #convert the series into dataframe
  dataset = pd.DataFrame(dataset)

  return dataset
'''

#get the possible locations of a particular read
def get_loc(seq_hash, key_len, read):
  #get the first characters of the read, as key for dictionary
  key = read[:key_len]
  #return the value stored in dictionary corresponding to that key
  return seq_hash.get(key)

#check where does the read fit in the alignment, from the pool of the possible locations
def align_reads(new_seq, mm_th, numb, align_list, read, locs):
  #for each location, calculate the number of mismatches present while
  #aligning the reads with its possible location
  for loc in locs:
    mm = 0
    for i in range(len(read)):
      if (read[i]!=new_seq[loc+i] ):
        mm+=1
    #if the number of mismatches is lesser than the threshold, then put the
    #base pair in the alignment-list (increment the count of the particular base pair in the list)
    if mm <= mm_th:
      for i in range(len(read)):
        try:
          align_list[numb[read[i]]][loc+i]+=1
        except:
          pass

#function to align
def h_alignment(dataset, seq, seq_hash, key_len, overlap_size, mm_th, align_list):
  #get the possible locations of the particular read from the dictionary
  print("Fetching the locations of possible occurances of the reads...")
  #create a new column with the list of possible locations of the reads
  dataset['locations'] = dataset.apply(lambda x: get_loc(seq_hash, key_len, x['t_reads']), axis=1)

  #remove all the reads which are not present in the dictionary
  f_dataset = dataset.dropna()
  f_data_shape = f_dataset.shape[0]
  print("Total number of reads : ",dataset.shape[0])
  print("Total number of remaining reads : ",f_dataset.shape[0])
  print("Percentage of reads remaining after filtering : ", (f_dataset.shape[0]/dataset.shape[0])*100,"%")

  #add some part of the begining of the reference sequence
  #this is done for facilitating the alignments as the chloroplast genome is circular in shape
  new_seq = seq + seq[overlap_size:overlap_size*2]

  #dictionary used for storing the base pair count in the alignment-list (for merging)
  numb={'A':0,'C':1,'G':2,'T':3}

  #function to align each read, place it in the align_list, at the locations where it can possibly belong
  print("Aligning the filtered short reads...\nThis may take a few minutes...")
  f_dataset.apply(lambda x: align_reads(new_seq, mm_th, numb, align_list, x['t_reads'], x['locations']), axis=1)

  #get the locations of the reference genome which are not covered after alignment
  print("Uncovered locations in the reference genome: ",end='')
  ul_count=0
  for i in range(overlap_size,len(seq)):
    s=0
    s = align_list[0][i]+align_list[1][i]+align_list[2][i]+align_list[3][i]
    #after adding up all the scores, if the score is still 0,
    #means the particular base pair havent been covered even once
    if s==0:
      ul_count+=1
  print(ul_count)

#function to merge the align_list - extract the final output of the program
def merge_align_list(align_list, seq, ref_w, min_overlap, overlap_size):
  print("Merging the aligned reads...")
  #dictionary used to grab which nucleotide has the greatest occurance
  n_char = {0:'A', 1:'C', 2:'G', 3:'T'}
  numb={'A':0,'C':1,'G':2,'T':3}
  ch_gen=''

  #the begining and the end part of the sequence is skipped
  #we overlaped those parts due to circular shape of the genome
  for i in range(overlap_size,len(seq)):
    #initial count and nucleotide
    max_count=0
    max_char ='N'
    #check all four lists for the maximum occurance
    for j in range(4):
      #grab the maximum occuring nucleotide
      if align_list[j][i] > max_count:
        max_count=align_list[j][i]
        max_char=n_char[j]
    #if the occurance-count is less than the threshold, then pick the nucleotide from the reference
    if max_count < min_overlap:
      max_char='x'
      #max_char=seq[i]

    #if the selected nucleotide is a mutation, to be sure, we check that the occurance of it is 1.5 times more than the
    #occurance of nucleotide from the reference genome
    #if max_char!=seq[i]:
    if max_char!='x':
      try:
          if max_count < (align_list[numb[seq[i]]][i])*ref_w:
            max_char=seq[i]
      except:
          #max_char=seq[i]
          pass
    #'''
    #append the picked nucleotide in the output sequence
    #if max_char!= 'x':
    ch_gen=ch_gen+max_char
  #count the number of alterations present in the output vs the reference genome
  alt_count=0
  for i in range(overlap_size,len(seq)):
    #if a mismatch is found, increment the variable
    if seq[i]!=ch_gen[i-overlap_size] and ch_gen[i-overlap_size]!='x':
    #try:
    #  if seq[i]!=ch_gen[i-overlap_size]:
      alt_count+=1
    #except:
    #  pass

  print("Mutation site count :",alt_count)

  return (ch_gen, alt_count)

#function to save the output to a .fasta file
def save_op_file(op_filepath, info, ch_gen):
  print("Writing the output to fasta file...")
  #compile the name of the output file
  o_file = open(op_filepath,'w+')

  #write the info in the file
  o_file.write(">"+info)

  #store the generated genome in the file, 70 characters in each line
  for i in range(len(ch_gen)):
    if i%70==0:
      o_file.write('\n')
    if ch_gen[i]!='x':
      o_file.write(ch_gen[i])

  #close the output file
  o_file.close()
  print("Output file saved by the name :"+op_filepath)

#function to implement dataframe-coversion, adapter-detection, quality-profiling, quality-trimming,
#read-filtering, and alignment of the reads. This function doesnt implement the merging of the reads
def C_miner(filepath, read_limit, mark_len, ada_th_per, qp, min_read_size, read_score_list_len, seq, seq_hash, key_len, overlap_size, mm_th, align_list, negative_r=False):
  #open the input file, with short reads
  reads_file = open_input_files(filepath)
  first_tag=True
  while True :
    #create dataset from the short reads present in the file
    (dataset,status) =create_dataframe(reads_file, read_limit, negative_r)
    #if first_tag:
      #detect adapter
      #(dataset,adapter_stat)=adapter_detection(dataset, mark_len, ada_th_per)
      #quality profiling
      #th_score = quality_profiling(dataset, qp, min_read_size, read_score_list_len, overlap_size)
      #first_tag=False
    #trimming and filtering the reads
    #dataset = quality_trimming_filtering(dataset, th_score, min_read_size)
    dataset = quality_trimming_filtering(dataset, min_read_size)
    #performing the h_alignment algorithm for alignment
    h_alignment(dataset, seq, seq_hash, key_len, overlap_size, mm_th, align_list)
    #if the file has ended then break the loop, otherwise fetch the next batch of reads from the file
    if not status:
      break

#function to generate output file name and the information that is to be written as the first line
def generate_op_file_info(file_name_1, mark_len, min_read_size, qp, key_len, overlap_size, mm_th,  alt_count, min_overlap, ref_w, op_filepath):
  out_file_name = file_name_1[:len(file_name_1)-2]+"_"+str(mark_len)+"_"+str(min_read_size)+"_"+str(round(qp,2))+"_"+str(key_len)+"_"+str(overlap_size)+"_"+str(mm_th)+"_"+str(alt_count)+"_"+str(min_overlap)+"_"+str(ref_w)
  op_file = str(op_filepath) + str(out_file_name) + ".fasta"
  info = "ada_mark_size_"+str(mark_len)+"__size_trim_"+str(min_read_size)+"__read_quality_"+str(round(qp,2))+"__hash_key_len_"+str(key_len)+"__ref_gen_overlap_"+str(overlap_size)+"__mismatch_thresh_"+str(mm_th)+"__mutation_sites_"+str(alt_count)+"__min_overlap_"+str(min_overlap)+"__reference_weight_"+str(ref_w)
  return (op_file, info)

#program's parameters

#file_name_1 = "sra_data"
#file_name_1 = "sra_data_1"
#file_name_1 = "SO_7472_Tea_Leaf_Ext1_WGS_R1"
#file_name_2 = "SO_7472_Tea_Leaf_Ext1_WGS_R2"

#file_name_1 = "1_paired"
#file_name_2 = "2_paired"

#read as many short reads from the file at once
read_limit=20000000
#check the first few base pairs of each read for adapter detection
mark_len = 10
#adapter detection minimum occurance percentage in dataset
ada_th_per=10
#calculate the quality threshold - for quality score trimming
qp=0.6
#minimum size of the read to filter
min_read_size = 20
#how many base pairs to consider for creating hash-map key
key_len=20
#overlap the reference cholorplast genome ends
overlap_size=300
#how many mismatches are allowed per read while alignment
mm_th = 5
#minmum number of accepted overlap
min_overlap=1
#single reads (single file input) or paired reads (two files input)?
paired_reads=True
#weightage of the reference genome nucleotide against a mutation
ref_w=20
#initial length of the quality-profile lists (when the short-read length is unknown)
read_score_list_len=500
#filename of the reference chloroplast genome
#ref_file_name = '/content/drive/MyDrive/Old Data/Research Backups/Chloroplast Genome Project/Dataset/new data/OS_chl.fasta'
#p_read_path = '/content/drive/MyDrive/Old Data/Research Backups/Chloroplast Genome Project/Dataset/new data/'+file_name_1+'.fq'
#n_read_path = '/content/drive/MyDrive/Old Data/Research Backups/Chloroplast Genome Project/Dataset/new data/'+file_name_2+'.fq'
#op_filepath = '/content/drive/MyDrive/Old Data/Research Backups/Chloroplast Genome Project/Dataset/new data/'
#out_name = ''

ref_file_name = sys.argv[1]
p_read_path = sys.argv[2]
n_read_path = sys.argv[3]
op_filepath = sys.argv[4]
out_name = sys.argv[5]

#display the starting time
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Starting Time =", current_time)

#create dictionary for alignment
(seq,seq_hash)=create_seq_dict(ref_file_name, overlap_size, key_len)

#create lists for short-read alignment
align_list = create_align_list(seq,overlap_size)

#perform preprocessing and alignment
C_miner(p_read_path, read_limit, mark_len, ada_th_per, qp, min_read_size, read_score_list_len, seq, seq_hash, key_len, overlap_size, mm_th, align_list, negative_r=not (paired_reads))

#perform preprocessing and alignment of the negative reads
if paired_reads:
  C_miner(n_read_path, read_limit, mark_len, ada_th_per, qp, min_read_size, read_score_list_len, seq, seq_hash, key_len, overlap_size, mm_th, align_list, negative_r= paired_reads)

#merge the resultant output (after alignment)
(ch_gen, alt_count) = merge_align_list(align_list, seq, ref_w, min_overlap, overlap_size)

#generate the strings using the parameters used, for name and info of the output file
#op_file, info = generate_op_file_info(file_name_1, mark_len, min_read_size, qp, key_len, overlap_size, mm_th,  alt_count, min_overlap, ref_w, op_filepath)
op_file, info = generate_op_file_info(out_name, mark_len, min_read_size, qp, key_len, overlap_size, mm_th,  alt_count, min_overlap, ref_w, op_filepath)

#save the output to a file
save_op_file(op_file, info, ch_gen)

#print the finishing time
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Finishing Time =", current_time)

