## extractMatrix_bedGraph.py

import re
import os
import argparse
import pandas as pd
import glob



def readandindex(path, pattern, output):
	
	files = glob.glob(path+"/*/*"+pattern, recursive=True)
	print(path+"/*/*"+pattern)
	index_dict = {}; value_list = []
	index_dict = index_dict.fromkeys(files)
	ID = [os.path.splitext(os.path.basename(i))[0] for i in files]
	print(files);
	print(ID)
	### create empty dictionary
	## start reading lines from bedgraph
	for file in files:
		file_dict={}; empty_list=[]
		with open(file) as inputf:
			rline = inputf.readlines()

		print(file)
		for line in rline:
			indexp = ":".join(line.split("\t")[0:3]) ## chr:pos1:pos2
			valuep = line.split("\t")[3][:-1]   ## density
			file_dict[indexp] = valuep
			empty_list.append(indexp)

		value_list.append(file_dict)
		index_dict[file] = set(empty_list)
	uniqueindex = set.intersection(*index_dict.values())
	new_dict = {}; 
	for index1 in uniqueindex:
		new_list = []
		for index2 in value_list:
			new_list.append(int(float(index2[index1])))
		new_dict[index1] = new_list
	df = pd.DataFrame.from_dict(new_dict, orient='index')
	df.columns = ID;
	df.insert(0, "Index", list(new_dict.keys()), True)
	df.to_csv(output, index=False)	



def findOverlap(dict1, dict2):
	z1 = dict1.intersection(dict2)
	return(z1)


def readandindex2(path, pattern, output):
	
	files = glob.glob(path+"/*/*"+pattern, recursive=True)
	print(path+"/*/*"+pattern)
	index_dict = {}; value_list = []
	index_dict = index_dict.fromkeys(files)
	ID = [os.path.splitext(os.path.basename(i))[0] for i in files]
	print(files);
	print(ID)
	### create empty dictionary
	## start reading lines from bedgraph
	for file in files:
		file_dict={}; empty_list=[]
		with open(file) as inputf:
			rline = inputf.readlines()
		print(file)
		for line in rline:
			indexp = ":".join(line.split("\t")[0:3]) ## chr:pos1:pos2
			valuep = line.split("\t")[3][:-1]   ## density
			file_dict[indexp] = valuep
			empty_list.append(indexp)
		value_list.append(file_dict)
		index_dict[file] = set(empty_list)
	fz1 = next(iter(index_dict))
	fz2 = index_dict[fz1]
	for ld1 in index_dict.keys():
		fz2 = findOverlap(fz2, index_dict[ld1])
	uniqueindex = fz2
	new_dict = {}; 
	for index1 in uniqueindex:
		new_list = []
		for index2 in value_list:
			new_list.append(int(float(index2[index1])))
		new_dict[index1] = new_list
	df = pd.DataFrame.from_dict(new_dict, orient='index')
	df.columns = ID;
	df.insert(0, "Index", list(new_dict.keys()), True)
	df.to_csv(output, index=False)	



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='offer path and pattern of bedgraph file and return shared index, i.e., chr:pos1:po2')
    parser.add_argument('--path', metavar='string', type=str, required=True, nargs="+", 
    	help='directory and will search bedgraph in all subdirectories')
    parser.add_argument('--pattern', metavar='string', type=str, required=True, nargs="+",
    	help='bedgraph pattern')
    parser.add_argument('--output', metavar='string', type=str, required=True, nargs="+", 
    	help='output file name')
    args = parser.parse_args()

    readandindex2(args.path[0], args.pattern[0], args.output[0])



    