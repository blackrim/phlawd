#!/usr/bin/env python
import sys,os,sqlite3

"""
very simple script to correct the merged names from merge.dmp in the sequence table
*DO NOT USE THIS FOR UPDATING* 
this does not have any tests to make sure that names aren't already in there

assumes that python sqlite library would be around
"""

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print("usage: merge_old_names_in_sequences database mergefile")
		sys.exit(0)
	con = sqlite3.connect(sys.argv[1])
	cur = con.cursor()
	#read in the merge.dmp file so we can get the old_id and new_id
	mergef = open(sys.argv[2],"r")
	idkey = {}#[old] = new
	for i in mergef:
		spls = i.strip().split("|")
		old = spls[0].strip()
		idkey[old] = spls[1].strip()
	mergef.close()
	#begin update
	for i in idkey:
		try:
			cur.execute("update sequence set ncbi_id = "+str(idkey[i])+" where ncbi_id = "+str(i)+";")
		except sqlite3.Error as e:
			print ("An error occurred:", e.args[0])
	try:
		con.commit()
	except sqlite3.Error as e:
		print ("An error occurred:",e.args[0])
	con.close()
