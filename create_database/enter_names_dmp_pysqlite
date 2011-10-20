#!/usr/bin/env python
import sys,os,sqlite3

"""
very simple script to put all the names from names.dmp into a sqlite database
*DO NOT USE THIS FOR UPDATING* 
this does not have any tests to make sure that names aren't already in there

assumes that python sqlite library would be around
"""

if __name__ == "__main__":
	if len(sys.argv) != 4:
		print("usage: enter_names_dmp_sqlite database namesfile nodesfile")
		sys.exit(0)
	con = sqlite3.connect(sys.argv[1])
	cur = con.cursor()
	#read in the nodes.dmp file so we can get the parent_id and rank
	nodesf = open(sys.argv[3],"r")
	rank = {}#listed by the ncbi_id
	pid = {}#listed by the ncbi_id
	for i in nodesf:
		spls = i.strip().split("|")
		tid = spls[0].strip()
		rank[tid] = spls[2].strip()
		pid[tid] = spls[1].strip()
	nodesf.close()
	#begin insertion
	filen = open(sys.argv[2],"r")
	count = 0
	cmdin = ""
	for i in filen:
		if count % 10000 == 0:
			print(count)
		spls = i.strip().split("|")
		for j in spls:
			j = j.strip()
		gin = spls[0].strip()
		nm = spls[1].strip()
		nm_c = spls[3].strip()
		ednm = spls[1].replace('\"',"_").strip()
		ednm = ednm.replace("\'","_")
		ednm = ednm.replace("\\","_")
		ednm = ednm.replace("/","_")
		ednm = ednm.replace("(","_")
		ednm = ednm.replace(")","_")
		ednm = ednm.replace(".","_")
		ednm = ednm.replace("&","_")
		ednm = ednm.replace(",","_")
		ednm = ednm.replace(" ","_")
		#print(cmdin)
		try:
			cur.execute("insert into taxonomy (ncbi_id,name,name_class,node_rank,parent_ncbi_id,edited_name) values (?,?,?,?,?,?);",(gin,nm,nm_c,rank[gin],pid[gin],ednm))
		except sqlite3.Error as e:
			print ("An error occurred:", e.args[0])
		count += 1
	try:
		con.commit()
	except sqlite3.Error as e:
		print ("An error occurred:",e.args[0])
	filen.close()
	con.close()
