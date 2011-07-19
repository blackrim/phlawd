#!/usr/bin/env python

import sys,sqlite3,os

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print ("usage: update_database database_name division")
		print ("\t\tdatabase_name = any database name that will be referred to later")
		print ("\t\tdivision = the division as recognized by NCBI (used for downloading)")
		print ("\t\t\texample: use pln for the plant division")
		sys.exit(0)
	database = sys.argv[1]
	div = sys.argv[2]
	
	if os.path.exists(database) == False:
		print("database file doesn't exists -- run load_database instead")
		sys.exit(0)
	
	con = sqlite3.connect(database)
	
	print("downloading taxonomy")
	os.system("wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz")
	os.system("tar -xzvf taxdump.tar.gz")
	os.system("./update_names_dmp_pysqlite "+database+" names.dmp nodes.dmp")
	os.system("./rebuild_tree_pysqlite "+database)
	
	print("downloading sequences")
	#os.system("wget ftp://ftp.ncbi.nih.gov/genbank/gb"+div+"*.seq.gz")
	os.system("gunzip -d gb"+div+"*.seq.gz")
	print("loading sequences")
	os.system("./ungz_send_to_update_all_gb_files "+database+" . "+div)

	#merge old ids with new ids in sequences
	print("merging old ids with new ids")
	os.system("./merge_old_names_in_sequences "+database+" merged.dmp")

	print("done updating "+database)

	con.close()
