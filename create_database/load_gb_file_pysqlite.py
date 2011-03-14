#!/usr/bin/python
from Bio import SeqIO
import sqlite3,sys

def load_seqs_into_seqlite(database,filen):
	handle = open(filen,"rU")
	con = sqlite3.connect(database)
	curup = con.cursor()
	for i in SeqIO.parse(handle,"gb"):
		acc = i.id
		iden = i.annotations['gi']
		seq = i.seq.tostring()
		desc = i.description
		ncbi_id = ""
		try:
			a = i.features[0].qualifiers['db_xref']
			for j in a:
				if 'taxon' in j:
					ncbi_id = j[6:]
			curup.execute("insert into sequence (ncbi_id,accession_id,identifier,description,seq) values (?,?,?,?,?);",(ncbi_id,acc,iden,desc,seq))
		except:
			continue
	con.commit()
	con.close()
	handle.close()

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print "usage: load_gb_files_pysqlite database gb_file"
		sys.exit(0)
	print("loading "+sys.argv[2])
	load_seqs_into_seqlite(sys.argv[1],sys.argv[2])
