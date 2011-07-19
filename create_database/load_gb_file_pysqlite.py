#!/usr/bin/env python
from Bio import SeqIO
import sqlite3,sys

def load_seqs_into_sqlite(database,filen):
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

def update_seqs_into_sqlite(database,filen):
	handle = open(filen,"rU")
	con = sqlite3.connect(database)
	curup = con.cursor()
	try:
		for i in SeqIO.parse(handle,"gb"):
			acc = i.id
			iden = i.annotations['gi']
			seq = i.seq.tostring()
			desc = i.description
			ncbi_id = ""
			try:
				curup.execute("select accession_id from sequence where accession_id = ?;",(acc,))
				b = curup.fetchall()
				if len(b) == 0:
					a = i.features[0].qualifiers['db_xref']
					for j in a:
						if 'taxon' in j:
							ncbi_id = j[6:]
					curup.execute("insert into sequence (ncbi_id,accession_id,identifier,description,seq) values (?,?,?,?,?);",(ncbi_id,acc,iden,desc,seq))
					print("adding",acc)
			except:
				continue
	except:
		print("issue with genbank file")
	con.commit()
	con.close()
	handle.close()

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print ("usage: load_gb_files_pysqlite database gb_file")
		sys.exit(0)
	print("loading "+sys.argv[2])
	load_seqs_into_sqlite(sys.argv[1],sys.argv[2])
