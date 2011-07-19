#!/usr/bin/env python
import sys,os,load_gb_file_pysqlite

if __name__ == "__main__":
	if len(sys.argv) != 4:
		print("usage: ungz_send_to_load_all_gb_files database dir prefix (like pln)")
		sys.exit(0)
	files = os.listdir(sys.argv[2])
	for i in files:
		if i[0:(2+len(sys.argv[3]))] == "gb"+sys.argv[3]:
			filen = i
			print (filen)
			if filen[-3:] == ".gz":
				os.system("gunzip -d "+filen)
				filen = filen[:-3]
			load_gb_file_pysqlite.load_seqs_into_sqlite(sys.argv[1],filen)
