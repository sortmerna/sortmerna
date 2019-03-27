#pragma once
/**
 * FILE: reader.hpp
 * Created: Nov 06, 2017 Mon
 *
 * Reads 'Reads file' and pushes the records to a shared queue for further processing
 */

#include <string>
#include <fstream> // std::ifstream

#include "readsqueue.hpp"
#include "kvdb.hpp"
#include "options.hpp"
#include "gzip.hpp"

 // forward
class Read;

/* 
 * reads Reads file and, generates Read objects
 */
class Reader {
public:
	Reader(std::string id, std::ifstream &ifs, bool is_gzipped, KeyValueDatabase & kvdb);
	~Reader();

	Read nextread(Runopts & opts);
	static bool loadReadByIdx(Runopts & opts, Read & read);
	static bool loadReadById(Runopts & opts, Read & read);

public:
	bool is_done = false; // flags end of reads stream

private:
	std::string id;
	bool is_gzipped;
	KeyValueDatabase & kvdb; // key-value database
	Gzip gzip;
	std::ifstream &ifs; // reads file
	unsigned int read_count; // count of reads
	unsigned int line_count; // count of non-empty lines in the reads file
	int last_count;
	int last_stat;
	bool isFastq; // flags the file is FASTQ
	bool isFasta; // flags the file is FASTA
};

// ~reader.hpp