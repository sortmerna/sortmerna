#pragma once
/**
 * FILE: read_control.hpp
 * Created: Mar 12, 2019 Tue
 *
 * @copyright 2016-19 Clarity Genomics BVBA
 */
#include <vector>

#include "options.hpp"
#include "reader.hpp"

// forward
class ReadsQueue;
class KeyValueDatabase;

class ReadControl
{
public:
	ReadControl(Runopts & opts, ReadsQueue & readQueue, KeyValueDatabase & kvdb);
	~ReadControl();

	void operator()() { run(); }
	void run();

private:
	Runopts &opts;
	ReadsQueue &readQueue;
	KeyValueDatabase &kvdb;
	std::vector<Reader> vreader;
};

