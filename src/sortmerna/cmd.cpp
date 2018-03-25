/* 
 * FILE: cmd.cpp
 * Created: Jan 17, 2018 Wed
 */

#include <iostream>
#include <sstream>
#include <string>

#include "cmd.hpp"
#include "options.hpp"
#include "kvdb.hpp"
#include "read.hpp"

void CmdSession::run(Runopts & opts)
{
	std::stringstream ss;
	std::string cmd;
	Read read;
	KeyValueDatabase kvdb(opts.kvdbPath);

	for (;;)
	{
		std::cout << "Enter command: [record number, exit]: ";
		read.clear();
		std::cin >> cmd;
		if ("exit" == cmd) break;
		read.init(opts, kvdb, std::stoi(cmd));
		ss << read.matchesToJson() << std::endl;
		std::cout << ss.str(); ss.str("");
	}
} // ~CmdSession::run