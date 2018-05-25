#pragma once
/* 
 * FILE: cmd.hpp
 * Created: Jan 17, 2018 Wed
 *
 * Interactive command-line session:
 *    - query records in Key-Value database
 */
#include <string>

// forward
struct Runopts;

enum CMD { EXIT, READ, INDEX };

class CmdSession
{
public:
	CmdSession(){}
	void run(Runopts & opts);
private:
	void cmdRead(Runopts & opts, std::string & cmd);
	void cmdIndex(Runopts & opts, std::string & cmd);
	void cmdTest(Runopts & opts, std::string & cmd);
	void cmd_max_ref_part(Runopts & opts, std::string & cmd); // ref idx=0 part=1
};