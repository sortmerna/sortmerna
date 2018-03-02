#pragma once
/* 
 * FILE: cmd.hpp
 * Created: Jan 17, 2018 Wed
 *
 * Interactive command-line session:
 *    - query records in Key-Value database
 */

// forward
struct Runopts;

class CmdSession
{
public:
	CmdSession(){}
	void run(Runopts & opts);
};