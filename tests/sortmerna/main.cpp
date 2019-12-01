/** 
 * FILE: main.cpp
 * Created: Apr 29, 2019 Mon
 */
#include "reader.hpp"

#include <string>
#include <vector>
#include <iomanip> // setprecision

// forward
void kvdb_clear();

/**
 * Case 1
 * Iterate and count all records in provided reads files
 *
 * test.exe 1 C:/a01_projects/clarity_genomics/data/reads/SRR1635864_1.fastq.gz \
 *            C:/a01_projects/clarity_genomics/data/reads/SRR1635864_2.fastq.gz  // 2 X 24,654,594  2 x ~150sec (Win10)
 *
 * @param filev array of full file names
 */
void reader_nextread(std::vector<std::string> &filev)
{
	std::chrono::duration<double> elapsed;

	for (auto rfile : filev)
	{
		auto starts = std::chrono::high_resolution_clock::now();

		std::ifstream ifs(rfile, std::ios_base::in | std::ios_base::binary);
		std::string seq;
		Reader reader("0", true);
		std::cout << "Reads file: " << rfile << std::endl;
		size_t count = 0;
		for (bool hasrec = true; hasrec;)
		{
			hasrec = reader.nextread(ifs, rfile, seq);
			if (hasrec) ++count;
			if (count % 1000000 == 0)
				std::cout << "Number of records processed: " << count << std::endl;
		}
		if (ifs.is_open())
			ifs.close();
		std::cout << "Number of records in file: " << count << std::endl;

		elapsed = std::chrono::high_resolution_clock::now() - starts;
		std::cout << "Time elapsed: [" << std::setprecision(2) << std::fixed << elapsed.count() << "] sec" << std::endl <<std::endl;
	}
} // ~reader_nextread

int main(int argc, char** argv)
{
	std::cout << STAMP << "Running with " << argc << " options" << std::endl;
	//Runopts opts(argc, argv, false);
	if (argc > 2)
	{
		std::cout << "argv[0]: " << argv[0] << std::endl;
		std::cout << "Case: " << argv[1] << std::endl;
		std::string scase = std::string(argv[1]);
		switch (std::stoi(scase))
		{
		case 0:
			kvdb_clear(); 
			break;
		case 1:
			if (argc < 3)
				std::cerr << "Case 1 takes at least 3 args: Test case (1) and one or more full file paths" << std::endl;
			else
			{
				std::vector<std::string> filev;
				for (int i = 2; i < argc; ++i)
					filev.push_back(std::string(argv[i]));
				reader_nextread(filev);
			}
			break;
		default:
			std::cout << "Unknown arg: " << scase << std::endl;
		}
	}
	else
	{
		std::cerr << "Expecting at least one argument: test case e.g. 0 | 1 | 2 etc." << std::endl;
	}

	return 0;
} // ~main