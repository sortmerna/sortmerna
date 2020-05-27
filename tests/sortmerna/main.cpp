/** 
 * FILE: main.cpp
 * Created: Apr 29, 2019 Mon
 */
#include <string>
#include <vector>
#include <iomanip> // setprecision

#include "reader.hpp"
#include "ThreadPool.hpp"
#include "common.hpp"

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
#if 0
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
#endif
} // ~reader_nextread

/* 
 * Memory leak test - push/pop reads on queue and measure the memory use
 * ----------------------------------------------------------------------
 * Reads              with     Memory          Time, sec
 * number of          Read   Start  End    total        per million reads
 * ----------------------------------------------------------------------
 *  2.5M              0      1924   2320     9.4851     3.8     SRR1635864_1_5M.fastq.gz + SRR1635864_2_5M.fastq.gz
 *  2.5M              1		 1892   2284    10.636      4.25
 * 10M                0      1884   2224    38.3424     3.8     SRR1635864_1_20M.fastq.gz + SRR1635864_2_20M.fastq.gz
 * 10M                1      1892   2276    43.1497     4.31
 * 49.3M (49309188)   0      1896   2216   191.04       3.87    SRR1635864_1.fastq.gz + SRR1635864_2.fastq.gz
 * 49.3M (49309188)   1      1888   2112   233.262      4.73
 */
void test_2(std::vector<std::string>& readfiles, size_t num_reads, bool with_read=false)
{
	size_t queue_size_max = 1000;
	ThreadPool tpool(2);
	ReadsQueue read_queue("queue_1", queue_size_max, num_reads);

	{
		std::stringstream ss;
		ss << STAMP << "Running simple memory leak test. 1 Push thread, 1 Pop thread. Memory KB: " << (get_memory() >> 10) << std::endl;
		std::cout << ss.str();
	}

	auto start = std::chrono::high_resolution_clock::now();

	// add Reader job
	tpool.addJob(Reader(read_queue, readfiles, true));
	tpool.addJob([&]() {
		std::string readstr;
		unsigned count = 0;
		for (;read_queue.pop(readstr); ++count) {
			if (with_read)
				Read read(readstr);
		}

		{
			std::stringstream ss;
			ss << "[test_2:Pop thread] Thread ID: " << std::this_thread::get_id() << " popped " << count << " reads" << std::endl;
			std::cout << ss.str();
		}
	});
	tpool.waitAll();

	std::chrono::duration<double> diff = std::chrono::high_resolution_clock::now() - start;

	{
		std::stringstream ss;
		ss << STAMP << "Time elapsed: " << diff.count() << " sec. Memory KB: " << (get_memory() >> 10) << std::endl;
		std::cout << ss.str();
	}
} // ~test_2

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
		case 2:
			if (argc < 4)
				std::cerr << "Case 2 takes at least 4 args: Test case (1), one or two full file paths, "
				"total number of reads, do(1)/don't(0) create Read objects" << std::endl;
			else {
				std::vector<std::string> readfiles;
				size_t num_reads = 0;
				bool with_read = false;
				for (auto i = 2; i < argc; ++i) {
					if (i == 4)
						num_reads = std::stoi(argv[i]);
					else if (i == 5)
						with_read = std::stoi(argv[i]) > 0 ? true : false;
					else
						readfiles.push_back(std::string(argv[i]));
				}
				test_2(readfiles, num_reads, with_read);
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