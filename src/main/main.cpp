#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <omp.h>

#include "StopClock.hpp"
#include <sys/resource.h>
#include <sys/time.h>

#include "parser.hpp"



using namespace seqan3;

// Source: https://github.com/JensUweUlrich/ReadBouncer/blob/master/src/main/main.cpp
double cputime(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

/**
 * Returns the peak (maximum so far) resident set size (physical
 * memory use) measured in bytes, or zero if the value cannot be
 * determined on this OS.
 */
long getPeakRSS(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_maxrss * 1024;
}


int main(int argc, char const **argv)
{

	StopClock megaXUsage;

	if (__cplusplus != 202002L) {
        std::cerr << "[INFO-DEV] MegaX supports only C++20, your version is: " << __cplusplus << std::endl;
        return 1;
    }

	else{

		std::cout << "MegaX passed C++20 version check! your version is: " << __cplusplus << std::endl;
	}

	seqan3::argument_parser MegaXParser("MegaX", argc, argv);
	CmdArguments args { };
	initializeMainArgumentParser(MegaXParser, args);

	initializeArgumentParser(MegaXParser, args);

    try
	{
		MegaXParser.parse();
	}
	catch (seqan3::argument_parser_error const& ext)
	{
		std::cerr << "[PARSER ERROR] " << ext.what() << '\n';
		return -1;
	}

	megaXUsage.start();
	runProgram(args);
	megaXUsage.stop();

	size_t peakSize = getPeakRSS();
	int peakSizeMByte = (int)(peakSize / (1024 * 1024));

	std::filesystem::path bin = std::filesystem::current_path();
	const std::string memoryUsageReport = bin.string() + "/memory.txt";

	std::ofstream memoryLog(memoryUsageReport);
	memoryLog << "*********************** MegaX Complete Usage Report ******" << '\n';
	memoryLog << "* Real time : " << megaXUsage.elapsed() << " sec         " << '\n';
	memoryLog << "* CPU time  : " << cputime() << " sec                    " << '\n';
	memoryLog << "* Peak RSS  : " << peakSizeMByte << " MByte              " << '\n';
	memoryLog << "**********************************************************" << '\n';

	std::cout << '\n';
	std::cout << "*********************** MegaX Complete Usage Report ******" << std::endl;
	std::cout << "* Real time : " << megaXUsage.elapsed() << " sec         " << std::endl;
	std::cout << "* CPU time  : " << cputime() << " sec                    " << std::endl;
	std::cout << "* Peak RSS  : " << peakSizeMByte << " MByte              " << std::endl;
	std::cout << "**********************************************************" << std::endl;
	
	return 0;
}


