/**
 * Contact: philipp.rescheneder@gmail.com
 */

#include <iostream>
#include <string>

#include <tclap/CmdLine.h>

using std::cerr;
using std::cout;
using std::endl;

class ArgParseOutput : public TCLAP::StdOutput
{
private:

	std::string usageStr;

	std::string versionStr;

public:

	ArgParseOutput(std::string usage, std::string version) {
		usageStr = usage;
		versionStr = version;
	}

	virtual ~ArgParseOutput() {

	}

	virtual void failure(TCLAP::CmdLineInterface& c, TCLAP::ArgException& e) {
		cerr << "Error:" << endl;
		cerr << "         " << e.error() << endl;
		cerr << endl;
		cerr << "Short usage:" << endl;
		cerr << "       sniffles [options] -m <sorted.bam> -v <output.vcf> " << endl;
		cerr << endl;
		cerr << "For complete USAGE and HELP type:" << endl;
		cerr << "    sniffles --help" << endl;
		cerr << endl;
		exit(1);
	}

	virtual void usage(TCLAP::CmdLineInterface& c) {
		cerr << usageStr << std::endl;
	}

	virtual void version(TCLAP::CmdLineInterface& c) {

	}
};
