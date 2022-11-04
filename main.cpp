#include "posegraph.h"
#include "utils.h"

int main(int argc, char** argv)
{
	std::string input_vfile = "P4-v.txt";
	std::string input_efile = "P4-e.txt";
	Posegraph pg = Posegraph();
	pg.readGraph(input_vfile, input_efile);
	pg.optimize();
}