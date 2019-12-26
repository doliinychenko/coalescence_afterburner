#include <getopt.h>

#include "smash/stringfunctions.h"

#include <iostream>
#include <string>
#include <vector>

namespace {
void usage(const int rc, const std::string &progname) {
  std::printf("\nUsage: %s [option]\n\n", progname.c_str());
  std::printf(
      "  -h, --help              usage information\n\n"
      "  -i, --inputfiles        <list of particle files>\n"
      "                          should be in SMASH extended binary format\n"
      "  -o, --outputfile        output file name, where the nuclei\n"
      "                          coordinates, momenta, and pdg ids\n"
      "                          will be printed out\n"
      "                          (default: ./nuclei.bin)\n\n");
  std::exit(rc);
}
};  // unnamed namespace

int main(int argc, char **argv) {
  std::string output_file;
  std::vector<std::string> input_files;

  constexpr option longopts[] = {
      {"help", no_argument, 0, 'h'},
      {"inputfiles", required_argument, 0, 'i'},
      {"outputfile", required_argument, 0, 'o'},
      {nullptr, 0, 0, 0}};

  const std::string full_progname = std::string(argv[0]);
  const int i1 = full_progname.find_last_of("\\/") + 1,
            i2 = full_progname.size();
  const std::string progname = full_progname.substr(i1, i2);
  int opt = 0;
  while ((opt = getopt_long(argc, argv, "hi:o:",
          longopts, nullptr)) != -1) {
    switch (opt) {
      case 'h':
        usage(EXIT_SUCCESS, progname);
        break;
      case 'i':
        {
          // A little hack from
          // stackoverflow.com/questions/3939157/c-getopt-multiple-value
          optind--;
          for( ;optind < argc && *argv[optind] != '-'; optind++){
            std::string s(argv[optind]);
            input_files.push_back(s);
          }
         break;
        }
      case 'o':
        output_file = optarg;
        break;
     default:
        usage(EXIT_FAILURE, progname);
    }
  }

  // Abort if there are unhandled arguments left.
  if (optind < argc) {
    std::cout << argv[0] << ": invalid argument -- '" << argv[optind] << "'\n";
    usage(EXIT_FAILURE, progname);
  }


  for (const std::string &s : input_files) {
    std::cout << s << std::endl;
  }
  std::cout << "Output file: " << output_file << std::endl;
}
