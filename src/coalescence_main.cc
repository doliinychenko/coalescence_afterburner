#include <getopt.h>

#include "coalescence/coalescence.h"

#include <iostream>
#include <string>
#include <vector>

namespace {
void usage(const int rc, const std::string &progname) {
  std::printf("\nUsage: %s [option]\n\n", progname.c_str());
  std::printf(
      "  -h, --help              usage information\n\n"
      "  -p, --dp                coalescence dp [GeV]\n"
      "  -r, --dr                coalescence dr [fm]\n"
      "  -w, --probabilistic     probabilistic coalescence, 3 exp(-dr2/d2 - dp2 * d2)\n"
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
  using namespace coalescence;
  std::string output_file("nuclei.dat");
  std::vector<std::string> input_files;

  constexpr option longopts[] = {
      {"help", no_argument, 0, 'h'},
      {"dp", required_argument, 0, 'p'},
      {"dr", required_argument, 0, 'r'},
      {"probabilistic", no_argument, 0, 'w'},
      {"inputfiles", required_argument, 0, 'i'},
      {"outputfile", required_argument, 0, 'o'},
      {nullptr, 0, 0, 0}};

  const std::string full_progname = std::string(argv[0]);
  const int i1 = full_progname.find_last_of("\\/") + 1,
            i2 = full_progname.size();
  const std::string progname = full_progname.substr(i1, i2);
  int opt = 0;
  double inputdp = 0.44;  // GeV
  double inputdr = 2.0 * M_PI * 0.19732 / inputdp;  // fm
  bool probabilistic = false;

  while ((opt = getopt_long(argc, argv, "hi:o:p:r:w",
          longopts, nullptr)) != -1) {
    switch (opt) {
      case 'h':
        usage(EXIT_SUCCESS, progname);
        break;
      case 'p':
        inputdp = std::stod(optarg);
        break;
      case 'r':
        inputdr = std::stod(optarg);
        break;
      case 'w':
        probabilistic = true;
        break;
      case 'i':
        {
          // A little hack from
          // stackoverflow.com/questions/3939157/c-getopt-multiple-value
          // to treat multiple arguments
          optind--;
          for( ;optind < argc && *argv[optind] != '-'; optind++){
            std::string s(argv[optind]);
            input_files.push_back(s);
          }
         break;
        }
      case 'o':
        if (optarg) { output_file = optarg; }
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

  std::cout << "Input files: ";
  for (const std::string &input_file : input_files) {
    std::cout << input_file << " ";
  }
  std::cout << "\nOutput file: " << output_file << std::endl;
  if (!probabilistic) {
    std::cout << "\n dp = " << inputdp << ", dr =  " << inputdr << std::endl;
  } else {
    std::cout << "Printing out coalescence weights"
              << " according to deuteron Wigner function." << std::endl;
  }
  Coalescence coalescence(output_file, inputdp, inputdr, probabilistic);
  for (const std::string &input_file : input_files) {
    coalescence.make_nuclei(input_file);
  }
}
