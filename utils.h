#ifndef UTILS_h
#define UTILS_h

#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>
#include "readUtils.h"
#include "api/BamAlignment.h"
#include "api/BamReader.h"

double medianOfVector(std::vector<int> intVector, int vecSize);
int align(const std::string &a, const std::string &b, std::string &a_aligned, std::string &b_aligned);
bool verifyCoordinates(std::string refFasta, std::string inputBamFile);
std::vector<std::string> split(const std::string &s, char delim);
std::string getPaddedIntervalList(std::string intervalList, std::string outputPrefix, int padding_length);
//bool readGermlineFile(std::string germlineList);

#endif

