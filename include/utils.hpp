#pragma once
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

void split (const std::string &s, char delim, std::vector<std::string> &elems);

std::vector<std::string> split(const std::string &s, char delim);
