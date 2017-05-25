#ifndef _H_LOADDATA
#define _H_LOADDATA
#include "rawdata.h"
#include <iostream>
#include <cstring>
#include <fstream>
#include "tool.h"
std::vector<raw_region> load_region();
std::vector<raw_nation> load_nation();
std::vector<raw_supplier> load_supplier();
std::vector<raw_partsupp> load_partsupp();
std::vector<raw_customer> load_customer();
std::vector<raw_orders> load_orders();
std::vector<raw_lineitem> load_lineitem();
std::vector<raw_part> load_part();
#endif
