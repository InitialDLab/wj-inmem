#include "loaddata.h"

std::vector<raw_region> load_region( ){
	std::vector<raw_region> result;
	std::ifstream fin;
	fin.open("../data/region.tbl", std::ifstream::in);
	if(!fin.good())
		return result;
	std::string line;
	while(getline(fin,line)){
		auto tmp = parse_raw_region(line);
		result.push_back(tmp);
	}
	std::cout<<"load region done..."<<std::endl;
	return result;
}


std::vector<raw_nation> load_nation( ){
	std::vector<raw_nation> result;
	std::ifstream fin;
	fin.open("../data/nation.tbl", std::ifstream::in);
	if(!fin.good())
		return result;
	std::string line;
	while(getline(fin,line)){
		auto tmp = parse_raw_nation(line);
		result.push_back(tmp);
	}
	std::cout<<"load nation done..."<<std::endl;
	return result;
}
std::vector<raw_supplier> load_supplier( ){
	std::vector<raw_supplier> result;
	std::ifstream fin;
	fin.open("../data/supplier.tbl", std::ifstream::in);
	if(!fin.good())
		return result;
	std::string line;
	while(getline(fin,line)){
		auto tmp = parse_raw_supplier(line);
		result.push_back(tmp);
	}
	std::cout<<"load supplier done..."<<std::endl;
	return result;
}

std::vector<raw_partsupp> load_partsupp( ){
	std::vector<raw_partsupp> result;
	std::ifstream fin;
	fin.open("../data/partsupp.tbl", std::ifstream::in);
	if(!fin.good())
		return result;
	std::string line;
	while(getline(fin,line)){
		auto tmp = parse_raw_partsupp(line);
		result.push_back(tmp);
	}
	std::cout<<"load partsupp done..."<<std::endl;
	return result;
}
std::vector<raw_part> load_part( ){
	std::vector<raw_part> result;
	std::ifstream fin;
	fin.open("../data/part.tbl", std::ifstream::in);
	if(!fin.good())
		return result;
	std::string line;
	while(getline(fin,line)){
		auto tmp = parse_raw_part(line);
		result.push_back(tmp);
	}
	std::cout<<"load part done..."<<std::endl;
	return result;
}
std::vector<raw_customer> load_customer( ){
	std::vector<raw_customer> result;
	std::ifstream fin;
	fin.open("../data/customer.tbl", std::ifstream::in);
	if(!fin.good())
		return result;
	std::string line;
	while(getline(fin,line)){
		auto tmp = parse_raw_customer(line);
		result.push_back(tmp);
	}
	std::cout<<"load customer done..."<<std::endl;
	return result;
}
std::vector<raw_orders> load_orders( ){
	std::vector<raw_orders> result;
	std::ifstream fin;
	fin.open("../data/orders.tbl", std::ifstream::in);
	if(!fin.good())
		return result;
	std::string line;
	while(getline(fin,line)){
		auto tmp = parse_raw_orders(line);
		result.push_back(tmp);
	}
	std::cout<<"load orders done..."<<std::endl;
	return result;
}
std::vector<raw_lineitem> load_lineitem( ){
	std::vector<raw_lineitem> result;
	std::ifstream fin;
	fin.open("../data/lineitem.tbl", std::ifstream::in);
	if(!fin.good())
		return result;
	std::string line;
	while(getline(fin,line)){
		auto tmp = parse_raw_lineitem(line);
		result.push_back(tmp);
	}
	std::cout<<"load lineitem done..."<<std::endl;
	return result;
}
