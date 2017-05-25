#include <iostream>
#include <vector>
#include "tool.h"
using namespace std;
int main(int argc, char **argv){
	std::vector<std::vector<size_t> > adj{{1}, {0,2}, {1, 3}, {2}};
	auto x = tool::gen_order(4);
	for(auto &i:x){
		for(auto &j:i){
			std::cout<<j<<" ";
		}
		std::cout<<std::endl;
	}
	std::cout<<"after back trace" <<std::endl;
	auto y = tool::back_trace(x, adj);
	for(auto &i:y){
		for(auto &j:i){
			std::cout<<j<<" ";
		}
		std::cout<<std::endl;
	}
	return 0;
}
