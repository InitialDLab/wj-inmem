#include "tool.h"
#include <iostream>
using namespace std;
int main(){
	cout<<"0.95: "<<tool::get_zp(0.95)<<std::endl;
	cout<<"0.96: "<<tool::get_zp(0.96)<<std::endl;
	cout<<"0.97: "<<tool::get_zp(0.97)<<std::endl;
	cout<<"0.98: "<<tool::get_zp(0.98)<<std::endl;
	cout<<"0.99: "<<tool::get_zp(0.99)<<std::endl;
	return 0;
}
