#include <iostream>
#include <thread>
using namespace std;

void a(const string &input){
	cout<<input<<endl;
};
void b(const string &input){
	cout<<input<<endl;
}

int main(){
	string input = "1";
	thread funcTest1(a, "hello");
	funcTest1.join();
	return 0;
}
