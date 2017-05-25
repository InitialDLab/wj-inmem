#ifndef _H_BASE
#define _H_BASE
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <functional>
#include <random>
#include <chrono>
const double root_2 = 1.4142135623;
using namespace std;
struct date_t {
	uint16_t year;
	uint16_t month;
	uint16_t day;
};
const date_t INVALID_DATE=date_t{9999,99,99};
inline uint32_t datetoi(const date_t &d) {
	return d.year * 500 + (d.month-1) * 40 + d.day;
}
inline uint32_t _C2N(const char &x){return x - 48;};

inline date_t strtodate(const std::string &s) {
	uint16_t max_day[13] = {0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
	if (s.length() != 10 ||
		s[4] != '-' ||
		s[7] != '-') return INVALID_DATE;

	date_t d;
	if (isdigit(s[0]) && isdigit(s[1]) && isdigit(s[2]) && isdigit(s[3])) {
		d.year = _C2N(s[0]) * 1000 + _C2N(s[1]) * 100
			+ _C2N(s[2]) * 10 + _C2N(s[3]);
	} else return INVALID_DATE;

	if (isdigit(s[5]) && isdigit(s[6])) {
		d.month = _C2N(s[5]) * 10 + _C2N(s[6]);
		if (d.month > 12) return INVALID_DATE;
	} else return INVALID_DATE;

	if (isdigit(s[8]) && isdigit(s[9])) {
		d.day = _C2N(s[8]) * 10 + _C2N(s[9]);
		if (d.day > max_day[d.month]) return INVALID_DATE;
	} else return INVALID_DATE;

	return d;
}
inline bool operator<(const date_t &t1, const date_t &t2){
	return datetoi(t1) < datetoi(t2);
}
inline std::ostream& operator<<(std::ostream& os, const date_t& d) {
	return os << d.year << '-' << d.month << '-' << d.day;
}
//! we should use two different structure to maintain range and equi-join selections

//! range query
template <class T>
struct range_node{
	T id; //key type to compare
	size_t num;
	range_node(){};
	range_node(const T&t1){
		id = t1;
	}
	range_node(T t1, size_t &t2){
		id = t1;
		num = t2;
	}
};
template <class T>
bool operator<(const range_node<T> &t1, const range_node<T> &t2){
	return t1.id < t2.id;
}
template <class T, class F, class K>
class range_tree{
	public:
		range_tree(F t, std::function<bool(const T&, const T&)> s, const std::vector<T> &input){
			key_func = t;
			cmp_func = s;
			seed = std::chrono::system_clock::now().time_since_epoch().count();
			gen.seed(seed);
			srand(time(NULL));
			data = input;
		}
		//>=
		size_t count_less(const K &t){
			if(data.empty() || index.empty())
				return -1;
			range_node<K> tmp(t);
			auto iter = std::lower_bound(index.begin(), index.end(), tmp);
			return iter->num;
		}
		//>
		size_t count_more(const K &t){
			if(data.empty() || index.empty())
				return -1;
			range_node<K> tmp(t);
			auto iter = std::upper_bound(index.begin(), index.end(), tmp);
			return data.size() - iter->num;
		};
		size_t size(){
			return data.size();
		}
		void build_index(){
			if(data.empty())
				return;
			std::sort(data.begin(), data.end(), cmp_func);
			size_t count = 0;
			auto p = data[0];
			auto q = p;
			index.push_back(range_node<K>(key_func(p), count));
			for(size_t i = 1; i < data.size();++i){
				count++;
				p = data[i];
				if(cmp_func(p, q) || cmp_func(q, p)){
					q = p;
					index.push_back(range_node<K>(key_func(p), count));
				}
			}
		}
		T get_sample(const size_t &start, const size_t &m){
			std::uniform_int_distribution<uint64_t> dis(0, m-1);
			auto ind = dis(gen);
			return data[start + ind];
		}
		~range_tree(){};
	private:
		unsigned seed;
		std::mt19937 gen;
		std::vector<T> data;
		std::vector<range_node<K> > index;
		F key_func;
		std::function<bool(const T&, const T&)> cmp_func;
};

//! equal query
template <class T, class F, class K>
class hash_tree{
	public:
		hash_tree(F input){
			key_func = input;
			seed = std::chrono::system_clock::now().time_since_epoch().count();
			gen.seed(seed);
		};
		void build_index(const std::vector<T> &input){
			size = 0;
			for(auto &x:input){
				K y = key_func(x);
				if(head.count(y) > 0){
					auto iter = head.find(y);
					auto ind = iter->second;
					data[ind].emplace_back(x);
				}
				else
				{
					size_t ind = data.size();
					std::vector<T> tmp;
					tmp.push_back(x);
					data.emplace_back(tmp);
					head.insert(std::make_pair(y, ind));
				}
				size++;
			}
		}
		size_t get_size(){
			return size;
		}
		size_t get_size(const K& t){
			if(head.count(t) == 0)
				return 0;
			else
			{
				auto iter = head.find(t);
				auto ind = iter->second;
				return data[ind].size();
			}
		}
		T get_sample(const K& t){
			auto iter = head.find(t);
			auto ind = iter->second;
			auto pos = get_rand(data[ind].size());
			return data[ind][pos];
		}
		uint64_t get_rand(const uint64_t m)
		{
			std::uniform_int_distribution<uint64_t> dis(0, m-1);
			return dis(gen);
		}
		~hash_tree(){};
	private:
		unsigned seed;
		std::mt19937 gen;
		std::unordered_map<K, size_t> head;
		std::vector<std::vector<T> >data;
		size_t size;
		F key_func;
};
#endif
