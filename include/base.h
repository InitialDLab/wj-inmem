#ifndef _H_BASE
#define _H_BASE
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <functional>
#include <random>
#include <iostream>
#include <chrono>

using namespace std;
struct date_t {
	uint16_t year;
	uint16_t month;
	uint16_t day;
	date_t operator=(const date_t &t){
		if(this == &t)
			return *this;
		year = t.year;
		month = t.month;
		day = t.day;
		return *this;
	}
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
const double root_2 = 1.4142135623;
//! base type
enum raw_tag_type{RAW_EMPTY, RAW_REGION, RAW_NATION, RAW_CUSTOMER, RAW_PARTSUPP, RAW_PART, RAW_LINEITEM, RAW_ORDERS, RAW_SUPPLIER};
enum table_tag_type{EMPTY, REGION, NATION, CUSTOMER, PARTSUPP, PART, LINEITEM, ORDERS, SUPPLIER};
class base_raw{
	public:
		raw_tag_type raw_tag;
		base_raw *pre;
		base_raw(){
			raw_tag = RAW_EMPTY;
			pre = nullptr;
		}
		virtual ~base_raw(){};
};
typedef std::function<uint32_t(base_raw *)>* key_func_type;
struct plan_result_type{
	size_t d;
	base_raw *p_raw;
	size_t len;
	raw_tag_type raw_type;
	std::vector<base_raw *> vec_p;
	plan_result_type(){
		p_raw = nullptr;
		raw_type = RAW_EMPTY;
		len = 0;
		d = 0;
	}
	plan_result_type &operator=(const plan_result_type &t){
		if(this == &t)
			return *this;
		d = t.d;
		p_raw = t.p_raw;
		vec_p = t.vec_p;
		raw_type = t.raw_type;
		len = t.len;
		return *this;
	};
};
struct head_node_type{
	size_t col_id;
	date_t start_date;
	date_t end_date;
	uint32_t key;
	std::string segment;
	unsigned char flag;
};
typedef struct walk_plan_params{
	double sum_y;
	double sum_y2;
	uint32_t n_samples;
	size_t decision_time;
	double decision_value;
	walk_plan_params(){
		sum_y = 0;
		sum_y2 = 0;
		n_samples = 0;
		decision_time = 0;
		decision_value = 0;
	}
}walk_plan_params;
typedef std::function<double(std::vector<plan_result_type> &)> result_func_type;
typedef std::function<bool(base_raw *)>* cond_func_type;
typedef std::function<double(plan_result_type)>* agg_func_type;
//! we should use two different structure to maintain range and equi-join selections

template <class T, class K>
class range_tree{
	public:
		range_tree(std::vector<T> &input, std::function<K(const T&)> s){
			key_func = s;
			seed = std::chrono::system_clock::now().time_since_epoch().count();
			gen.seed(seed);
			srand(time(NULL));
			data.resize(input.size());
			for(size_t i = 0; i < input.size(); ++i)
				data[i] = &input[i];
			cmp_func =[&](const T* x, const T* y)->bool{
				return key_func(*x) <  key_func(*y);
			};
			if(data.empty())
				cout<<"range tree size 0"<<endl;
			else{
				std::sort(data.begin(), data.end(), cmp_func);
			}
		}
		T* get_sample(const K &p, const K &q, size_t &d){
			//auto cmp = [&](const K&x, const T*y)->bool{return key_func(*y) <x;};
			auto iter_p = std::lower_bound(data.begin(), data.end(), p);
			auto iter_q = std::upper_bound(data.begin(), data.end(), q);
			size_t pos_p, pos_q;
			pos_p = iter_p - data.begin();
			pos_q = iter_q - data.begin();
			d = pos_q - pos_p;
			if(d == 0){
				return nullptr;
			}
			std::uniform_int_distribution<uint64_t> dis(0, d-1);
			return data[pos_p + dis(gen)];
		}
		std::pair<base_raw*, size_t> get_all(const K &p, const K&q){
			auto iter_p = std::lower_bound(data.begin(), data.end(), p);
			auto iter_q = std::upper_bound(data.begin(), data.end(), q);
			return std::make_pair(*(*iter_p), iter_q - iter_p);
		}
		size_t get_size(const K &p, const K& q){
			auto iter_p = std::lower_bound(data.begin(), data.end(), p);
			auto iter_q = std::upper_bound(data.begin(), data.end(), q);
			size_t pos_p, pos_q;
			pos_p = iter_p - data.begin();
			pos_q = iter_q - data.begin();
			size_t d = pos_q - pos_p;
			return d;
		}
		size_t size(){
			return data.size();
		}
		~range_tree(){};
	private:
		unsigned seed;
		std::mt19937 gen;
		std::vector<T *> data;
		std::function<K(const T&)> key_func;
		std::function<bool(const T*, const T*)> cmp_func;
};

//! equal query
template <class T, class K>
class hash_tree{
	public:
		hash_tree(std::vector<T> &input, std::function<K(const T&)> s){
			key_func = s;
			seed = std::chrono::system_clock::now().time_since_epoch().count();
			head.reserve(input.size());
			gen.seed(seed);
			count = 0;
			for(auto &x:input){
				++count;
				K y = key_func(x);
				if(head.count(y) > 0){
					auto iter = head.find(y);
					auto ind = iter->second;
					data[ind].emplace_back(&x);
				}
				else
				{
					size_t ind = data.size();
					std::vector<T*> tmp;
					tmp.push_back(&x);
					data.emplace_back(tmp);
					head.insert(std::make_pair(y, ind));
				}
			}
		}
		T* get_sample(const K& t, size_t &d){
			auto iter = head.find(t);
			if(iter == head.end()){
				d = 0;
				return nullptr;
			}
			auto ind = iter->second;
			auto pos = get_rand(data[ind].size());
			d = data[ind].size();
			return data[ind][pos];
		}
	 	std::pair<base_raw*, size_t> get_all(const K&t){
			auto iter = head.find(t);
			if(iter == head.end()){
				return std::make_pair(nullptr, 0);
			}
			auto ind = iter->second;
			base_raw* p = dynamic_cast<base_raw*>(data[ind][0]);
			return std::make_pair(p, data[ind].size());
		}
		size_t get_size(const K&t){
			auto iter = head.find(t);
			if(iter == head.end()){
				return 0;
			}
			auto ind = iter->second;
			return data[ind].size();
		}
		uint64_t get_rand(const uint64_t &m)
		{
			std::uniform_int_distribution<uint64_t> dis(0, m-1);
			return dis(gen);
		}
		size_t size(){
			return count;
		}
		~hash_tree(){};
	private:
		unsigned seed;
		std::mt19937 gen;
		std::unordered_map<K, size_t> head;
		std::vector<std::vector<T*> > data;
		size_t count;
		std::function<K(const T&)> key_func;
};
#endif
