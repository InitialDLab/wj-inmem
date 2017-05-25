#ifndef _H_SAMPLEJOIN
#define _H_SAMPLEJOIN
#include <iostream>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <functional>
#include <random>
#include <chrono>
#include "timer.h"
#include "base.h"
using namespace std;
//! abstract method for samplejoin

//! 3 table sample join, first table is a range structure
template <class T>
struct result_node{
	size_t size;
	T data;
	result_node(){
		size = 0;
	}
};
//! R is the original tuple type
//! K is the key type
//! S is the store type
//! T is the tuple type we want
template <class R, class K, class S, class T>
result_node<T> equal_sample(const R& tuple, std::function<K(const R &)> & key_func, S &store){
	result_node<T> tmp;
	auto key = key_func(tuple);
	tmp.size = store.get_size(key);
	if(tmp.size == 0)
		return tmp;
	tmp.data = store.get_sample(key);
	return tmp;
}
//! get a tuple from given store
template <class R, class K, class S>
result_node<R> equal_sample(const K &key, S &store){
	result_node<R> tmp;
	tmp.size = store.get_size(key);
	if(tmp.size == 0)
		return tmp;
	tmp.data = store.get_sample(key);
	return tmp;
} 
//! tuple is the tuple from previous table
template <class R, class K, class S>
result_node<R> range_sample_more(const K &key, S &store){
	result_node<R> tmp;
	tmp.size = store.count_more(key);
	if(tmp.size == 0)
		return tmp;
	tmp.data = store.get_sample(store.size() - tmp.size, tmp.size);
	return tmp;
}
template <class R, class K, class S>
result_node<R> range_sample_more(const K &key, S &store, const size_t &num){
	result_node<R> tmp;
	tmp.size = num;
	tmp.data = store.get_sample(store.size() - num, tmp.size);
	return tmp;
}
template <class R, class K, class S>
result_node<R> range_sample_less(const K&key, S &store){
	result_node<R> tmp;
	tmp.size = store.count_less(key);
	if(tmp.size == 0)
		return tmp;
	tmp.data = store.get_sample(0, tmp.size);
	return tmp;
}
template <class R, class K, class S>
result_node<R> range_sample_less(const K&key, S &store, const size_t &num){
	result_node<R> tmp;
	tmp.size = num;
	if(tmp.size == 0)
		return tmp;
	tmp.data = store.get_sample(0, num);
	return tmp;
}
//3 table join without cond
template <class R1, class R2, class R3, class F1, class F2, class F3, class K1, class K2, class K3, class S>
void samplejoin(const std::vector<R1> &t1, \
		const std::vector<R2> &t2, \
		const std::vector<R3> &t3, \
		F2 key_func2, \
		F3 key_func3, \
		std::function<K2(const R1 &)> key_cmp1,\
		std::function<K3(const R2 &)> key_cmp2,\
		S result_func, \
		const size_t &STEP, \
		const size_t &MAX,\
		const double &prob){
	auto seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::mt19937 gen;
	gen.seed(seed);
	std::uniform_int_distribution<uint64_t> dis(0, t1.size()-1);
	double sum_y = 0;
	double sum_y2 = 0;
	hash_tree<R2, std::function<K2(const R2 &)>, K2> table2(key_func2);
	table2.build_index(t2);
	hash_tree<R3, std::function<K3(const R3 &)>, K3> table3(key_func3);
	table3.build_index(t3);
	if(t1.size() == 0)
		return;
	tool::start_report();
	size_t n_rejected = 0;
	timer maxtimer;
	timer steptimer;
	size_t round = 0;
	size_t report_round = 1;
	while(maxtimer.get_elapsed() < MAX){
		round++;
		auto ind1 = dis(gen);
		result_node<R1> tuple1;
		tuple1.size = t1.size();
		tuple1.data = t1[ind1];
		auto tuple2 = equal_sample<R1, K2, decltype(table2),R2>(tuple1.data, key_cmp1, table2); 
		if(tuple2.size == 0){
			n_rejected++;
			continue;
		}
		auto tuple3 = equal_sample<R2, K3, decltype(table3), R3>(tuple2.data, key_cmp2, table3);
		if(tuple3.size == 0){
			n_rejected++;
			continue;
		}
		auto tmp = tuple1.size*tuple2.size*tuple3.size*result_func(tuple3.data);
		sum_y += tmp;
		auto y_hat = sum_y / round;
		sum_y2 += std::pow(tmp, 2);
		if(steptimer.get_elapsed()> STEP){
			auto ci = tool::calc_ci(sum_y2, sum_y, y_hat, round, prob);
			tool::report(STEP*report_round, report_round, n_rejected, 0, double(n_rejected)/ round, y_hat, ci, prob);
			steptimer.restart();
			report_round++;
		}
	};
}
template <class R1, class R2, class R3, class F1, class F2, class F3, class K1, class K2, class K3, class S>
void samplejoin(const std::vector<R1> &t1, \
		const std::vector<R2> &t2, \
		const std::vector<R3> &t3, \
		F1 key_func1, \
		F2 key_func2, \
		F3 key_func3, \
		std::function<K2(const R1 &)> key_cmp1,\
		std::function<K3(const R2 &)> key_cmp2,\
		const K1 &key, \
		std::function<bool(const R2 &)> cond2, \
		std::function<bool(const R3&)> cond3, \
		S result_func, \
		const size_t &STEP, \
		const size_t &MAX,\
		const double &prob){
	double sum_y = 0;
	double sum_y2 = 0;
	hash_tree<R1, std::function<K1(const R1 &)>, K1> table1(key_func1);
	table1.build_index(t1);
	hash_tree<R2, std::function<K2(const R2 &)>, K2> table2(key_func2);
	table2.build_index(t2);
	hash_tree<R3, std::function<K3(const R3 &)>, K3> table3(key_func3);
	table3.build_index(t3);
	auto tuple1= equal_sample<R1, K1, decltype(table1)>(key, table1);
	if(tuple1.size == 0)
		return;
	tool::start_report();
	size_t n_rejected_join = 0;
	size_t n_rejected_cond = 0;
	timer maxtimer;
	timer steptimer;
	size_t round = 0;
	size_t report_round = 1;
	while(maxtimer.get_elapsed() < MAX){
		round++;
		auto tuple1 = equal_sample<R1, K1, decltype(table1)>(key, table1);
		auto tuple2 = equal_sample<R1, K2, decltype(table2),R2>(tuple1.data, key_cmp1, table2); 
		if(tuple2.size == 0){
			n_rejected_join++;
			continue;
		}
		if(!cond2(tuple2.data)){
			n_rejected_cond++;
			continue;
		}
		auto tuple3 = equal_sample<R2, K3, decltype(table3), R3>(tuple2.data, key_cmp2, table3);
		if(tuple3.size == 0){
			n_rejected_join++;
			continue;
		}
		if(!cond3(tuple3.data)){
			n_rejected_cond++;
			continue;
		}
		auto tmp = tuple1.size*tuple2.size*tuple3.size*result_func(tuple3.data);
		sum_y += tmp;
		auto y_hat = sum_y / round;
		sum_y2 += std::pow(tmp, 2);
		if(steptimer.get_elapsed()> STEP){
			auto ci = tool::calc_ci(sum_y2, sum_y, y_hat, round, prob);
			tool::report(STEP*report_round, report_round, n_rejected_join, n_rejected_cond, double(n_rejected_join +
						n_rejected_cond)/ round, y_hat, ci, prob);
			steptimer.restart();
			report_round++;
		}
	};
}


template <class R1, class R2, class R3, class F1, class F2, class F3, class K1, class K2, class K3, class S>
void samplejoin_range_less(const std::vector<R1> &t1, \
		const std::vector<R2> &t2, \
		const std::vector<R3> &t3, \
		F1 key_func1, \
		F2 key_func2, \
		F3 key_func3, \
		std::function<K2(const R1 &)> key_cmp1,\
		std::function<K3(const R1 &)> key_cmp2,\
		const K1 &key, \
		std::function<bool(const R2 &)> cond2, \
		std::function<bool(const R3&)> cond3, \
		std::function<bool(const R1&, const R1 &)> cmp_func,\
		S result_func, \
		const size_t &STEP, \
		const size_t &MAX,\
		const double &prob){
	double sum_y = 0;
	double sum_y2 = 0;
	range_tree<R1, std::function<K1(const R1 &)>, K1> table1(key_func1, cmp_func, t1);
	table1.build_index();
	hash_tree<R2, std::function<K2(const R2 &)>, K2> table2(key_func2);
	table2.build_index(t2);
	hash_tree<R3, std::function<K3(const R3 &)>, K3> table3(key_func3);
	table3.build_index(t3);
	auto tuple1= range_sample_less<R1, K1, decltype(table1)>(key, table1);
	auto less_size = tuple1.size;
	if(tuple1.size == 0)
		return;
	cout<<"order hit:"<<double(less_size) / t1.size()<<endl;
	tool::start_report();
	size_t n_rejected_join = 0;
	size_t n_rejected_cond = 0;
	timer maxtimer;
	timer steptimer;
	size_t round = 0;
	size_t report_round = 1;
	while(maxtimer.get_elapsed() < MAX){
		round++;
		auto y_hat = sum_y / round;
		if(steptimer.get_elapsed()> STEP){
			auto ci = tool::calc_ci(sum_y2, sum_y, y_hat, round, prob);
			tool::report(STEP*report_round, report_round, n_rejected_join, n_rejected_cond, double(n_rejected_join +
						n_rejected_cond)/ round, y_hat, ci, prob);
			steptimer.restart();
			report_round++;
		}
		auto tuple1= range_sample_less<R1, K1, decltype(table1)>(key, table1, less_size);
		auto tuple2 = equal_sample<R1, K2, decltype(table2),R2>(tuple1.data, key_cmp1, table2); 
		if(tuple2.size == 0){
			n_rejected_join++;
			continue;
		}
		if(!cond2(tuple2.data)){
			n_rejected_cond++;
			continue;
		}
		auto tuple3 = equal_sample<R1, K3, decltype(table3), R3>(tuple1.data, key_cmp2, table3);
		if(tuple3.size == 0){
			n_rejected_join++;
			continue;
		}
		if(!cond3(tuple3.data)){
			n_rejected_cond++;
			continue;
		}
		auto tmp = tuple1.size*tuple2.size*tuple3.size*result_func(tuple3.data);
		sum_y += tmp;
		sum_y2 += std::pow(tmp, 2);
	}
}



//! first sample is range query 
template <class R1, class R2, class R3, class F1, class F2, class F3, class K1, class K2, class K3, class S>
void samplejoin_range_more(const std::vector<R1> &t1, \
		const std::vector<R2> &t2, \
		const std::vector<R3> &t3, \
		F1 key_func1, \
		F2 key_func2, \
		F3 key_func3, \
		std::function<K2(const R1 &)> key_cmp1,\
		std::function<K3(const R2 &)> key_cmp2,\
		const K1 &key, \
		std::function<bool(const R2 &)> cond2, \
		std::function<bool(const R3&)> cond3, \
		std::function<bool(const R1&, const R1 &)> cmp_func,\
		S result_func, \
		const size_t &STEP, \
		const size_t &MAX,\
		const double &prob){
	double sum_y = 0;
	double sum_y2 = 0;
	range_tree<R1, std::function<K1(const R1 &)>, K1> table1(key_func1, cmp_func, t1);
	table1.build_index();
	hash_tree<R2, std::function<K2(const R2 &)>, K2> table2(key_func2);
	table2.build_index(t2);
	hash_tree<R3, std::function<K3(const R3 &)>, K3> table3(key_func3);
	table3.build_index(t3);
	auto tuple1= range_sample_more<R1, K1, decltype(table1)>(key, table1);
	auto more_size = tuple1.size;
	if(tuple1.size == 0)
		return;
	cout<<"lineitem hit: "<<double(more_size) / t1.size()<<endl;
	tool::start_report();
	size_t n_rejected_join = 0;
	size_t n_rejected_cond = 0;
	timer maxtimer;
	timer steptimer;
	size_t round = 0;
	size_t report_round = 1;
	while(maxtimer.get_elapsed() < MAX){
		round++;
		auto y_hat = sum_y / round;
		if(steptimer.get_elapsed()> STEP){
			auto ci = tool::calc_ci(sum_y2, sum_y, y_hat, round, prob);
			tool::report(STEP*report_round, report_round, n_rejected_join, n_rejected_cond, double(n_rejected_join +
						n_rejected_cond)/ round, y_hat, ci, prob);
			steptimer.restart();
			report_round++;
		}
		auto tuple1= range_sample_more<R1, K1, decltype(table1)>(key, table1, more_size);
		auto tuple2 = equal_sample<R1, K2, decltype(table2),R2>(tuple1.data, key_cmp1, table2); 
		if(tuple2.size == 0){
			n_rejected_join++;
			continue;
		}
		if(!cond2(tuple2.data)){
			n_rejected_cond++;
			continue;
		}
		auto tuple3 = equal_sample<R2, K3, decltype(table3), R3>(tuple2.data, key_cmp2, table3);
		if(tuple3.size == 0){
			n_rejected_join++;
			continue;
		}
		if(!cond3(tuple3.data)){
			n_rejected_cond++;
			continue;
		}
		auto tmp = tuple1.size*tuple2.size*tuple3.size*result_func(tuple1.data);
		sum_y += tmp;
		sum_y2 += std::pow(tmp, 2);
	}
}
#endif
