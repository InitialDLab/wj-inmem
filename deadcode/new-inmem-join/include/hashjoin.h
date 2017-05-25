#ifndef _H_HASHJOIN
#define _H_HASHJOIN
/*
   internal hash join method
   R1 raw type1
   R2 raw type2
   C1 conditional function pointer table1
   C2 conditional function pointer table2
   keyfunc1 extractkey from table1
   keyfunc2 extract key from table2
   S rawtype of result
*/
#include <iostream>
#include <unordered_map>
#include <vector>
#include <memory>
#include <functional>
using namespace std;
template <class R1, class R2, class S, class F1, class F2, class K1, class K2>
std::vector<S> hash_join(std::vector<R1> &table1, std::vector<R2> &table2, F1 key_func1,
		F2 key_func2){
	//! internal hashjoin with no condition on two vecs
	
	//! R[x] is rawtype of x-th vec
	//! S is the type of result
	//! F[] is the key_extract func
	//! K[] is the key type for each table
	std::vector<S> result;
	if(table1.size() == 0 || table2.size() == 0)
		return result;
	std::vector<R1> *p_table1;
	std::vector<R2> *p_table2;
	std::vector<R1> table1_new;
	std::vector<R2> table2_new;
	p_table1 = &table1;
	p_table2 = &table2;
	if(p_table1->size() <= p_table2->size()){
		std::unordered_multimap<K1, R1> table1_map;
		for(size_t i = 0; i < p_table1->size(); ++i){
			table1_map.insert(std::make_pair(key_func1((*p_table1)[i]), (*p_table1)[i]));
		}
		for(auto &x:table2){
			if(table1_map.count(key_func2(x)) > 0){
				auto range = table1_map.equal_range(key_func2(x));
				auto iter = range.first;
				for(; iter != range.second; ++iter){
					result.push_back(iter->second + x);
				}
			}
		}
		return result;
	}
	else
	{
		std::unordered_multimap<K2, R2> table2_map;
		for(size_t i = 0; i < p_table2->size(); ++i){
			table2_map.insert(std::make_pair(key_func2((*p_table2)[i]), (*p_table2)[i]));
		}
		for(auto &x:table1){
			if(table2_map.count(key_func1(x)) > 0){
				auto range = table2_map.equal_range(key_func1(x));
				auto iter = range.first;
				for(; iter != range.second; ++iter){
					result.push_back(x + iter->second);
				}
			}
		}
		return result;
	}
};
template <class R1, class R2, class S, class C1, class F1, class F2, class K1, class K2>
std::vector<S> hash_join_left(std::vector<R1> &table1, std::vector<R2> &table2, C1 cond1, F1 key_func1,
		F2 key_func2){
	std::vector<S> result;
	if(table1.size() == 0 || table2.size() == 0)
		return result;
	std::vector<R1> *p_table1;
	std::vector<R2> *p_table2;
	std::vector<R1> table1_new;
	for(auto &x:table1){
		if(cond1(x)){
			table1_new.push_back(x);
		}
	}
	p_table1 = &table1_new;
	p_table2 = &table2;
	if(p_table1->size() <= p_table2->size()){
		std::unordered_multimap<K1, R1> table1_map;
		for(size_t i = 0; i < p_table1->size(); ++i){
			table1_map.insert(std::make_pair(key_func1((*p_table1)[i]), (*p_table1)[i]));
		}
		for(auto &x:table2){
			if(table1_map.count(key_func2(x)) > 0){
				auto range = table1_map.equal_range(key_func2(x));
				auto iter = range.first;
				for(; iter != range.second; ++iter){
					result.emplace_back(iter->second + x);
				}
			}
		}
		return result;
	}
	else
	{
		std::unordered_multimap<K2, R2> table2_map;
		for(size_t i = 0; i < p_table2->size(); ++i){
			table2_map.insert(std::make_pair(key_func1((*p_table2)[i]), (*p_table2)[i]));
		}
		for(size_t i = 0; i < p_table1->size(); ++i){
			auto x = (*p_table1)[i];
			if(table2_map.count(key_func1(x)) > 0){
				auto range = table2_map.equal_range(key_func1(x));
				auto iter = range.first;
				for(; iter != range.second; ++iter){
					result.emplace_back(x + iter->second);
				}
			}
		}
		return result;
	}
};
template <class R1, class R2, class S, class C2, class F1, class F2, class K1, class K2>
std::vector<S> hash_join_right(std::vector<R1> &table1, std::vector<R2> &table2, C2 cond2, F1 key_func1,
		F2 key_func2){
	std::vector<S> result;
	if(table1.size() == 0 || table2.size() == 0)
		return result;
	std::vector<R1> *p_table1;
	std::vector<R2> *p_table2;
	std::vector<R2> table2_new;
	for(auto &x:table2){
		if(cond2(x)){
			table2_new.push_back(x);
		}
	}
	p_table2 = &table2_new;
	p_table1 = &table1;
	if(p_table1->size() <= p_table2->size()){
		std::unordered_multimap<K1, R1> table1_map;
		for(size_t i = 0; i < p_table1->size(); ++i){
			table1_map.insert(std::make_pair(key_func1((*p_table1)[i]), (*p_table1)[i]));
		}
		for(size_t i = 0; i < p_table2->size(); ++i){
			auto ele = (*p_table2)[i];
			if(table1_map.count(key_func2(ele)) > 0){
				auto range = table1_map.equal_range(key_func2(ele));
				auto iter = range.first;
				for(; iter != range.second; ++iter){
					result.emplace_back(iter->second + ele);
				}
			}
		}
		return result;
	}
	else
	{
		std::unordered_multimap<K2, R2> table2_map;
		for(size_t i = 0; i < p_table2->size(); ++i){
			table2_map.insert(std::make_pair(key_func2((*p_table2)[i]), (*p_table2)[i]));
		}
		for(size_t i = 0; i < p_table1->size(); ++i){
			auto ele = (*p_table1)[i];
			if(table2_map.count(key_func1(ele)) > 0){
				auto range = table2_map.equal_range(key_func1(ele));
				auto iter = range.first;
				for(; iter != range.second; ++iter){
					result.emplace_back(ele + iter->second );
				}
			}
		}
		return result;
	}
};
template <class R1, class R2, class S, class C1, class C2, class F1, class F2, class K1, class K2>
std::vector<S> hash_join(std::vector<R1> &table1, std::vector<R2> &table2, C1 cond1, C2 cond2, F1 key_func1,
		F2 key_func2){
	std::vector<S> result;
	if(table1.size() == 0 || table2.size() == 0)
		return result;
	std::vector<R1> *p_table1;
	std::vector<R2> *p_table2;
	std::vector<R2> table2_new;
	std::vector<R1> table1_new;
	for(auto &x:table1){
		if(cond1(x)){
			table1_new.push_back(x);
		}
	}
	p_table1 = &table1_new;

	for(auto &x:table2){
		if(cond2(x)){
			table2_new.push_back(x);
		}
	}
	p_table2 = &table2_new;
	if(p_table1->size() <= p_table2->size()){
		std::unordered_multimap<K1, R1> table1_map;
		for(size_t i = 0; i < p_table1->size(); ++i){
			table1_map.insert(std::make_pair(key_func1((*p_table1)[i]), (*p_table1)[i]));
		}
		for(size_t i = 0; i < p_table2->size(); ++i){
			auto ele = (*p_table2)[i];
			if(table1_map.count(key_func2(ele)) > 0){
				auto range = table1_map.equal_range(key_func2(ele));
				auto iter = range.first;
				for(; iter != range.second; ++iter){
					auto q = raw_customer_orders(iter->second, ele);
					result.push_back(q);
				}
			}
		}
		return result;
	}
	else
	{
		std::unordered_multimap<K2, R2> table2_map;
		for(size_t i = 0; i < p_table2->size(); ++i){
			table2_map.insert(std::make_pair(key_func2((*p_table2)[i]), (*p_table2)[i]));
		}
		for(size_t i = 0; i < p_table1->size(); ++i){
			auto ele = (*p_table1)[i];
			if(table2_map.count(key_func1(ele)) > 0){
				auto range = table2_map.equal_range(key_func1(ele));
				auto iter = range.first;
				for(; iter != range.second; ++iter){
					result.push_back(ele + iter->second);
				}
			}
		}
		return result;
	}
};
#endif
