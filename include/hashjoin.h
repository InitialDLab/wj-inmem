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
#include "rawdata.h"
using namespace std;
typedef std::function<bool(base_raw *)>* cond_func_type;
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
std::vector<base_raw *> hash_join_run(std::vector<base_raw*> &L, std::vector<base_raw *> &R, bool flag){
	std::vector<base_raw *> result;
	if(L.size() == 0 || R.size() == 0){
		std::cerr<<"EMPTY RELATION"<<std::endl;
		return result;
	}
	auto p = L[0];
	auto q = R[0];
	std::unordered_map<uint32_t, base_raw *> h_map;
	//! region
	if(p->raw_tag == RAW_REGION){
		if(q->raw_tag == RAW_NATION){
			for(auto &x:L){
				auto ele = dynamic_cast<raw_region *>(x);
				h_map.insert(std::make_pair(ele->r_regionkey, x));
			}
			for(auto &x:R){
				auto ele = dynamic_cast<raw_nation *>(x);
				auto range = h_map.equal_range(ele->n_regionkey);
				for(auto iter = range.first; iter != range.second; ++iter){
					if(!flag){
						x->pre = iter->second;
						result.push_back(x);
					}
					else
					{
						iter->second->pre = x;
						result.push_back(iter->second);
					}
				}
			}
		}
		else
			std::cout<<"TWO RELATION CANNOT MATCH"<<std::endl;
	}
	//! left part is nation
	if(p->raw_tag == RAW_NATION){
		std::unordered_map<uint32_t, base_raw *> h_map;
		if(q->raw_tag == RAW_REGION){
			for(auto &x:L){
				auto ele = dynamic_cast<raw_nation *>(x);
				h_map.insert(std::make_pair(ele->n_regionkey, x));
			}
			for(auto &x:R){
				auto ele = dynamic_cast<raw_region *>(x);
				auto range = h_map.equal_range(ele->r_regionkey);
				for(auto iter = range.first; iter != range.second; ++iter){
					if(!flag){
						x->pre = iter->second;
						result.push_back(x);
					}
					else
					{
						iter->second->pre = x;
						result.push_back(iter->second);
					}
				}
			}
		}
		else if(q->raw_tag == RAW_CUSTOMER){
			for(auto &x:L){
				auto ele = dynamic_cast<raw_nation *>(x);
				h_map.insert(std::make_pair(ele->n_nationkey, x));
			}
			for(auto &x:R){
				auto ele = dynamic_cast<raw_customer *>(x);
				auto range = h_map.equal_range(ele->c_nationkey);
				for(auto iter = range.first; iter != range.second; ++iter){
					if(!flag){
						x->pre = iter->second;
						result.push_back(x);
					}
					else
					{
						iter->second->pre = x;
						result.push_back(iter->second);
					}
				}
			}
		}
		else if(q->raw_tag == RAW_SUPPLIER){
			for(auto &x:L){
				auto ele = dynamic_cast<raw_nation *>(x);
				h_map.insert(std::make_pair(ele->n_nationkey, x));
			}
			for(auto &x:R){
				auto ele = dynamic_cast<raw_supplier *>(x);
				auto range = h_map.equal_range(ele->s_nationkey);
				for(auto iter = range.first; iter != range.second; ++iter){
					if(!flag){
						x->pre = iter->second;
						result.push_back(x);
					}
					else
					{
						iter->second->pre = x;
						result.push_back(iter->second);
					}
				}
			}
		}
		else
			std::cout<<"TWO RELATION CANNOT MATCH"<<std::endl;
	}
	//! supplier
	if(p->raw_tag == RAW_SUPPLIER){
		std::unordered_map<uint32_t, base_raw *> h_map;
		if(q->raw_tag == RAW_NATION){
			for(auto &x:L){
				auto ele = dynamic_cast<raw_supplier *>(x);
				h_map.insert(std::make_pair(ele->s_nationkey, x));
			}
			for(auto &x:R){
				auto ele = dynamic_cast<raw_nation *>(x);
				auto range = h_map.equal_range(ele->n_nationkey);
				for(auto iter = range.first; iter != range.second; ++iter){
					if(!flag){
						x->pre = iter->second;
						result.push_back(x);
					}
					else
					{
						iter->second->pre = x;
						result.push_back(iter->second);
					}
				}
			}
		}
		else if(q->raw_tag == RAW_PARTSUPP){
			for(auto &x:L){
				auto ele = dynamic_cast<raw_supplier *>(x);
				h_map.insert(std::make_pair(ele->s_suppkey, x));
			}
			for(auto &x:R){
				auto ele = dynamic_cast<raw_partsupp *>(x);
				auto range = h_map.equal_range(ele->ps_suppkey);
				for(auto iter = range.first; iter != range.second; ++iter){
					if(!flag){
						x->pre = iter->second;
						result.push_back(x);
					}
					else
					{
						iter->second->pre = x;
						result.push_back(iter->second);
					}
				}
			}
		}
		else if(q->raw_tag == RAW_LINEITEM){
			for(auto &x:L){
				auto ele = dynamic_cast<raw_supplier *>(x);
				h_map.insert(std::make_pair(ele->s_suppkey, x));
			}
			for(auto &x:R){
				auto ele = dynamic_cast<raw_lineitem *>(x);
				auto range = h_map.equal_range(ele->l_suppkey);
				for(auto iter = range.first; iter != range.second; ++iter){
					if(!flag){
						x->pre = iter->second;
						result.push_back(x);
					}
					else
					{
						iter->second->pre = x;
						result.push_back(iter->second);
					}
				}
			}
		}
		else
			std::cout<<"TWO RELATION CANNOT MATCH"<<std::endl;
	}
	//! customer
	if(p->raw_tag == RAW_CUSTOMER){
		std::unordered_map<uint32_t, base_raw *>h_map;
		if(q->raw_tag == RAW_NATION){
			for(auto &x:L){
				auto ele = dynamic_cast<raw_customer *>(x);
				h_map.insert(std::make_pair(ele->c_nationkey, x));
			}
			for(auto &x:R){
				auto ele = dynamic_cast<raw_nation *>(x);
				auto range = h_map.equal_range(ele->n_nationkey);
				for(auto iter = range.first; iter != range.second; ++iter){
					if(!flag){
						x->pre = iter->second;
						result.push_back(x);
					}
					else
					{
						iter->second->pre = x;
						result.push_back(iter->second);
					}
				}
			}
		}
		else if(q->raw_tag == RAW_ORDERS){
			for(auto &x:L){
				auto ele = dynamic_cast<raw_customer *>(x);
				h_map.insert(std::make_pair(ele->c_custkey, x));
			}
			for(auto &x:R){
				auto ele = dynamic_cast<raw_orders *>(x);
				auto range = h_map.equal_range(ele->o_custkey);
				for(auto iter = range.first; iter != range.second; ++iter){
					if(!flag){
						x->pre = iter->second;
						result.push_back(x);
					}
					else
					{
						iter->second->pre = x;
						result.push_back(iter->second);
					}
				}
			}
		}
		else
			std::cout<<"TWO RELATION CANNOT MATCH"<<std::endl;
	}
	//! part
	if(p->raw_tag == RAW_PART){
		std::unordered_map<uint32_t, base_raw *>h_map;
		if(q->raw_tag == RAW_PARTSUPP){
			for(auto &x:L){
				auto ele = dynamic_cast<raw_part *>(x);
				h_map.insert(std::make_pair(ele->p_partkey, x));
			}
			for(auto &x:R){
				auto ele = dynamic_cast<raw_partsupp *>(x);
				auto range = h_map.equal_range(ele->ps_partkey);
				for(auto iter = range.first; iter != range.second; ++iter){
					if(!flag){
						x->pre = iter->second;
						result.push_back(x);
					}
					else
					{
						iter->second->pre = x;
						result.push_back(iter->second);
					}
				}
			}
		}
		else if(q->raw_tag == RAW_LINEITEM){
			for(auto &x:L){
				auto ele = dynamic_cast<raw_part *>(x);
				h_map.insert(std::make_pair(ele->p_partkey, x));
			}
			for(auto &x:R){
				auto ele = dynamic_cast<raw_lineitem *>(x);
				auto range = h_map.equal_range(ele->l_partkey);
				for(auto iter = range.first; iter != range.second; ++iter){
					if(!flag){
						x->pre = iter->second;
						result.push_back(x);
					}
					else
					{
						iter->second->pre = x;
						result.push_back(iter->second);
					}
				}
			}
		}
		else
			std::cout<<"TWO RELATION CANNOT MATCH"<<std::endl;
	}
	//! partsupp
	if(p->raw_tag == RAW_PARTSUPP){
		std::unordered_map<uint32_t, base_raw *>h_map;
		if(q->raw_tag == RAW_PART){
			for(auto &x:L){
				auto ele = dynamic_cast<raw_partsupp *>(x);
				h_map.insert(std::make_pair(ele->ps_partkey, x));
			}
			for(auto &x:R){
				auto ele = dynamic_cast<raw_part *>(x);
				auto range = h_map.equal_range(ele->p_partkey);
				for(auto iter = range.first; iter != range.second; ++iter){
					if(!flag){
						x->pre = iter->second;
						result.push_back(x);
					}
					else
					{
						iter->second->pre = x;
						result.push_back(iter->second);
					}
				}
			}
		}
		else if(q->raw_tag == RAW_SUPPLIER){
			for(auto &x:L){
				auto ele = dynamic_cast<raw_partsupp *>(x);
				h_map.insert(std::make_pair(ele->ps_suppkey, x));
			}
			for(auto &x:R){
				auto ele = dynamic_cast<raw_supplier *>(x);
				auto range = h_map.equal_range(ele->s_suppkey);
				for(auto iter = range.first; iter != range.second; ++iter){
					if(!flag){
						x->pre = iter->second;
						result.push_back(x);
					}
					else
					{
						iter->second->pre = x;
						result.push_back(iter->second);
					}
				}
			}
		}
		else if(q->raw_tag == RAW_LINEITEM){
			std::unordered_map<uint32_t, std::unordered_map<uint32_t, base_raw*>* >h_multi_map;
			for(auto &x:L){
				auto ele = dynamic_cast<raw_partsupp *>(x);
				if(h_multi_map.count(ele->ps_partkey) > 0){
					auto range = h_multi_map.equal_range(ele->ps_partkey);
					for(auto iter = range.first; iter != range.second; ++iter){
						auto tmp = iter->second;
						tmp->insert(std::make_pair(ele->ps_suppkey, x));
					}
				}
				else{
					auto tmp = new std::unordered_map<uint32_t, base_raw *>();
					tmp->insert(std::make_pair(ele->ps_suppkey, x));
					h_multi_map.insert(std::make_pair(ele->ps_partkey, tmp));
				}
			}
			for(auto &x:R){
				auto ele = dynamic_cast<raw_lineitem *>(x);
				auto range = h_multi_map.equal_range(ele->l_partkey);
				for(auto iter = range.first; iter != range.second; ++iter){
					auto tmp = iter->second;
					auto range_suppkey = tmp->equal_range(ele->l_suppkey);
					for(auto iter_suppkey = range_suppkey.first; iter_suppkey != range_suppkey.second; ++iter_suppkey){
						if(!flag){
							x->pre = iter_suppkey->second;
							result.push_back(x);
						}
						else
						{
							iter_suppkey->second->pre = x;
							result.push_back(iter_suppkey->second);
						}
					}
				}
			}
			for(auto iter = h_multi_map.begin(); iter != h_multi_map.end(); ++iter)
				delete iter->second;
		}
		else
			std::cout<<"TWO RELATION CANNOT MATCH"<<std::endl;
	}
	//! orders
	if(p->raw_tag == RAW_ORDERS){
		std::unordered_map<uint32_t, base_raw *>h_map;
		if(q->raw_tag == RAW_CUSTOMER){
			for(auto &x:L){
				auto ele = dynamic_cast<raw_orders *>(x);
				h_map.insert(std::make_pair(ele->o_custkey, x));
			}
			for(auto &x:R){
				auto ele = dynamic_cast<raw_customer *>(x);
				auto range = h_map.equal_range(ele->c_custkey);
				for(auto iter = range.first; iter != range.second; ++iter){
					if(!flag){
						x->pre = iter->second;
						result.push_back(x);
					}
					else
					{
						iter->second->pre = x;
						result.push_back(iter->second);
					}
				}
			}
		}
		else if(q->raw_tag == RAW_LINEITEM){
			for(auto &x:L){
				auto ele = dynamic_cast<raw_orders *>(x);
				h_map.insert(std::make_pair(ele->o_orderkey, x));
			}
			for(auto &x:R){
				auto ele = dynamic_cast<raw_lineitem *>(x);
				auto range = h_map.equal_range(ele->l_orderkey);
				for(auto iter = range.first; iter != range.second; ++iter){
					if(!flag){
						x->pre = iter->second;
						result.push_back(x);
					}
					else
					{
						iter->second->pre = x;
						result.push_back(iter->second);
					}
				}
			}
		}
		else
			std::cout<<"TWO RELATION CANNOT MATCH"<<std::endl;
	}
	if(p->raw_tag == RAW_LINEITEM){
		std::unordered_map<uint32_t, base_raw *>h_map;
		if(q->raw_tag == RAW_PARTSUPP){
			std::unordered_map<uint32_t, std::unordered_map<uint32_t, base_raw*>* >h_multi_map;
			for(auto &x:L){
				auto ele = dynamic_cast<raw_lineitem *>(x);
				if(h_multi_map.count(ele->l_partkey) > 0){
					auto range = h_multi_map.equal_range(ele->l_partkey);
					for(auto iter = range.first; iter != range.second; ++iter){
						auto tmp = iter->second;
						tmp->insert(std::make_pair(ele->l_suppkey, x));
					}
				}
				else{
					auto tmp = new std::unordered_map<uint32_t, base_raw *>();
					tmp->insert(std::make_pair(ele->l_suppkey, x));
					h_multi_map.insert(std::make_pair(ele->l_partkey, tmp));
				}
			}
			for(auto &x:R){
				auto ele = dynamic_cast<raw_partsupp *>(x);
				auto range = h_multi_map.equal_range(ele->ps_partkey);
				for(auto iter = range.first; iter != range.second; ++iter){
					auto tmp = iter->second;
					auto range_suppkey = tmp->equal_range(ele->ps_suppkey);
					for(auto iter_suppkey = range_suppkey.first; iter_suppkey != range_suppkey.second; ++iter_suppkey){
						if(!flag){
							x->pre = iter_suppkey->second;
							result.push_back(x);
						}
						else
						{
							iter_suppkey->second->pre = x;
							result.push_back(iter_suppkey->second);
						}
					}
				}
			}
			for(auto iter = h_multi_map.begin(); iter != h_multi_map.end(); ++iter)
				delete iter->second;
		}
		else if(q->raw_tag == RAW_ORDERS){
			for(auto &x:L){
				auto ele = dynamic_cast<raw_lineitem *>(x);
				h_map.insert(std::make_pair(ele->l_orderkey, x));
			}
			for(auto &x:R){
				auto ele = dynamic_cast<raw_orders *>(x);
				auto range = h_map.equal_range(ele->o_orderkey);
				for(auto iter = range.first; iter != range.second; ++iter){
					if(!flag){
						x->pre = iter->second;
						result.push_back(x);
					}
					else
					{
						iter->second->pre = x;
						result.push_back(iter->second);
					}
				}
			}
		}
		else if(q->raw_tag == RAW_PART){
			for(auto &x:L){
				auto ele = dynamic_cast<raw_lineitem *>(x);
				h_map.insert(std::make_pair(ele->l_partkey, x));
			}
			for(auto &x:R){
				auto ele = dynamic_cast<raw_part *>(x);
				auto range = h_map.equal_range(ele->p_partkey);
				for(auto iter = range.first; iter != range.second; ++iter){
					if(!flag){
						x->pre = iter->second;
						result.push_back(x);
					}
					else
					{
						iter->second->pre = x;
						result.push_back(iter->second);
					}
				}
			}
		}
		else if(q->raw_tag == RAW_SUPPLIER){
			for(auto &x:L){
				auto ele = dynamic_cast<raw_lineitem *>(x);
				h_map.insert(std::make_pair(ele->l_suppkey, x));
			}
			for(auto &x:R){
				auto ele = dynamic_cast<raw_supplier *>(x);
				auto range = h_map.equal_range(ele->s_suppkey);
				for(auto iter = range.first; iter != range.second; ++iter){
					if(!flag){
						x->pre = iter->second;
						result.push_back(x);
					}
					else
					{
						iter->second->pre = x;
						result.push_back(iter->second);
					}
				}
			}
		}
		else
			std::cout<<"TWO RELATION CANNOT MATCH"<<std::endl;
	}
	return result;
}
std::vector<base_raw *> hash_join(std::vector<base_raw *> &L, std::vector<base_raw *> &R, cond_func_type func_l,
		cond_func_type func_r){
	std::vector<base_raw *> NL;
	std::vector<base_raw *> NR;
	std::vector<base_raw *> result;
	std::vector<base_raw *>* p_small, *p_large, *p_l, *p_r;
	if(func_l != nullptr){
		for(auto &x:L){
			if((*func_l)(x))
				NL.push_back(x);
		}
		p_l = &NL;
	}
	else
		p_l = &L;
	if(func_r != nullptr){
		for(auto &x:R){
			if((*func_r)(x))
				NR.push_back(x);
		}
		p_r = &NR;
	}
	else
		p_r = &R;
	p_small = p_l->size() <= p_r->size()? p_l:p_r;
	p_large = p_r->size() >= p_l->size()? p_r:p_l;
	//! save whether the position is changed
	bool flag = false;
	if(p_small != p_l)
		flag = true;
	result = hash_join_run(*p_small, *p_large, flag);
	return result;
}
//! Here are filtering, the left part is the small part. build unordered_hash map on the left part.
#endif
