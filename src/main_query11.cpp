#include "plan.h"
#include <iostream>
#include <vector>
#include "table.h"
#include "hashjoin.h"
#include "timer.h"
#include "base.h"
#include <cstring>
using namespace std;
static std::function<bool(base_raw *)> cond_func_nation = [&](base_raw *input)->bool{\
			if(input->raw_tag == RAW_NATION){
				auto p = (raw_nation*) input;
				if(strcmp(p->n_name, "GERMANY")== 0 ){
					return true;
				}
			}
			return false;
		};
void do_ripple_join(){
}
void test_dynamic(){
	const size_t n_rels = 3;
	//dynamic_method
	base_table *p[n_rels];
	p[0] = new table_nation();
	p[1] = new table_supplier();
	p[2] = new table_partsupp();
	plan_node node[n_rels];
	//! here to define cond function
	std::array<head_node_type, n_rels> head_node;
	head_node[0].col_id = 2;
	head_node[0].segment = "GERMANY";
	for(size_t i = 0; i < n_rels; ++i){
		p[i]->load_data();
		p[i]->build_hash();
		node[i].p_table = p[i];
		//node[i].p_head = &head_node[i]; 
	}
	//node[0].p_head = &head_node[0];
	//node[2].p_head = nullptr;
	//node[0].cond_func = &cond_func_nation;
	//node[1].cond_func = &cond_func_orders;
	//node[2].cond_func = &cond_func_lineitem;
	//! here to define the head node

	//hash_join
	std::vector<base_raw *> vec_n_p;
	std::vector<base_raw *> vec_s_p;
	std::vector<base_raw *> vec_ps_p;
	auto p0 = dynamic_cast<table_nation*>(p[0]);
	auto p1 = dynamic_cast<table_supplier *>(p[1]);
	auto p2 = dynamic_cast<table_partsupp *>(p[2]);
	for(auto &x:p0->data)
		vec_n_p.push_back(&x);
	for(auto &x:p1->data)
		vec_s_p.push_back(&x);
	for(auto &x:p2->data)
		vec_ps_p.push_back(&x);
	//! hash join no cond
	size_t fullsize = 0;
	{
	timer mytimer;
	auto a = hash_join(vec_n_p, vec_s_p, nullptr, nullptr);
	auto b = hash_join(a, vec_ps_p, nullptr, nullptr);
	double count = 0;
	for(auto &x:b){
		auto tmp = dynamic_cast<raw_partsupp* >(x);
		count += tmp->ps_supplycost * (tmp->ps_availqty);
	}
	fullsize = b.size();
	std::cout<<"size: "<<fullsize<<"\tresult: "<<count << "\trunning time is: "<< mytimer.get_elapsed()<<endl;
	}
	//! hash join with cond
	{
	timer mytimer;
	auto a = hash_join(vec_n_p, vec_s_p, &cond_func_nation, nullptr);
	auto b = hash_join(a, vec_ps_p, nullptr, nullptr);
	double count = 0;
	for(auto &x:b){
		auto tmp = dynamic_cast<raw_partsupp* >(x);
		count += tmp->ps_supplycost * (tmp->ps_availqty);
	}
	std::cout<<"Selectivity: "<<1.0-b.size()*1.0/fullsize<<"\tresult: "<<count << "\trunning time is: "<< mytimer.get_elapsed()<<endl;
	}

	
	//! here to define the result function
	auto result_func = [&](std::vector<plan_result_type> &input)->double{
		double result = 1;
		for(auto &x:input){
			if(x.p_raw->raw_tag == RAW_PARTSUPP)
			{
				raw_partsupp *q = (raw_partsupp *)x.p_raw;
				result *= x.d * q->ps_supplycost * (q->ps_availqty);
			}
			else
				result *= x.d;
		}
		return result;
	};
	//! here to define the adjacency list of relations
	std::vector<std::vector<size_t> > adjacency_list{\
		{1}, \
		{0,2}, \
		{1}\
	};
	//! here to define the orders of execution
	std::vector<std::vector<size_t> > order_list{\
		{0, 1, 2},\
		{1, 0, 2},\
		{1, 2, 0},\
		{2, 1, 0}
	};
	//! define plan
	plan query(result_func);
	query.set_decision(100);
	std::vector<exec_params> result(order_list.size());
	for(size_t i = 0; i < 3; ++i){
		query.insert(node[i]);
	}
	double min_value;
	size_t min_id;
	size_t p_count = 0;
	size_t n_count = 0;
	query.set_adj(adjacency_list);
	for(size_t n = 0; n < 100; ++n)
	{
		for(size_t r = 0; r < order_list.size(); ++r){
			query.set_order(order_list[r]);	
			result[r] = query.run_time_limit(1000000);
			auto tmp = tool::calc_ci(result[r].var, 1.0 * result[r].n_samples, .95);
			if(r == 0){
				min_value = tmp;
				min_id = r;
			}
			else
			{
				if(tmp < min_value){
					min_value = tmp;
					min_id = r;
				}
			}
		}
		size_t d_id;
		double d_value;
		d_value = result[0].decision_value;
		d_id = 0;
		for(size_t i = 1; i < order_list.size(); ++i){
			if(result[i].decision_value < d_value){
				d_value = result[i].decision_value;
				d_id = i;
			}
		}
		if(d_id == min_id)
			p_count++;
		else
			n_count++;
	}
	std::cout<<"positive: "<<p_count<<"\tnegtive: "<<n_count<<std::endl;
	for(size_t i = 0; i < n_rels; ++i){
		delete p[i];
	}
	return;
}

int main(){
	test_dynamic();
	return 0;
}
