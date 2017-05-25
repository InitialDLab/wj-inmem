#include "plan.h"
#include <iostream>
#include <vector>
#include "table.h"
#include "hashjoin.h"
#include "timer.h"
#include "base.h"
#include "ripplejoin.h"
#include <cstring>
static size_t max_time = 1000000 * 30;
static size_t step_time = 10000;
const size_t total_round = 1;
void test_dynamic(date_t &cond_date1, date_t &cond_date2){
	std::function<bool(base_raw *)> cond_func_lineitem = [&](base_raw *input)->bool{\
			if(input->raw_tag == RAW_LINEITEM){
				auto p = (raw_lineitem*) input;
				if(p->l_shipdate < cond_date2 && cond_date1 < p->l_shipdate){
					return true;
				}
			}
			return false;
		};
	const size_t n_rels = 4;
	//dynamic_method
	base_table *p[n_rels];
	p[0] = new table_customer();
	p[1] = new table_orders();
	p[2] = new table_lineitem();
	p[3] = new table_supplier();
	plan_node node[n_rels];
	//! here to define cond function
	std::array<head_node_type, n_rels> head_node;
	head_node[0].col_id = 2;
	head_node[0].key = 18;
	head_node[2].col_id = 5;
	head_node[2].start_date = cond_date1;
	head_node[2].end_date = cond_date2;
	head_node[3].col_id = 2;
	head_node[3].key = 18;
	for(size_t i = 0; i < n_rels; ++i){
		p[i]->load_data();
		p[i]->build_hash();
		node[i].p_table = p[i];
		node[i].p_head = nullptr;
		//node[i].p_head = &head_node[i]; 
	}
	//node[2].p_head = &head_node[2];
	//node[3].p_head = &head_node[3];
	//node[2].cond_func = &cond_func_lineitem;
	//node[3].cond_func = &cond_func_lineitem;
	//! here to define the head node

	//hash_join
	std::vector<base_raw *> vec_nc_p;
	std::vector<base_raw *> vec_c_p;
	std::vector<base_raw *> vec_o_p;
	std::vector<base_raw *> vec_l_p;
	std::vector<base_raw *> vec_s_p;
	std::vector<base_raw *> vec_ns_p;
	auto p0 = dynamic_cast<table_customer*>(p[0]);
	auto p1 = dynamic_cast<table_orders *>(p[1]);
	auto p2 = dynamic_cast<table_lineitem *>(p[2]);
	auto p3 = dynamic_cast<table_supplier *>(p[3]);
	auto vec_nc = load_nation();
	auto vec_ns = load_nation();
	for(auto &x:p0->data)
		vec_c_p.push_back(&x);
	for(auto &x:p1->data)
		vec_o_p.push_back(&x);
	for(auto &x:p2->data)
		vec_l_p.push_back(&x);
	for(auto &x:p3->data)
		vec_s_p.push_back(&x);
	for(auto &x:vec_nc)
		vec_nc_p.push_back(&x);
	for(auto &x:vec_ns)
		vec_ns_p.push_back(&x);
	//! hash join no cond
	size_t fullsize = 0;
	{
	timer mytimer;
	auto b = hash_join(vec_c_p, vec_o_p, nullptr, nullptr);
	auto c = hash_join(b, vec_l_p, nullptr, nullptr);
	auto d = hash_join(c, vec_s_p, nullptr, nullptr);
	double count = 0;
	for(auto &x:d){
		auto p = x;
		while(p->raw_tag != RAW_LINEITEM)
			p = p->pre;
		auto tmp = dynamic_cast<raw_lineitem* >(p);
		count += tmp->l_extendedprice * (1.0 - tmp->l_discount);
	}
	fullsize = d.size();
	std::cout<<"size: "<<fullsize<<"\tresult: "<<count << "\trunning time is: "<< mytimer.get_elapsed()<<endl;
	}
	//! hash join with cond
	double value = 0;
	{
	timer mytimer;
	auto b = hash_join(vec_c_p, vec_o_p, nullptr, nullptr);
	auto c = hash_join(b, vec_l_p, nullptr, nullptr);
	auto d = hash_join(c, vec_s_p, nullptr, nullptr);
	double count = 0;
	for(auto &x:d){
		auto p = x;
		while(p->raw_tag != RAW_LINEITEM)
			p = p->pre;
		auto tmp = dynamic_cast<raw_lineitem* >(p);
		count += tmp->l_extendedprice * (1.0 - tmp->l_discount);
	}
	value = count;
	std::cout<<"Selectivity: "<<1.0-d.size()*1.0/fullsize<<"\tresult: "<<count << "\trunning time is: "<< mytimer.get_elapsed()<<endl;
	}
	
	//! here to define the result function
	auto result_func = [&](std::vector<plan_result_type> &input)->double{
		double result = 1;
		for(auto &x:input){
			if(x.p_raw->raw_tag == RAW_LINEITEM)
			{
				raw_lineitem *q = (raw_lineitem *)x.p_raw;
				result *= x.d * q->l_extendedprice * (1 - q->l_discount);
			}
			else
				result *= x.d;
		}
		return result;
	};
	//! here to define the adjacency list of relations
	//! here to define the graph of relations
	std::vector<std::vector<size_t> > adjacency_list{\
		{1}, \
		{0, 2}, \
		{1,3}, \
		{2}
	};
	//! here to define the orders of execution
	auto full_order = tool::gen_order(4);
	auto order_list = tool::back_trace(full_order, adjacency_list);
	//! define plan
	plan query(result_func);
	query.set_decision(100);
	std::vector<exec_params> result(order_list.size());
	for(size_t i = 0; i < n_rels; ++i){
		query.insert(node[i]);
	}
	query.set_adj(adjacency_list);
	for(size_t n = 0; n < total_round; ++n)
	{
		std::cout<<"round: "<<n<<std::endl;
		timer walk_timer;
		walk_plan_params opt;
		double d_value = 0;
		size_t best_order = 0;
		for(size_t r = 0; r < order_list.size(); ++r){
			query.set_order(order_list[r]);	
			auto tmp = query.walk_plan_optimizer();	
			opt.sum_y += tmp.sum_y;
			opt.sum_y2 += tmp.sum_y2;
			opt.n_samples += tmp.n_samples;
			if(r == 0){
				d_value = tmp.decision_value;
			}
			else
			{
				if(tmp.decision_value < d_value){
					d_value = tmp.decision_value;
					best_order = r;
				}
			}
		}
		std::cout<<"optimizer time: "<<walk_timer.get_elapsed()<<std::endl;
		for(size_t r = 0; r < order_list.size(); ++r){
			query.set_order(order_list[r]);	
			query.run_ci_limit(value*0.01);
		}
		query.set_order(order_list[best_order]);
		query.run_ci_limit(value * 0.01, opt);
	}
	for(size_t i = 0; i < n_rels; ++i){
		delete p[i];
	}
	return;
}
int main(int argc, char* argv[]){
	if(argc < 1)
		return 0;
	std::vector<date_t> cond_date{{1998,12,1}, {1998,6,1}, {1997,12,1}, {1997,6,1},\
		{1996,12,1}, {1996, 6, 1},{1995,12,1}, {1995,6,1}, {1994,12,1}, {1994,6,1}, {1993,12,1}, {1993, 6,1},\
		{1992,12,1}};
	auto i = std::strtod(argv[1], NULL);
	auto x = cond_date[i];
	date_t date_start = {1992,1,1};
	std::cout<<x.year<<":"<<x.month<<":"<<x.day<<std::endl;
	test_dynamic(date_start, x);
	return 0;
}
