#include "plan.h"
#include <iostream>
#include <vector>
#include "table.h"
#include "hashjoin.h"
#include "ripplejoin.h"
#include "timer.h"
#include "base.h"
#include <cstring>
using namespace std;
static size_t max_time = 1000000 * 300;
static size_t step_time = 100000;
static size_t total_round = 1;
void test_dynamic(date_t &cond_date1, date_t &cond_date2){
	double value = 0;
	double error_limit = 0;
	std::function<bool(base_raw *)> cond_func_customer= [&](base_raw *input)->bool{\
			if(input->raw_tag == RAW_CUSTOMER){
				auto p = (raw_customer*)input;
				if(std::strcmp(p->c_mktsegment, "BUILDING") == 0)
					return true;
			}
			return false;
		};
	std::function<bool(base_raw *)> cond_func_orders = [&](base_raw *input)->bool{\
			if(input->raw_tag == RAW_ORDERS){
				auto p = (raw_orders*) input;
				if(p->o_orderdate < cond_date1){
					return true;
				}
			}
			return false;
		};
	std::function<bool(base_raw *)> cond_func_lineitem = [&](base_raw *input)->bool{\
			if(input->raw_tag == RAW_LINEITEM){
				auto p = (raw_lineitem*) input;
				if(cond_date2 < p->l_shipdate){
					return true;
				}
			}
			return false;
		};
	const size_t n_rels = 3;
	//dynamic_method
	base_table *p[n_rels];
	p[0] = new table_customer();
	p[1] = new table_orders();
	p[2] = new table_lineitem();
	plan_node node[n_rels];
	//! here to define cond function
	std::array<head_node_type, n_rels> head_node;
	head_node[0].col_id = 3;
	head_node[0].segment = "BUILDING";
	head_node[1].col_id = 3;
	head_node[1].start_date = date_t{1992,1,1};
	head_node[1].end_date = cond_date1;
	head_node[2].col_id = 5;
	head_node[2].start_date = cond_date2;
	head_node[2].end_date = date_t{9999,9,9};
	for(size_t i = 0; i < n_rels; ++i){
		p[i]->load_data();
		p[i]->build_hash();
		node[i].p_table = p[i];
		node[i].cond_func = nullptr;
		node[i].p_head = nullptr; 
	}
	node[0].p_head = &head_node[0];
	node[1].p_head = &head_node[1];
	node[0].cond_func = &cond_func_customer;
	//node[1].cond_func = &cond_func_orders;
	node[1].cond_func = &cond_func_orders;
	//! here to define the head node

	//hash_join
	std::vector<base_raw *> vec_c_p;
	std::vector<base_raw *> vec_o_p;
	std::vector<base_raw *> vec_l_p;
	auto p0 = dynamic_cast<table_customer*>(p[0]);
	auto p1 = dynamic_cast<table_orders *>(p[1]);
	auto p2 = dynamic_cast<table_lineitem *>(p[2]);
	for(auto &x:p0->data)
		vec_c_p.push_back(&x);
	for(auto &x:p1->data)
		vec_o_p.push_back(&x);
	for(auto &x:p2->data)
		vec_l_p.push_back(&x);
	//! hash join no cond
	size_t fullsize = 0;
	{
	timer mytimer;
	auto a = hash_join(vec_c_p, vec_o_p, nullptr, nullptr);
	auto b = hash_join(a, vec_l_p, nullptr, nullptr);
	double count = 0;
	for(auto &x:b){
		auto tmp = static_cast<raw_lineitem* >(x);
		count += tmp->l_extendedprice * (1.0 - tmp->l_discount);
	}
	fullsize = b.size();
	std::cout<<"size: "<<fullsize<<"\tresult: "<<count << "\trunning time is: "<< mytimer.get_elapsed()<<endl;
	}
	//! hash join with cond
	{
	timer mytimer;
	auto a = hash_join(vec_c_p, vec_o_p, nullptr, nullptr);
	auto b = hash_join(a, vec_l_p, nullptr, nullptr);
	double count = 0;
	for(auto &x:b){
		auto tmp = dynamic_cast<raw_lineitem* >(x);
		count += tmp->l_extendedprice * (1.0 - tmp->l_discount);
	}
	std::cout<<"Selectivity: "<<1.0-b.size()*1.0/fullsize<<"\tresult: "<<count << "\trunning time is: "<< mytimer.get_elapsed()<<endl;
	value = count;
	}
	//! here to define the result function
	auto result_func = [&](std::vector<plan_result_type> &input)->double{
		double result = 1;
		for(auto &x:input){
			if(x.p_raw->raw_tag == RAW_CUSTOMER || x.p_raw->raw_tag == RAW_ORDERS)
				result *= x.d;
			else
			{
				raw_lineitem *q = (raw_lineitem *)x.p_raw;
				result *= x.d * q->l_extendedprice * (1 - q->l_discount);
			}
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
	auto full_order = tool::gen_order(3);
	auto order_list = tool::back_trace(full_order, adjacency_list);
	//! define plan
	plan query(result_func);
	query.set_decision(100);
	std::vector<exec_params> result(order_list.size());
	for(size_t i = 0; i < 3; ++i){
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
	std::cout<<x.year<<":"<<x.month<<":"<<x.day<<std::endl;
	test_dynamic(x, x);
	return 0;
}
