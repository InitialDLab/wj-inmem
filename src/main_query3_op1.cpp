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
static size_t max_time = 1000000 * 3;
static size_t step_time = 100000;
void test_dynamic(date_t &cond_date1, date_t &cond_date2){
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
	base_table *p[3];
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
	//! here to define the head node
	node[0].p_head = &head_node[0];
	node[1].p_head = &head_node[1];
	//node[2].p_head = &head_node[2];
	node[0].cond_func = &cond_func_customer;
	node[1].cond_func = &cond_func_orders;
	//node[2].cond_func = &cond_func_lineitem;

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
	double value = 0;
	size_t fullsize = 0;
	{
	timer mytimer;
	auto a = hash_join(vec_c_p, vec_o_p, nullptr, nullptr);
	auto b = hash_join(a, vec_l_p, nullptr, nullptr);
	double count = 0;
	for(auto &x:b){
		auto tmp = dynamic_cast<raw_lineitem* >(x);
		count += tmp->l_extendedprice * (1.0 - tmp->l_discount);
	}
	fullsize = b.size();
	value = count;
	std::cout<<"size: "<<fullsize<<"\tresult: "<<count << "\trunning time is: "<< mytimer.get_elapsed()<<endl;
	}
	//! hash join with cond
	{
	timer mytimer;
	auto a = hash_join(vec_c_p, vec_o_p, &cond_func_customer, &cond_func_orders);
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
	query.set_adj(adjacency_list);
	//no op part
	size_t total_round = 1;
	for(size_t n = 0; n < total_round; ++n)
	{
		std::cout<<"round: "<<n<<std::endl;
		timer walk_timer;
		walk_plan_params opt;
		double d_value = 0;
		size_t best_order = 0;
		for(size_t r = 0; r < order_list.size(); ++r){
			//std::cout<<"order: "<<r<<std::endl;
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
			std::cout<<"order: "<<r<<std::endl;
			//query.run_time_limit(max_time, step_time);
		}
		query.set_order(order_list[best_order]);
		query.run_ci_limit(value * 0.01, opt);
	}
	//op part
	query.set_agg_id(2);
	query.set_agg_limit(4 * vec_l_p.size()/vec_o_p.size());
	//query.set_agg_limit(10);
	std::cout<<"exact value: "<<value<<std::endl;
	for(size_t n = 0; n < 1; ++n)
	{
		std::cout<<"round: "<<n<<std::endl; 
		for(size_t r = 0; r < order_list.size(); ++r){
			std::cout<<"order:  "<<r<<"\t";
			query.set_order(order_list[r]);	
			//result[r] = query.run_time_limit_op1(max_time, step_time);
			result[r] = query.run_ci_limit_op1(0.01*value);
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
		std::cout<<"decision id: "<<d_id<<std::endl;
		std::cout<<"best id: "<<min_id<<std::endl;
	}
	for(size_t i = 0; i < 3; ++i)
		delete p[i];
	return;
}
int main(int argc, char* argv[]){
	date_t x, y;
	tool::read_date_t(x, y);
	std::cout<<x.year<<":"<<x.month<<":"<<x.day<<std::endl;
	std::cout<<y.year<<":"<<y.month<<":"<<y.day<<std::endl;
	test_dynamic(x, y);
	return 0;
}
