#ifndef _H_PLAN
#define _H_PLAN
#include "table.h" 
#include "base.h"
#include "rawdata.h"
#include "tool.h"
#include "timer.h"
#include <unordered_set>
#include <tuple>
//! value, var, step, #samples, #rej_join, #rej_cond
struct exec_params{
	double est;
	double var;
	size_t sum_step;
	size_t n_samples;
	size_t n_rejected_join;
	size_t n_rejected_cond;
	double decision_value;
	size_t decision_time;
	exec_params(){};
	exec_params(double e, double v, size_t s, size_t n_s, size_t n_r_j, size_t n_r_c){
		est = e;
		var = v;
		sum_step = s;
		n_samples = n_s;
		n_rejected_join = n_r_j;
		n_rejected_cond = n_r_c;
	}
	exec_params(double e, double v, size_t s, size_t n_s, size_t n_r_j, size_t n_r_c, double d_value){
		est = e;
		var = v;
		sum_step = s;
		n_samples = n_s;
		n_rejected_join = n_r_j;
		n_rejected_cond = n_r_c;
		decision_value = d_value;
	}
	exec_params(double e, double v, size_t s, size_t n_s, size_t n_r_j, size_t n_r_c, double d_value, size_t d_time){
		est = e;
		var = v;
		sum_step = s;
		n_samples = n_s;
		n_rejected_join = n_r_j;
		n_rejected_cond = n_r_c;
		decision_value = d_value;
		decision_time = d_time;
	}
};
inline std::ostream &operator<<(std::ostream &out, exec_params &x){
	out<<x.est<<" "<<x.var<<" "<<x.sum_step<<" "<<x.n_samples<<" "<<x.n_rejected_join<<" \
		"<<x.n_rejected_cond<<std::endl;
	return out;
}
struct plan_node{
	base_table *p_table;
	cond_func_type cond_func;
	agg_func_type agg_func;
	head_node_type *p_head;
	plan_node(){
		p_head = nullptr;
		cond_func = nullptr;
	};
};
class plan{
	public:
		plan(result_func_type func){
			n_rejected_join = 0;
			n_rejected_cond = 0;
			n_samples = 0;
			result_func = func;
			prob = 0.95;
			sum_step = 0;
		}
		void insert(plan_node &t){
			data.push_back(t);
		}
		void set_adj(std::vector<std::vector<size_t>> &t){
			adjacency_list.clear();
			for(auto &x:t)
				adjacency_list.emplace_back(x);
		}
		std::vector<std::vector<size_t> > gen_order(){
			std::vector<std::vector<size_t> > result;
			size_t num = data.size();
			auto all = tool::full_permutation(num);
			for(auto &x : all){
				std::vector<bool> tag(num, false);
				for(size_t i = 0; i < x.size(); ++i){
					if(i == 0){
						tag[x[i]] = true;
						continue;
					}
					for(auto &y: adjacency_list[x[i]]){
						if(tag[y]){
							break;
						}
					}
				}
			}
			return result;
		}
		void set_order(std::vector<size_t> &t){
			order.clear();
			for(auto &x:t){
				order.push_back(x);
			}
		}
		double single_run(){
			auto num_rels = data.size();
			std::unordered_set<size_t> rel_map;
			std::vector<plan_result_type> raw_list(num_rels);
			auto ind = order[0];
			if(data[ind].p_head == nullptr){
				first_node_sample_all(data[ind], raw_list, ind);
				rel_map.insert(ind);
			}
			else{
				first_node_sample_range(data[ind], raw_list, ind);
				rel_map.insert(ind);
			}
			sum_step++;
			for(size_t i = 1; i < num_rels; ++i){
				ind = order[i];
				base_table *p = data[ind].p_table;
				for(auto &j:adjacency_list[ind]){
					if(rel_map.count(j) > 0){
						plan_result_type tmp = p->search(raw_list[j].p_raw);
						if(tmp.d == 0){
							n_rejected_join++;
							return 0.0;
						}
						if(data[ind].cond_func != nullptr)
						{
							if((*data[ind].cond_func)(tmp.p_raw))
								raw_list[ind] = tmp;
							else
							{
								n_rejected_cond++;
								return 0.0;
							}
						}
						else
							raw_list[ind] = tmp;
						rel_map.insert(ind);
						break;
					}
				}
				sum_step++;
			}
			return result_func(raw_list);
		}
		double single_run_op1(){
			auto num_rels = data.size();
			double est = 0;
			std::unordered_set<size_t> rel_map;
			std::vector<plan_result_type> raw_list(num_rels);
			auto ind = order[0];
			if(data[ind].p_head == nullptr){
				first_node_sample_all(data[ind], raw_list, ind);
				if(ind == agg_id){
					raw_lineitem *ele = static_cast<raw_lineitem*>(raw_list[ind].p_raw);
					est =  ele->l_extendedprice * (1 - ele->l_discount) * raw_list[ind].d;
				}
				else
					est = raw_list[ind].d;
				rel_map.insert(ind);
			}
			else{
				first_node_sample_range(data[ind], raw_list, ind);
				if(ind == agg_id){
					raw_lineitem *ele = static_cast<raw_lineitem*>(raw_list[ind].p_raw);
					est =  ele->l_extendedprice * (1 - ele->l_discount) * raw_list[ind].d;
				}
				else
					est = raw_list[ind].d;
				rel_map.insert(ind);
			}
			sum_step++;
			for(size_t i = 1; i < num_rels; ++i){
				ind = order[i];
				base_table *p = data[ind].p_table;
				for(auto &j:adjacency_list[ind]){
					if(rel_map.count(j) > 0){
						if(ind != agg_id){
							plan_result_type tmp = p->search(raw_list[j].p_raw);
							if(tmp.d == 0){
								n_rejected_join++;
								return 0.0;
							}
							if(data[ind].cond_func != nullptr)
							{
								if((*data[ind].cond_func)(tmp.p_raw))
									raw_list[ind] = tmp;
								else
								{
									n_rejected_cond++;
									return 0.0;
								}
							}
							else
								raw_list[ind] = tmp;
							rel_map.insert(ind);
							sum_step++;
							est *= raw_list[ind].d;
							break;
						}
						else{
							auto q = static_cast<table_lineitem *>(p);
							if(q->get_size(raw_list[j].p_raw) > agg_limit){
								plan_result_type tmp = p->search(raw_list[j].p_raw);
								raw_list[ind] = tmp;
								if(data[ind].cond_func != nullptr)
								{
									if((*data[ind].cond_func)(tmp.p_raw))
										raw_list[ind] = tmp;
									else
									{
										n_rejected_cond++;
										return 0.0;
									}
								}
								else
									raw_list[ind] = tmp;
								raw_lineitem *ele = static_cast<raw_lineitem*>(raw_list[ind].p_raw);
								est *=  ele->l_extendedprice * (1 - ele->l_discount) * raw_list[ind].d;
								rel_map.insert(ind);
								sum_step++;
								break;
							}
							else{
								plan_result_type tmp = q->search_all(raw_list[j].p_raw);
								if(tmp.d == 0){
									n_rejected_join++;
									return 0.0;
								}
								auto tp = tmp.p_raw;
								double acm = 0;
								if(data[ind].cond_func != nullptr)
								{
									for(size_t k = 0; k < tmp.len; ++k){
										if((*data[ind].cond_func)(tp)){
											raw_list[ind].p_raw = tp;
											raw_lineitem *ele = static_cast<raw_lineitem*>(raw_list[ind].p_raw);
											acm +=  ele->l_extendedprice * (1 - ele->l_discount);
										}
										tp = increment(tp);
									}
									est *= acm;
								}
								else
								{
									for(size_t k = 0; k < tmp.len; ++k){
										raw_list[ind].p_raw = tp;
										raw_lineitem *ele = static_cast<raw_lineitem*>(raw_list[ind].p_raw);
										acm +=  ele->l_extendedprice * (1 - ele->l_discount);
										tp = increment(tp);
									}
									est *= acm;
								}
								rel_map.insert(ind);
								sum_step += q->get_size(raw_list[j].p_raw);
								break;
							}
						}
					}
				}
			}
			return est;
		}
		double single_run_op2(){
			auto num_rels = data.size();
			double pho = 0;
			double local_d = 0;
			std::unordered_set<size_t> rel_map;
			std::vector<plan_result_type> raw_list(num_rels);
			auto ind = order[0];
			timer level_timer();
			if(data[ind].p_head == nullptr){
				first_node_sample_all(data[ind], raw_list, ind);
				num_sel[ind] = raw_list[ind].d;
				rel_map.insert(ind);
			}
			else{
				first_node_sample_range(data[ind], raw_list, ind);
				num_sel[ind] = raw_list[ind].d;
				rel_map.insert(ind);
			}
			++num_level_sample[ind];
			plan_result_type tmp, tmp_op;
			sum_step++;
			for(size_t i = 1; i < num_rels; ++i){
				ind = order[i];
				++num_level_sample[ind];
				base_table *p = data[ind].p_table;
				for(auto &j:adjacency_list[ind]){
					if(rel_map.count(j) > 0){
						if(ind != agg_id){
							//here to do decide whether we can do optimization 2
							if(data[ind].cond_func != nullptr)
							{
								local_d = p->get_size(raw_list[j].p_raw);
								if(local_d == 0){
									return 0.0;
								}
								//op2
								/*
								pho = 1;
								for(size_t k = 1; k < i; ++k){
									pho *= static_cast<double>(num_level_success[order[k]]) / num_level_sample[order[k]];
								}
								*/
								pho = static_cast<double>(num_level_success[ind]) / num_level_sample[ind];
								//std::cout<<pho<<"|"<<local_d<<"|"<<local_d*pho<<"|"<<i<<std::endl;
								if(pho > 0 && local_d * pho < (i+1)){
									tmp_op = p->search_cond(raw_list[j].p_raw, data[ind].cond_func);
									num_level_sample[ind] += local_d  -1;
									if(tmp_op.d == 0){
										n_rejected_join++;
										return 0.0;
									}
									raw_list[ind] = tmp_op;
									num_level_success[ind] += tmp_op.d;
									num_level_sample[ind] += tmp_op.d - 1;
									rel_map.insert(ind);
									break;
								}
								else
								{
									++num_no_op2[ind];
									plan_result_type tmp = p->search(raw_list[j].p_raw);
									if(tmp.d == 0){
										n_rejected_join++;
										return 0.0;
									}
									if((*data[ind].cond_func)(tmp.p_raw)){
										++num_level_success[ind];
										raw_list[ind] = tmp;
									}
									else
									{
										n_rejected_cond++;
										return 0.0;
									}
									rel_map.insert(ind);
									break;
								}
							}
							else
							{
								plan_result_type tmp = p->search(raw_list[j].p_raw);
								if(tmp.d == 0){
									n_rejected_join++;
									return 0.0;
								}
								raw_list[ind] = tmp;
								rel_map.insert(ind);
								break;
							}
							rel_map.insert(ind);
							break;
						}
						else
						{
							plan_result_type tmp = p->search(raw_list[j].p_raw);
							if(tmp.d == 0){
								n_rejected_join++;
								return 0.0;
							}
							if(data[ind].cond_func != nullptr)
							{
								if((*data[ind].cond_func)(tmp.p_raw))
									raw_list[ind] = tmp;
								else
								{
									n_rejected_cond++;
									return 0.0;
								}
							}
							else
								raw_list[ind] = tmp;
							rel_map.insert(ind);
							break;
						};
					}
				}
				sum_step++;
			}
			return result_func(raw_list);
		}
		double single_run_op(){
			auto num_rels = data.size();
			double est = 0;
			double pho = 0;
			double local_d = 0;
			std::unordered_set<size_t> rel_map;
			std::vector<plan_result_type> raw_list(num_rels);
			auto ind = order[0];
			timer level_timer();
			if(data[ind].p_head == nullptr){
				first_node_sample_all(data[ind], raw_list, ind);
				num_sel[ind] = raw_list[ind].d;
				if(ind == agg_id){
					raw_lineitem *ele = static_cast<raw_lineitem*>(raw_list[ind].p_raw);
					est =  ele->l_extendedprice * (1 - ele->l_discount) * raw_list[ind].d;
				}
				else
					est = raw_list[ind].d;
				rel_map.insert(ind);
			}
			else{
				first_node_sample_range(data[ind], raw_list, ind);
				num_sel[ind] = raw_list[ind].d;
				if(ind == agg_id){
					raw_lineitem *ele = static_cast<raw_lineitem*>(raw_list[ind].p_raw);
					est =  ele->l_extendedprice * (1 - ele->l_discount) * raw_list[ind].d;
				}
				else
					est = raw_list[ind].d;
				rel_map.insert(ind);
			}
			++num_level_sample[ind];
			plan_result_type tmp, tmp_op;
			sum_step++;
			for(size_t i = 1; i < num_rels; ++i){
				ind = order[i];
				++num_level_sample[ind];
				base_table *p = data[ind].p_table;
				for(auto &j:adjacency_list[ind]){
					if(rel_map.count(j) > 0){
						if(ind != agg_id){
							//here to do decide whether we can do optimization 2
							if(data[ind].cond_func != nullptr){
								local_d = p->get_size(raw_list[j].p_raw);
								if(local_d == 0){
									return 0.0;
								}
								//op2
								pho = static_cast<double>(num_level_success[ind]) / num_level_sample[ind];
								if(pho > 0 && local_d * pho < (i+1) && local_d > 1){
									tmp_op = p->search_cond(raw_list[j].p_raw, data[ind].cond_func);
									num_level_sample[ind] += local_d  -1;
									if(tmp_op.d == 0){
										n_rejected_join++;
										return 0.0;
									}
									raw_list[ind] = tmp_op;
									num_level_success[ind] += tmp_op.d;
									num_level_sample[ind] += tmp_op.d - 1;
									est *= tmp_op.d;
								}
								else
								{
									++num_no_op2[ind];
									tmp = p->search(raw_list[j].p_raw);
									if(tmp.d == 0){
										n_rejected_join++;
										return 0.0;
									}
									if((*data[ind].cond_func)(tmp.p_raw)){
										++num_level_success[ind];
										raw_list[ind] = tmp;
										est *= tmp.d;
									}
									else
									{
										n_rejected_cond++;
										return 0.0;
									}
								}
							}
							else
							{
								tmp = p->search(raw_list[j].p_raw);
								if(tmp.d == 0){
									n_rejected_join++;
									return 0.0;
								}
								raw_list[ind] = tmp;
								est *= tmp.d;
							}
							rel_map.insert(ind);
							break;
						}
						//ind == agg_id
						else
						{
							//agg_id last table 
							//use op1
							if(agg_id == order[num_rels - 1]){
								auto q = static_cast<table_lineitem *>(p);
								if(q->get_size(raw_list[j].p_raw) > agg_limit){
									plan_result_type tmp = p->search(raw_list[j].p_raw);
									raw_list[ind] = tmp;
									if(data[ind].cond_func != nullptr)
									{
										if((*data[ind].cond_func)(tmp.p_raw))
											raw_list[ind] = tmp;
										else
										{
											n_rejected_cond++;
											return 0.0;
										}
									}
									else
										raw_list[ind] = tmp;
									raw_lineitem *ele = static_cast<raw_lineitem*>(raw_list[ind].p_raw);
									est *=  ele->l_extendedprice * (1 - ele->l_discount) * raw_list[ind].d;
									rel_map.insert(ind);
									sum_step++;
									break;
								}
								else{
									plan_result_type tmp = q->search_all(raw_list[j].p_raw);
									if(tmp.d == 0){
										n_rejected_join++;
										return 0.0;
									}
									auto tp = tmp.p_raw;
									double acm = 0;
									if(data[ind].cond_func != nullptr)
									{
										for(size_t k = 0; k < tmp.len; ++k){
											if((*data[ind].cond_func)(tp)){
												raw_list[ind].p_raw = tp;
												raw_lineitem *ele = static_cast<raw_lineitem*>(raw_list[ind].p_raw);
												acm +=  ele->l_extendedprice * (1 - ele->l_discount);
											}
											tp = increment(tp);
										}
										est *= acm;
									}
									else
									{
										for(size_t k = 0; k < tmp.len; ++k){
											raw_list[ind].p_raw = tp;
											raw_lineitem *ele = static_cast<raw_lineitem*>(raw_list[ind].p_raw);
											acm +=  ele->l_extendedprice * (1 - ele->l_discount);
											tp = increment(tp);
										}
										est *= acm;
									}
									rel_map.insert(ind);
									//sum_step += q->get_size(raw_list[j].p_raw);
									break;
								}
							}
							else{
									plan_result_type tmp = p->search(raw_list[j].p_raw);
									raw_list[ind] = tmp;
									if(data[ind].cond_func != nullptr)
									{
										if((*data[ind].cond_func)(tmp.p_raw))
											raw_list[ind] = tmp;
										else
										{
											n_rejected_cond++;
											return 0.0;
										}
									}
									else
										raw_list[ind] = tmp;
									raw_lineitem *ele = static_cast<raw_lineitem*>(raw_list[ind].p_raw);
									est *=  ele->l_extendedprice * (1 - ele->l_discount) * raw_list[ind].d;
									rel_map.insert(ind);
									break;
							}
						};
					}
				}
				sum_step++;
			}
			return est;
		}
		exec_params run(const size_t &round){
			auto num_rels = data.size();
			if(num_rels == 0){
				std::cerr<<"empty plan"<<std::endl;
				return exec_params(0.0, 0.0, 0, 0, 0, 0);
			}
			double y = 0;
			double sum_y = 0;
			double sum_y2 = 0;
			sum_step = 0;
			double ci = 0;
			n_samples = n_rejected_join = n_rejected_cond = sum_step = 0;
			tool::start_report();
			exec_timer.restart();
			size_t cur_time = 0;
			size_t old_time = 0;
			while(n_samples < round){
				n_samples++;
				auto est = single_run();
				sum_y += est;
				y = sum_y / (n_samples);
				sum_y2 += est * est;
				cur_time = exec_timer.get_elapsed();
				if((cur_time - old_time) > 10000 && cur_time > 0){
					ci = tool::calc_ci(sum_y2, sum_y, y, n_samples, prob);
					tool::report(cur_time, n_samples, n_rejected_join, n_rejected_cond, double(n_rejected_join + n_rejected_cond)/ n_samples, y, ci, prob);
					old_time = cur_time;
				}
			}
			ci = tool::calc_ci(sum_y2, sum_y, y, n_samples, prob);
			auto var = tool::calc_var(sum_y2, sum_y, y, n_samples);
			return exec_params(y, var, sum_step, n_samples, n_rejected_join, n_rejected_cond);
		};
		exec_params run_round_limit(const size_t &round){
			auto num_rels = data.size();
			if(num_rels == 0){
				std::cerr<<"empty plan"<<std::endl;
				return exec_params(0.0, 0.0, 0, 0, 0, 0);
			}
			double y = 0;
			double sum_y = 0;
			double sum_y2 = 0;
			sum_step = 0;
			double ci = 0;
			n_samples = n_rejected_join = n_rejected_cond = sum_step = 0;
			exec_timer.restart();
			size_t cur_time = 0;
			size_t old_time = 0;
			double decision_value = 0;
			size_t d_step = 0;
			size_t decision_time = 0;
			while(n_samples < round){
				n_samples++;
				auto est = single_run();
				sum_y += est;
				y = sum_y / (n_samples);
				sum_y2 += est * est;
				if(est > 0 && d_step < decision_round){
					d_step++;
					if(d_step == decision_round){
						auto var = tool::calc_var(sum_y2, sum_y, y, n_samples);
						decision_value = var * sum_step;
						decision_time = exec_timer.get_elapsed();
					}
				}
				cur_time = exec_timer.get_elapsed();
				if((cur_time - old_time) > 10000 && cur_time > 0){
					ci = tool::calc_ci(sum_y2, sum_y, y, n_samples, prob);
					tool::report(cur_time, n_samples, n_rejected_join, n_rejected_cond, double(n_rejected_join + n_rejected_cond)/ n_samples, y, ci, prob);
					old_time = cur_time;
				}
			}
			ci = tool::calc_ci(sum_y2, sum_y, y, n_samples, prob);
			auto var = tool::calc_var(sum_y2, sum_y, y, n_samples);
			return exec_params(y, var, sum_step, n_samples, n_rejected_join, n_rejected_cond, decision_value,
					decision_time);
		};
		exec_params run_ci_limit(const double &ci_limit, const walk_plan_params &input){
			auto num_rels = data.size();
			if(num_rels == 0){
				std::cerr<<"empty plan"<<std::endl;
				return exec_params(0.0, 0.0, 0, 0, 0, 0);
			}
			double y = 0;
			double sum_y = input.sum_y;
			double sum_y2 = input.sum_y2;
			sum_step = 0;
			double ci = 0;
			size_t d_step = 0;
			n_rejected_join = n_rejected_cond = 0;
			n_samples = input.n_samples;
			exec_timer.restart();
			while(true){
				n_samples++;
				auto est = single_run();
				sum_y += est;
				y = sum_y / (n_samples);
				sum_y2 += est * est;
				ci = tool::calc_ci(sum_y2, sum_y, y, n_samples, prob);
				if(est > 0 && ci > 0 && ci < ci_limit)
					break;
			}
			ci = tool::calc_ci(sum_y2, sum_y, y, n_samples, prob);
			auto var = tool::calc_var(sum_y2, sum_y, y, n_samples);
			std::cout<<"execution time: " <<exec_timer.get_elapsed()<<std::endl;
			return exec_params(y, var, sum_step, n_samples, n_rejected_join, n_rejected_cond);
		};
		exec_params run_ci_limit(const double &ci_limit){
			walk_plan_params tmp;
			tmp.sum_y = 0;
			tmp.sum_y2 = 0;
			tmp.n_samples = 0;
			return run_ci_limit(ci_limit, tmp);
		};
		walk_plan_params walk_plan_optimizer(){
			sum_step = 0;
			walk_plan_params result;
			double est = 0;
			size_t decision_time = 0;
			double decision_value = 0;
			double y = 0;
			timer mytimer;
			while(result.n_samples < decision_round){
				result.n_samples++;
				est = single_run();
				if(est == 0)
					continue;
				result.sum_y += est;
				y = result.sum_y / result.n_samples;
				result.sum_y2 += est * est;
			}
			double var = tool::calc_var(result.sum_y2, result.sum_y, y, result.n_samples);
			result.decision_value = var * sum_step;
			result.decision_time = mytimer.get_elapsed();
			return result;
		}
		exec_params run_error_limit(const double value, const double &error_limit){
			auto num_rels = data.size();
			if(num_rels == 0){
				std::cerr<<"empty plan"<<std::endl;
				return exec_params(0.0, 0.0, 0, 0, 0, 0);
			}
			double y = 0;
			double sum_y = 0;
			double sum_y2 = 0;
			sum_step = 0;
			n_samples = n_rejected_join = n_rejected_cond = sum_step = 0;
			exec_timer.restart();
			size_t cur_time = 0;
			double decision_value = 0;
			size_t d_step = 0;
			size_t decision_time = 0;
			while(std::abs(y - value) > error_limit){
				n_samples++;
				auto est = single_run();
				sum_y += est;
				y = sum_y / (n_samples);
				sum_y2 += est * est;
				if(est > 0 && d_step < decision_round){
					d_step++;
					if(d_step == decision_round){
						auto var = tool::calc_var(sum_y2, sum_y, y, n_samples);
						decision_value = var * sum_step;
						decision_time = exec_timer.get_elapsed();
					}
				}
				cur_time = exec_timer.get_elapsed();
			}
			auto var = tool::calc_var(sum_y2, sum_y, y, n_samples);
			return exec_params(y, var, cur_time, n_samples, n_rejected_join, n_rejected_cond, decision_value,
					decision_time);
		};
		exec_params run_time_limit(const size_t &max_time){
			auto num_rels = data.size();
			if(num_rels == 0){
				std::cerr<<"empty plan"<<std::endl;
				return exec_params(0.0, 0.0, 0, 0, 0, 0);
			}
			double y = 0;
			double sum_y = 0;
			double sum_y2 = 0;
			sum_step = 0;
			n_samples = n_rejected_join = n_rejected_cond = sum_step = 0;
			exec_timer.restart();
			size_t cur_time = 0;
			size_t old_time = 0;
			double decision_value = 0;
			double ci = 0;
			size_t d_step = 0;
			size_t decision_time = 0;
			while(cur_time < max_time){
				n_samples++;
				auto est = single_run();
				sum_y += est;
				y = sum_y / (n_samples);
				sum_y2 += est * est;
				if(est > 0 && d_step < decision_round){
					d_step++;
					if(d_step == decision_round){
						auto var = tool::calc_var(sum_y2, sum_y, y, n_samples);
						decision_value = var * sum_step;
						decision_time = exec_timer.get_elapsed();
					}
				}
				if((cur_time - old_time) > 10000 && cur_time > 0){
					ci = tool::calc_ci(sum_y2, sum_y, y, n_samples, prob);
					tool::report(cur_time, n_samples, n_rejected_join, n_rejected_cond, double(n_rejected_join + n_rejected_cond)/ n_samples, y, ci, prob);
					old_time = cur_time;
				}
				cur_time = exec_timer.get_elapsed();
			}
			auto var = tool::calc_var(sum_y2, sum_y, y, n_samples);
			return exec_params(y, var, sum_step, n_samples, n_rejected_join, n_rejected_cond, decision_value,
					decision_time);
		};
		exec_params run_time_limit(const size_t &max_time, const size_t &step_time, const double &VALUE){
			auto num_rels = data.size();
			if(num_rels == 0){
				std::cerr<<"empty plan"<<std::endl;
				return exec_params(0.0, 0.0, 0, 0, 0, 0);
			}
			double y = 0;
			double sum_y = 0;
			double sum_y2 = 0;
			sum_step = 0;
			n_samples = n_rejected_join = n_rejected_cond = sum_step = 0;
			exec_timer.restart();
			size_t cur_time = 0;
			size_t old_time = 0;
			double decision_value = 0;
			double ci = 0;
			size_t d_step = 0;
			size_t decision_time = 0;
			while(cur_time < max_time){
				n_samples++;
				auto est = single_run();
				sum_y += est;
				y = sum_y / (n_samples);
				sum_y2 += est * est;
				if(est > 0 && d_step < decision_round){
					d_step++;
					if(d_step == decision_round){
						auto var = tool::calc_var(sum_y2, sum_y, y, n_samples);
						decision_value = var * sum_step;
						decision_time = exec_timer.get_elapsed();
					}
				}
				if((cur_time - old_time) > step_time && cur_time > 0){
					ci = tool::calc_ci(sum_y2, sum_y, y, n_samples, prob);
					if(ci > 0 && ci < VALUE)
						break;
					tool::report(cur_time, n_samples, n_rejected_join, n_rejected_cond, double(n_rejected_join + n_rejected_cond)/ n_samples, y, ci, prob);
					old_time = cur_time;
				}
				cur_time = exec_timer.get_elapsed();
			}
			auto var = tool::calc_var(sum_y2, sum_y, y, n_samples);
			return exec_params(y, var, sum_step, n_samples, n_rejected_join, n_rejected_cond, decision_value,
					decision_time);
		};
		exec_params run_ci_limit_factory(const double &ci_limit, const walk_plan_params &input, const std::function<double(void)> func){
			auto num_rels = data.size();
			if(num_rels == 0){
				std::cerr<<"empty plan"<<std::endl;
				return exec_params(0.0, 0.0, 0, 0, 0, 0);
			}
			double y = 0;
			double sum_y = input.sum_y;
			double sum_y2 = input.sum_y2;
			sum_step = 0;
			double ci = 0;
			size_t d_step = 0;
			n_rejected_join = n_rejected_cond = 0;
			n_samples = input.n_samples;
			exec_timer.restart();
			while(true){
				n_samples++;
				double est = func();
				sum_y += est;
				y = sum_y / (n_samples);
				sum_y2 += est * est;
				ci = tool::calc_ci(sum_y2, sum_y, y, n_samples, prob);
				if(est > 0 && ci > 0 && ci < ci_limit)
					break;
			}
			ci = tool::calc_ci(sum_y2, sum_y, y, n_samples, prob);
			auto var = tool::calc_var(sum_y2, sum_y, y, n_samples);
			std::cout<<"execution time: " <<exec_timer.get_elapsed()<<std::endl;
			return exec_params(y, var, sum_step, n_samples, n_rejected_join, n_rejected_cond);
		};
		exec_params run_time_limit_factory(const size_t &max_time, const size_t &step_time, const
				std::function<double(void)> func){
			auto num_rels = data.size();
			if(num_rels == 0){
				std::cerr<<"empty plan"<<std::endl;
				return exec_params(0.0, 0.0, 0, 0, 0, 0);
			}
			double y = 0;
			double sum_y = 0;
			double sum_y2 = 0;
			sum_step = 0;
			n_samples = n_rejected_join = n_rejected_cond = sum_step = 0;
			exec_timer.restart();
			size_t cur_time = 0;
			size_t old_time = 0;
			double decision_value = 0;
			double ci = 0;
			size_t d_step = 0;
			size_t decision_time = 0;
			size_t cnt = 0;
			while(cur_time < max_time){
				n_samples++;
				double est = func();
				if(est) ++cnt;
				sum_y += est;
				y = sum_y / (n_samples);
				sum_y2 += est * est;
				if(est > 0 && d_step < decision_round){
					d_step++;
					if(d_step == decision_round){
						auto var = tool::calc_var(sum_y2, sum_y, y, n_samples);
						decision_value = var * sum_step;
						decision_time = exec_timer.get_elapsed();
					}
				}
				if((cur_time - old_time) > step_time && cur_time > 0){
					ci = tool::calc_ci(sum_y2, sum_y, y, n_samples, prob);
					tool::report(cur_time, n_samples, n_rejected_join, n_rejected_cond, double(n_rejected_join + n_rejected_cond)/ n_samples, y, ci, prob);
					old_time = cur_time;
				}
				cur_time = exec_timer.get_elapsed();
			}
			auto var = tool::calc_var(sum_y2, sum_y, y, n_samples);
			std::cout<<"# of success: "<<cnt<<std::endl;
			return exec_params(y, var, sum_step, n_samples, n_rejected_join, n_rejected_cond, decision_value,
					decision_time);
		};
		exec_params run_time_limit(const size_t &max_time, const size_t &step_time){
			return run_time_limit_factory(max_time, step_time, std::bind(&plan::single_run, this));
			/*
			auto num_rels = data.size();
			if(num_rels == 0){
				std::cerr<<"empty plan"<<std::endl;
				return exec_params(0.0, 0.0, 0, 0, 0, 0);
			}
			double y = 0;
			double sum_y = 0;
			double sum_y2 = 0;
			sum_step = 0;
			n_samples = n_rejected_join = n_rejected_cond = sum_step = 0;
			exec_timer.restart();
			size_t cur_time = 0;
			size_t old_time = 0;
			double decision_value = 0;
			double ci = 0;
			size_t d_step = 0;
			size_t decision_time = 0;
			size_t cnt = 0;
			while(cur_time < max_time){
				n_samples++;
				auto est = single_run();
				if(est) ++cnt;
				sum_y += est;
				y = sum_y / (n_samples);
				sum_y2 += est * est;
				if(est > 0 && d_step < decision_round){
					d_step++;
					if(d_step == decision_round){
						auto var = tool::calc_var(sum_y2, sum_y, y, n_samples);
						decision_value = var * sum_step;
						decision_time = exec_timer.get_elapsed();
					}
				}
				if((cur_time - old_time) > step_time && cur_time > 0){
					ci = tool::calc_ci(sum_y2, sum_y, y, n_samples, prob);
					tool::report(cur_time, n_samples, n_rejected_join, n_rejected_cond, double(n_rejected_join + n_rejected_cond)/ n_samples, y, ci, prob);
					old_time = cur_time;
				}
				cur_time = exec_timer.get_elapsed();
			}
			auto var = tool::calc_var(sum_y2, sum_y, y, n_samples);
			std::cout<<"# of success: "<<cnt<<std::endl;
			return exec_params(y, var, sum_step, n_samples, n_rejected_join, n_rejected_cond, decision_value,
					decision_time);
					*/
		};
		exec_params run_time_limit_op1(const size_t &max_time, const size_t &step_time){
			return run_time_limit_factory(max_time, step_time, std::bind(&plan::single_run_op1, this));
			/*
			auto num_rels = data.size();
			if(num_rels == 0){
				std::cerr<<"empty plan"<<std::endl;
				return exec_params(0.0, 0.0, 0, 0, 0, 0);
			}
			double y = 0;
			double sum_y = 0;
			double sum_y2 = 0;
			sum_step = 0;
			n_samples = n_rejected_join = n_rejected_cond = sum_step = 0;
			exec_timer.restart();
			size_t cur_time = 0;
			size_t old_time = 0;
			double decision_value = 0;
			double ci = 0;
			size_t d_step = 0;
			size_t decision_time = 0;
			while(cur_time < max_time){
				n_samples++;
				auto est = single_run_op1();
				sum_y += est;
				y = sum_y / (n_samples);
				sum_y2 += est * est;
				if(est > 0 && d_step < decision_round){
					d_step++;
					if(d_step == decision_round){
						auto var = tool::calc_var(sum_y2, sum_y, y, n_samples);
						decision_value = var * sum_step;
						decision_time = exec_timer.get_elapsed();
					}
				}
				if((cur_time - old_time) > step_time && cur_time > 0){
					ci = tool::calc_ci(sum_y2, sum_y, y, n_samples, prob);
					tool::report(cur_time, n_samples, n_rejected_join, n_rejected_cond, double(n_rejected_join + n_rejected_cond)/ n_samples, y, ci, prob);
					old_time = cur_time;
				}
				cur_time = exec_timer.get_elapsed();
			}
			auto var = tool::calc_var(sum_y2, sum_y, y, n_samples);
			return exec_params(y, var, sum_step, n_samples, n_rejected_join, n_rejected_cond, decision_value,
					decision_time);
			*/
		};
		exec_params run_ci_limit_op1(const double &ci_limit){
			auto num_rels = data.size();
			if(num_rels == 0){
				std::cerr<<"empty plan"<<std::endl;
				return exec_params(0.0, 0.0, 0, 0, 0, 0);
			}
			double y = 0;
			double sum_y = 0;
			double sum_y2 = 0;
			sum_step = 0;
			n_samples = n_rejected_join = n_rejected_cond = sum_step = 0;
			exec_timer.restart();
			size_t cur_time = 0;
			size_t old_time = 0;
			double decision_value = 0;
			double ci = 0;
			size_t d_step = 0;
			size_t decision_time = 0;
			while(true){
				n_samples++;
				auto est = single_run_op1();
				sum_y += est;
				y = sum_y / (n_samples);
				sum_y2 += est * est;
				if(est > 0 && d_step < decision_round){
					d_step++;
					if(d_step == decision_round){
						auto var = tool::calc_var(sum_y2, sum_y, y, n_samples);
						decision_time = exec_timer.get_elapsed();
						decision_value = var * decision_time;
					}
				}
				ci = tool::calc_ci(sum_y2, sum_y, y, n_samples, prob);
				if(est > 0 && ci > 0 && ci < ci_limit)
				{
					break;
				}
			}
			auto var = tool::calc_var(sum_y2, sum_y, y, n_samples);
			std::cout<<"execution time: " <<exec_timer.get_elapsed()<<std::endl;
			return exec_params(y, var, sum_step, n_samples, n_rejected_join, n_rejected_cond, decision_value,
					decision_time);
		};
		void first_node_sample_all(plan_node &t, std::vector<plan_result_type> &raw_list, size_t ind){
		base_table *q = t.p_table;
		switch(q->table_tag){
			case EMPTY:
			{
				return;
			}
			case REGION:
			{
				auto p = static_cast<table_region *>(q);
				raw_list[ind].p_raw = p->sample_all();
				raw_list[ind].d = p->data.size();
				return;
			}
			case NATION:
			{
				auto p = static_cast<table_nation *>(q);
				raw_list[ind].p_raw = p->sample_all();
				raw_list[ind].d = p->data.size();
				return;
			}
			case CUSTOMER:
			{
				auto p = static_cast<table_customer *>(q);
				raw_list[ind].p_raw = p->sample_all();
				raw_list[ind].d = p->data.size();
				return;
			}
			case PARTSUPP:
			{
				auto p = static_cast<table_customer *>(q);
				raw_list[ind].p_raw = p->sample_all();
				raw_list[ind].d = p->data.size();
				return;
			}
			case PART:
			{
				auto p = static_cast<table_part *>(q);
				raw_list[ind].p_raw = p->sample_all();
				raw_list[ind].d = p->data.size();
				return;
			}
			case LINEITEM:
			{
				auto p = static_cast<table_lineitem *>(q);
				raw_list[ind].p_raw = p->sample_all();
				raw_list[ind].d = p->data.size();
				return;
			}
			case ORDERS:
			{
				auto p = static_cast<table_orders *>(q);
				raw_list[ind].p_raw = p->sample_all();
				raw_list[ind].d = p->data.size();
				return;
			}
			case SUPPLIER:
			{
				auto p = static_cast<table_supplier *>(q);
				raw_list[ind].p_raw = p->sample_all();
				raw_list[ind].d = p->data.size();
				return;
			}
		}
		}

		exec_params run_time_limit_op2(const size_t &max_time, const size_t &step_time){
			return run_time_limit_factory(max_time, step_time, std::bind(&plan::single_run_op2, this));
			/*
			auto num_rels = data.size();
			if(num_rels == 0){
				std::cerr<<"empty plan"<<std::endl;
				return exec_params(0.0, 0.0, 0, 0, 0, 0);
			}
			double y = 0;
			double sum_y = 0;
			double sum_y2 = 0;
			sum_step = 0;
			n_samples = n_rejected_join = n_rejected_cond = sum_step = 0;
			exec_timer.restart();
			size_t cur_time = 0;
			size_t old_time = 0;
			double decision_value = 0;
			double ci = 0;
			size_t d_step = 0;
			size_t decision_time = 0;
			size_t cnt = 0;
			while(cur_time < max_time){
				n_samples++;
				auto est = single_run_op2();
				if(est) cnt++;
				sum_y += est;
				y = sum_y / (n_samples);
				sum_y2 += est * est;
				if(est > 0 && d_step < decision_round){
					d_step++;
					if(d_step == decision_round){
						auto var = tool::calc_var(sum_y2, sum_y, y, n_samples);
						decision_value = var * sum_step;
						decision_time = exec_timer.get_elapsed();
					}
				}
				if((cur_time - old_time) > step_time && cur_time > 0){
					ci = tool::calc_ci(sum_y2, sum_y, y, n_samples, prob);
					tool::report(cur_time, n_samples, n_rejected_join, n_rejected_cond, double(n_rejected_join + n_rejected_cond)/ n_samples, y, ci, prob);
					old_time = cur_time;
				}
				cur_time = exec_timer.get_elapsed();
			}
			auto var = tool::calc_var(sum_y2, sum_y, y, n_samples);
			std::cout<<"# of success: "<<cnt<<std::endl;
			return exec_params(y, var, sum_step, n_samples, n_rejected_join, n_rejected_cond, decision_value,
					decision_time);
					*/
		};
		exec_params run_ci_limit_op(const double &ci_limit, const walk_plan_params &input){
			return run_ci_limit_factory(ci_limit, input, std::bind(&plan::single_run_op, this));
		}
		exec_params run_time_limit_op(const size_t &max_time, const size_t &step_time){
			return run_time_limit_factory(max_time, step_time, std::bind(&plan::single_run_op, this));
			/*
			auto num_rels = data.size();
			if(num_rels == 0){
				std::cerr<<"empty plan"<<std::endl;
				return exec_params(0.0, 0.0, 0, 0, 0, 0);
			}
			double y = 0;
			double sum_y = 0;
			double sum_y2 = 0;
			sum_step = 0;
			n_samples = n_rejected_join = n_rejected_cond = sum_step = 0;
			exec_timer.restart();
			size_t cur_time = 0;
			size_t old_time = 0;
			double decision_value = 0;
			double ci = 0;
			size_t d_step = 0;
			size_t decision_time = 0;
			size_t cnt = 0;
			while(cur_time < max_time){
				n_samples++;
				auto est = single_run_op();
				if(est) cnt++;
				sum_y += est;
				y = sum_y / (n_samples);
				sum_y2 += est * est;
				if(est > 0 && d_step < decision_round){
					d_step++;
					if(d_step == decision_round){
						auto var = tool::calc_var(sum_y2, sum_y, y, n_samples);
						decision_value = var * sum_step;
						decision_time = exec_timer.get_elapsed();
					}
				}
				if((cur_time - old_time) > step_time && cur_time > 0){
					ci = tool::calc_ci(sum_y2, sum_y, y, n_samples, prob);
					tool::report(cur_time, n_samples, n_rejected_join, n_rejected_cond, double(n_rejected_join + n_rejected_cond)/ n_samples, y, ci, prob);
					old_time = cur_time;
				}
				cur_time = exec_timer.get_elapsed();
			}
			auto var = tool::calc_var(sum_y2, sum_y, y, n_samples);
			std::cout<<"# of success: "<<cnt<<std::endl;
			return exec_params(y, var, sum_step, n_samples, n_rejected_join, n_rejected_cond, decision_value,
					decision_time);
					*/
		};
		void first_node_sample_range(plan_node &t, std::vector<plan_result_type> &raw_list, size_t ind){
		base_table *q = t.p_table;
		switch(q->table_tag){
			case EMPTY:
			{
				return;
			}
			case REGION:
			{
				auto p = static_cast<table_region *>(q);
				raw_list[ind] = p->sample_range(*t.p_head);
				return;
			}
			case NATION:
			{
				auto p = static_cast<table_nation *>(q);
				raw_list[ind] = p->sample_range(*t.p_head);
				return;
			}
			case CUSTOMER:
			{
				auto p = static_cast<table_customer *>(q);
				raw_list[ind] = p->sample_range(*t.p_head);
				return;
			}
			case PARTSUPP:
			{
				auto p = static_cast<table_customer *>(q);
				raw_list[ind] = p->sample_range(*t.p_head);
				return;
			}
			case PART:
			{
				auto p = static_cast<table_part *>(q);
				raw_list[ind] = p->sample_range(*t.p_head);
				return;
			}
			case LINEITEM:
			{
				auto p = static_cast<table_lineitem *>(q);
				raw_list[ind] = p->sample_range(*t.p_head);
				return;
			}
			case ORDERS:
			{
				auto p = static_cast<table_orders *>(q);
				raw_list[ind] = p->sample_range(*t.p_head);
				return;
			}
			case SUPPLIER:
			{
				auto p = static_cast<table_supplier *>(q);
				raw_list[ind] = p->sample_range(*t.p_head);
				return;
			}
		}
		}
		plan_node get_element(size_t i){
			return data[i];
		}
		void set_decision(size_t t){
			decision_round = t;
		}
		void set_agg_id(const size_t &t){
			agg_id = t;
		}
		void set_agg_limit(const size_t &t){
			agg_limit = t;
		}
		void init_op2(){
			size_t n = data.size();
			num_level_sample.assign(n, 0);
			num_level_success.assign(n, 0);
			num_sel.assign(n, 0);
			num_no_op2.assign(n , 0);
		}
		void print_op2(){
			size_t n = data.size();
			for(size_t i = 0; i < n; ++i){
				std::cout<<i<<std::endl;
				std::cout<<"level ratio "<<i<<" : "<<static_cast<double>(num_level_success[i])/num_level_sample[i]<<std::endl;
				std::cout<<"level sample "<<i<<" : "<<num_level_sample[i]<<std::endl;
				std::cout<<"level success "<<i<<" : "<<num_level_success[i]<<std::endl;
				std::cout<<"level sel "<<i<<" : "<<num_sel[i]<<std::endl;
				std::cout<<"no op2 "<<i<<" : "<<num_no_op2[i]<<std::endl;
			}
		}
	private:
		std::vector<plan_node> data;
		std::vector<std::vector<size_t>> adjacency_list;
		std::vector<size_t> order;
		std::vector<size_t> num_level_sample;
		std::vector<size_t> num_level_success;
		std::vector<size_t> num_sel;
		std::vector<size_t> num_no_op2;
		result_func_type result_func;
		size_t n_rejected_join;
		size_t n_rejected_cond;
		size_t n_samples;
		size_t sum_step;
		size_t decision_round;
		size_t agg_id;
		size_t agg_limit;
		timer exec_timer;
		double prob;
};
#endif
