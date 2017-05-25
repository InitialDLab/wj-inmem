#ifndef _H_RIPPLEJOIN
#define _H_RIPPLEJOIN
#include "rawdata.h"
#include "base.h"
#include "timer.h"
#include "tool.h"
#include <vector>
#include <unordered_map>
#include <functional>
#include <climits>
//! template hash join algorithm

//!3 talbe hash join with user defined aspect ratio
template<class R1, class R2, class R3, class K1, class K2, class K3>
void ripple_simple(std::vector<R1> &table1, \
		std::vector<R2> &table2, \
		std::vector<R3> &table3,\
		std::function<K1(const R1 &)> key_func1, \
		std::function<K2(const R2&)> key_func2,\
		std::function<K3(const R3&)>key_func3,\
		std::function<K2(const R1&)> cmp_func1,\
		std::function<K3(const R2&)> cmp_func2,\
		std::function<double(const R1&, const R2&, const R3&)> result_func,
		const double STEP,\
		const size_t &MAX,\
		const double &prob)
{
	constexpr int K = 3;
	size_t max_round = table1.size();;
	max_round = std::min(max_round, table2.size());
	max_round = std::min(max_round, table3.size());
	typedef std::vector<pair<uint32_t, uint32_t>> sample_container;
	sample_container sc1, sc2, sc3;
	sc1.reserve(max_round);
	sc2.reserve(max_round);
	sc3.reserve(max_round);
	std::vector<double> sample_values;
	sample_values.reserve(max_round);

	tool::start_report();
	size_t test_count = 0;
	typedef unordered_multimap<uint32_t, size_t> index_type;
	index_type l_index_table1_key2;
	index_type l_index_table2_key2;
	index_type l_index_table2_key3;
	index_type l_index_table3_key2;
	std::vector<double> q[K]; 
	std::vector<double> w(3, 0);
	std::vector<double> x(3, 0);
	for (int k = 0; k < K; ++k) q[k].resize(max_round);

	const double total_size = (double) table1.size() * table2.size() * table3.size();
	cout<<total_size<<endl;
	const double zp = tool::xql_erf_inv(prob);

	double y = 0.0;
	double ci = 1.0 / 0.0;
	size_t round = 0;
	/* round 1 */
	sc1.push_back(std::make_pair(cmp_func1(table1[round]), round)); 
	l_index_table1_key2.insert(std::make_pair(cmp_func1(table1[round]), round)); 
	sc2.push_back(std::make_pair(key_func2(table2[round]), round)); 
	l_index_table2_key2.insert(std::make_pair(key_func2(table2[round]), sc2.size() - 1));
	l_index_table2_key3.insert(std::make_pair(cmp_func2(table2[round]), sc2.size() - 1));
	sc3.push_back(std::make_pair(key_func3(table3[round]), round));
	l_index_table3_key2.insert(std::make_pair(key_func3(table3[round]), sc3.size() - 1));
	sample_values.push_back(result_func(table1[round], table2[round], table3[round]));	
	auto e1 = table1[sc1.back().second];
	auto e2 = table2[sc2.back().second];
	auto e3 = table3[sc3.back().second];
	double sum = 0;
	if(cmp_func1(e1) == key_func2(e2) && cmp_func2(e2) == key_func3(e3)){
		double v = sample_values.back();
		y = v *total_size;
		sum+=v;
		test_count++;
		for (int k = 0; k < K; ++k) {
			q[k][0] = v;
			w[k] = 0;
			x[k] = y;
		}
	}
	size_t n = 1;
	double z = 0;
	uint32_t n_rejected = 0;
	size_t report_round = 1;
	timer maxtimer;
	timer steptimer;
	while(maxtimer.get_elapsed() < MAX && n < max_round){
		n++;
		sc1.push_back(std::make_pair(cmp_func1(table1[n-1]), n-1)); 
		l_index_table1_key2.insert(std::make_pair(cmp_func1(table1[n-1]), n-1)); 
		sc2.push_back(std::make_pair(key_func2(table2[n-1]), n-1)); 
		l_index_table2_key2.insert(std::make_pair(key_func2(table2[n-1]), sc2.size() - 1));
		l_index_table2_key3.insert(std::make_pair(cmp_func2(table2[n-1]), sc2.size() - 1));
		sc3.push_back(std::make_pair(key_func3(table3[n-1]), n-1));
		l_index_table3_key2.insert(std::make_pair(key_func3(table3[n-1]), sc3.size() - 1));
		sample_values.push_back(result_func(table1[n-1], table2[n-1], table3[n-1]));	
		y = y * std::pow(double(n-1)/n, 2);
		auto e1 = table1[sc1.back().second];
		auto e2 = table2[sc2.back().second];
		auto e3 = table3[sc3.back().second];
		if(cmp_func1(e1) == key_func2(e2) && cmp_func2(e2) == key_func3(e3))
		{
			sum += sample_values.back();
			test_count++;
		}

		for(size_t k = 0; k < K; ++k){
			q[k][n-1] = 0;
			w[k] = w[k] * std::pow(double(n-1) / n, 4) + y * y /n/(n-1);
			x[k] = y;
		}
		//ik == n from table 1
		{
		std::vector<size_t> ind(3);
		ind[0] = sc1[n-1].second;
		auto range2 = l_index_table2_key2.equal_range(sc1[n-1].first);
		auto tuple1 = table1[ind[0]];
		for(auto iter = range2.first; iter != range2.second; iter++){
			ind[1] = iter->second;
			if(ind[1] >= n - 1)
				continue;
			auto tuple2 = table2[ind[1]];
			auto range3 = l_index_table3_key2.equal_range(cmp_func2(tuple2));
			for(auto iter3 = range3.first; iter3 != range3.second; ++iter3){
				ind[2] = iter3->second;
				if(ind[2] >= n-1)
					continue;
				auto tuple3 = table3[ind[2]];
				for(size_t k = 0; k < K; ++k){
					x[k] = x[k] - q[k][ind[k]]/n/n;
					w[k] = w[k] - (x[k] - double(n-1)/n/n*q[k][ind[k]])/n*(x[k]-double(n-1)/n/n*q[k][ind[k]])/(n-1);
					q[k][ind[k]] = q[k][ind[k]] + sample_values[ind[2]]*total_size;
					w[k] = w[k] + (x[k] - double(n-1)/n/n*q[k][ind[k]])/n*(x[k]-double(n-1)/n/n*q[k][ind[k]])/(n-1);
					x[k] = x[k] + q[k][ind[k]]/n/n;
				}
				sum += sample_values[ind[2]];
			}
		}
		}
		//ik == n from table2
		{
			std::vector<size_t> ind(3);
			ind[1] = sc2.back().second;
			auto tuple2 = table2[ind[1]];
			auto range1 = l_index_table1_key2.equal_range(key_func2(tuple2));
			for(auto iter1 = range1.first; iter1 != range1.second; ++iter1){
				ind[0] = iter1->second;
				if(ind[0] >= n-1)
					continue;
				auto tuple1 = table1[ind[0]];
				auto range3 = l_index_table3_key2.equal_range(cmp_func2(tuple2));
				for(auto iter3 = range3.first; iter3 != range3.second; ++iter3){
					ind[2] = iter3->second;
					if(ind[2] >= n-1)
						continue;
					auto tuple3 = table3[ind[2]];
					for(size_t k = 0; k < K; ++k){
						x[k] = x[k] - q[k][ind[k]]/n/n;
						w[k] = w[k] - (x[k] - double(n-1)/n/n*q[k][ind[k]])/n*(x[k]-double(n-1)/n/n*q[k][ind[k]])/(n-1);
						q[k][ind[k]] = q[k][ind[k]] + sample_values[ind[2]] * total_size;
						w[k] = w[k] + (x[k] - double(n-1)/n/n*q[k][ind[k]])/n*(x[k]-double(n-1)/n/n*q[k][ind[k]])/(n-1);
						x[k] = x[k] + q[k][ind[k]]/n/n;
					}
					sum += sample_values[ind[2]];
				}
			}
		}
		{
			std::vector<size_t> ind(3);
			ind[2] = sc3.back().second;
			auto tuple3 = table3[ind[2]];
			auto range2 = l_index_table2_key3.equal_range(key_func3(tuple3));
			for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
				ind[1] = iter2->second;
				if(ind[1] >= n-1)
					continue;
				auto tuple2 = table2[ind[1]];
				auto range1 = l_index_table1_key2.equal_range(key_func2(tuple2));
				for(auto iter1 = range1.first; iter1 != range1.second; ++iter1){
					ind[0] = iter1->second;
					if(ind[0] >= n-1)
						continue;
					auto tuple1 = table1[ind[0]];
					for(size_t k = 0; k < K; ++k){
						x[k] = x[k] - q[k][ind[k]]/n/n;
						w[k] = w[k] - (x[k] - double(n-1)/n/n*q[k][ind[k]])/n*(x[k]-double(n-1)/n/n*q[k][ind[k]])/(n-1);
						q[k][ind[k]] = q[k][ind[k]] + sample_values[ind[2]] * total_size;
						w[k] = w[k] + (x[k] - double(n-1)/n/n*q[k][ind[k]])/n*(x[k]-double(n-1)/n/n*q[k][ind[k]])/(n-1);
						x[k] = x[k] + q[k][ind[k]]/n/n;
					}
					sum += sample_values[ind[2]];
				}
			}
		}
		y = x[0];
		z = 0;
		double nk=std::pow(n, K);
		for(size_t i = 0; i < K; ++i)
			z += w[i] / (n-1);
		if(steptimer.get_elapsed()> STEP){
			ci = zp*std::sqrt(z) / sqrt(n);
			tool::report(STEP*report_round, n, 0, 0, 0, sum*total_size/nk, ci, prob);
			steptimer.restart();
			report_round++;
		}
	}
}
//with conditions
template<class R1, class R2, class R3, class K1, class K2, class K3>
void ripple_simple(std::vector<R1> &table1, \
		std::vector<R2> &table2, \
		std::vector<R3> &table3,\
		std::function<K1(const R1 &)> key_func1, \
		std::function<K2(const R2&)> key_func2,\
		std::function<K3(const R3&)>key_func3,\
		std::function<K2(const R1&)> cmp_func1,\
		std::function<K3(const R2&)> cmp_func2,\
		std::function<bool(const R1&)> cond1,\
		std::function<bool(const R2&)> cond2,\
		std::function<bool(const R3&)> cond3,\
		std::function<double(const R1&, const R2&, const R3&)> result_func,
		const double STEP,\
		const size_t &MAX,\
		const double &prob)
{
	constexpr int K = 3;
	size_t max_round = table1.size();;
	max_round = std::min(max_round, table2.size());
	max_round = std::min(max_round, table3.size());
	typedef std::vector<pair<uint32_t, uint32_t>> sample_container;
	sample_container sc1, sc2, sc3;
	sc1.reserve(max_round);
	sc2.reserve(max_round);
	sc3.reserve(max_round);
	std::vector<double> sample_values;
	sample_values.reserve(max_round);

	tool::start_report();
	typedef unordered_multimap<uint32_t, size_t> index_type;
	index_type l_index_table1_key2;
	index_type l_index_table2_key2;
	index_type l_index_table2_key3;
	index_type l_index_table3_key2;
	std::vector<double> q[K]; 
	std::vector<double> w(3, 0);
	std::vector<double> x(3, 0);
	for (int k = 0; k < K; ++k) q[k].resize(max_round);

	const double total_size = 
			(double) table1.size() * 
			table2.size() * 
			table3.size();
	const double zp = tool::xql_erf_inv(prob);

	double y = 0.0;
	double ci = 1.0 / 0.0;
	size_t round = 0;
	/* round 1 */
	sc1.push_back(std::make_pair(cmp_func1(table1[round]), round)); 
	l_index_table1_key2.insert(std::make_pair(cmp_func1(table1[round]), round)); 
	sc2.push_back(std::make_pair(key_func2(table2[round]), round)); 
	l_index_table2_key2.insert(std::make_pair(key_func2(table2[round]), sc2.size() - 1));
	l_index_table2_key3.insert(std::make_pair(cmp_func2(table2[round]), sc2.size() - 1));
	sc3.push_back(std::make_pair(key_func3(table3[round]), round));
	l_index_table3_key2.insert(std::make_pair(key_func3(table3[round]), sc3.size() - 1));
	sample_values.push_back(result_func(table1[round], table2[round], table3[round]));	
	auto e1 = table1[sc1.back().second];
	auto e2 = table2[sc2.back().second];
	auto e3 = table3[sc3.back().second];
	if(cmp_func1(e1) == key_func2(e2) && cmp_func2(e2) == key_func3(e3) && cond1(e1) && cond2(e2) && cond3(e3)){
		double v = sample_values.back();
		y = v;
		for (int k = 0; k < K; ++k) {
			q[k][0] = v;
			w[k] = 0;
			x[k] = y;
		}
	}
	size_t n = 1;
	double z = 0;
	size_t n_rejected_join = 0;
	size_t n_rejected_cond = 0;
	
	size_t report_round = 1;
	timer maxtimer;
	timer steptimer;
	while(maxtimer.get_elapsed() < MAX && n < max_round){
		n++;
		sc1.push_back(std::make_pair(cmp_func1(table1[n-1]), n-1)); 
		l_index_table1_key2.insert(std::make_pair(cmp_func1(table1[n-1]), n-1)); 
		sc2.push_back(std::make_pair(key_func2(table2[n-1]), n-1)); 
		l_index_table2_key2.insert(std::make_pair(key_func2(table2[n-1]), sc2.size() - 1));
		l_index_table2_key3.insert(std::make_pair(cmp_func2(table2[n-1]), sc2.size() - 1));
		sc3.push_back(std::make_pair(key_func3(table3[n-1]), n-1));
		l_index_table3_key2.insert(std::make_pair(key_func3(table3[n-1]), sc3.size() - 1));
		sample_values.push_back(result_func(table1[n-1], table2[n-1], table3[n-1]));	
		y = y * std::pow(double(n-1)/n, 2);
		for(size_t k = 0; k < K; ++k){
			q[k][n-1] = 0;
			w[k] = w[k] * std::pow(double(n-1) / n, 4) + y * y /n/(n-1);
			x[k] = y;
		}
		//ik == n from table 1
		{
		std::vector<size_t> ind(3);
		ind[0] = sc1[n-1].second;
		auto range2 = l_index_table2_key2.equal_range(sc1[n-1].first);
		auto tuple1 = table1[ind[0]];
		if(!cond1(tuple1)){
			n_rejected_cond++;
			continue;
		}
		for(auto iter = range2.first; iter != range2.second; iter++){
			ind[1] = iter->second;
			if(ind[1] == n - 1){
				continue;
			}
			auto tuple2 = table2[ind[1]];
			if(!cond2(tuple2)){
				n_rejected_cond++;
				continue;
			}
			auto range3 = l_index_table3_key2.equal_range(cmp_func2(tuple2));
			for(auto iter3 = range3.first; iter3 != range3.second; ++iter3){
				ind[2] = iter3->second;
				if(ind[2] == n-1){
					continue;
				}
				auto tuple3 = table3[ind[2]];
				if(!cond3(tuple3))
				{
					n_rejected_cond++;
					continue;
				}
				for(size_t k = 0; k < K; ++k){
					x[k] = x[k] - q[k][ind[k]]/n/n;
					w[k] = w[k] - (x[k] - double(n-1)/n/n*q[k][ind[k]])/n*(x[k]-double(n-1)/n/n*q[k][ind[k]])/(n-1);
					q[k][ind[k]] = q[k][ind[k]] + sample_values[ind[2]];
					w[k] = w[k] + (x[k] - double(n-1)/n/n*q[k][ind[k]])/n*(x[k]-double(n-1)/n/n*q[k][ind[k]])/(n-1);
					x[k] = x[k] + q[k][ind[k]]/n/n;
				}
			}
		}
		}
		//ik == n from table2
		{
			std::vector<size_t> ind(3);
			ind[1] = sc2.back().second;
			auto tuple2 = table2[ind[1]];
			if(!cond2(tuple2)){
				n_rejected_cond++;
				continue;
			}
			auto range1 = l_index_table1_key2.equal_range(key_func2(tuple2));
			for(auto iter1 = range1.first; iter1 != range1.second; ++iter1){
				ind[0] = iter1->second;
				if(ind[0] == n-1)
				{
					continue;
				}
				auto tuple1 = table1[ind[0]];
				if(!cond1(tuple1)){
					n_rejected_cond++;
					continue;
				}
				auto range3 = l_index_table3_key2.equal_range(cmp_func2(tuple2));
				for(auto iter3 = range3.first; iter3 != range3.second; ++iter3){
					ind[2] = iter3->second;
					if(ind[2] == n-1){
						continue;
					}
					auto tuple3 = table3[ind[2]];
					if(!cond3(tuple3)){
						++n_rejected_cond;
						continue;
					}
					for(size_t k = 0; k < K; ++k){
						x[k] = x[k] - q[k][ind[k]]/n/n;
						w[k] = w[k] - (x[k] - double(n-1)/n/n*q[k][ind[k]])/n*(x[k]-double(n-1)/n/n*q[k][ind[k]])/(n-1);
						q[k][ind[k]] = q[k][ind[k]] + sample_values[ind[2]];
						w[k] = w[k] + (x[k] - double(n-1)/n/n*q[k][ind[k]])/n*(x[k]-double(n-1)/n/n*q[k][ind[k]])/(n-1);
						x[k] = x[k] + q[k][ind[k]]/n/n;
					}
				}
			}
		}
		{
			std::vector<size_t> ind(3);
			ind[2] = sc3.back().second;
			auto tuple3 = table3[ind[2]];
			if(!cond3(tuple3)){
				++n_rejected_cond;
				continue;
			}
			auto range2 = l_index_table2_key3.equal_range(key_func3(tuple3));
			for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
				ind[1] = iter2->second;
				if(ind[1] == n-1){
					continue;
				}
				auto tuple2 = table2[ind[1]];
				if(!cond2(tuple2)){
					++n_rejected_cond;
					continue;
				}
				auto range1 = l_index_table1_key2.equal_range(key_func2(tuple2));
				for(auto iter1 = range1.first; iter1 != range1.second; ++iter1){
					ind[0] = iter1->second;
					if(ind[0] == n-1){
						continue;
					}
					auto tuple1 = table1[ind[0]];
					if(!cond1(tuple1)){
						++n_rejected_cond;
						continue;
					}
					for(size_t k = 0; k < K; ++k){
						x[k] = x[k] - q[k][ind[k]]/n/n;
						w[k] = w[k] - (x[k] - double(n-1)/n/n*q[k][ind[k]])/n*(x[k]-double(n-1)/n/n*q[k][ind[k]])/(n-1);
						q[k][ind[k]] = q[k][ind[k]] + sample_values[ind[2]];
						w[k] = w[k] + (x[k] - double(n-1)/n/n*q[k][ind[k]])/n*(x[k]-double(n-1)/n/n*q[k][ind[k]])/(n-1);
						x[k] = x[k] + q[k][ind[k]]/n/n;
					}
				}
			}
		}
		y = x[0];
		z = 0;
		for(size_t i = 0; i < K; ++i)
			z += w[i] / (n-1);
		if(steptimer.get_elapsed()> STEP){
			ci = zp*std::sqrt(z) / sqrt(n);
			double tmp = y*table1.size() *table2.size() * table3.size() /n /n /n ;
			tool::report(STEP*report_round, n, n_rejected_join, n_rejected_cond, double(n_rejected_join +
						n_rejected_cond) /n /n /n, y/n, ci, prob);
			steptimer.restart();
			report_round++;
		}
	}
}
#endif
