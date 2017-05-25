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
#include <random>
#include <chrono>
//calc ci for ripple join
void calc_ripple_ci(const size_t &K, std::vector<size_t> &ind, std::vector<double> &x, std::vector<double> &w, std::vector<double> *q,
		std::vector<double> &sample_values, size_t &n, std::vector<size_t> &a, size_t total_size){
		std::vector<double> pow_n(K);
		double all = 1;
		for(size_t i = 0; i < K; ++i)
			all *= n *a[i];
		for(size_t i = 0; i < K; ++i)
			pow_n[i] = all / n/ a[i];
		for(size_t k = 0; k < K; ++k){
			x[k] = x[k] - q[k][ind[k]/a[k]]/pow_n[k];
			w[k] = w[k] - (x[k] - double(n*a[k]-1)/pow_n[k]*q[k][ind[k]/a[k]])/n/a[k]*(x[k]-double(n*a[k]-1)/pow_n[k]*q[k][ind[k]/a[k]])/(n*a[k]-1);
			q[k][ind[k]/a[k]] = q[k][ind[k]/a[k]] + sample_values[ind[2]]*total_size/a[k];
			w[k] = w[k] + (x[k] - double(n*a[k]-1)/pow_n[k]*q[k][ind[k]/a[k]])/n/a[k]*(x[k]-double(n-1)/pow_n[k]*q[k][ind[k]/a[k]])/(n*a[k]-1);
			x[k] = x[k] + q[k][ind[k]/a[k]]/pow_n[k];
		}
	/*
		for(size_t k = 0; k < K; ++k){
			x[k] = x[k] - q[k][ind[k]]/std::pow(n, K-1);
			w[k] = w[k] - (x[k] - double(n-1)/std::pow(n, K-1)*q[k][ind[k]])/n*(x[k]-double(n-1)/std::pow(n, K-1)*q[k][ind[k]])/(n-1);
			q[k][ind[k]] = q[k][ind[k]] + sample_values[ind[2]]*total_size;
			w[k] = w[k] + (x[k] - double(n-1)/std::pow(n, K-1)*q[k][ind[k]])/n*(x[k]-double(n-1)/std::pow(n, K-1)*q[k][ind[k]])/(n-1);
			x[k] = x[k] + q[k][ind[k]]/std::pow(n, K-1);
		}
		*/
}
void calc_ripple_ci_test(const size_t &K, std::vector<size_t> &ind, std::vector<double> &x, std::vector<double> &w, std::vector<double> *q,
		std::vector<double> &sample_values, size_t &n, std::vector<size_t> &a, size_t total_size){
		std::vector<double> pow_n(K);
		double all = 1;
		for(size_t i = 0; i < K; ++i)
			all *= n *a[i];
		for(size_t i = 0; i < K; ++i)
			pow_n[i] = all / n/ a[i];
		for(size_t k = 0; k < K; ++k){
			x[k] = x[k] - q[k][ind[k]]/pow_n[k]/a[k];
			w[k] = w[k] - (x[k] -
					double(n*a[k]-1)/pow_n[k]/a[k]*q[k][ind[k]])/n/a[k]*(x[k]-double(n*a[k]-1)/pow_n[k]/a[k]*q[k][ind[k]])/(n*a[k]-1);
			q[k][ind[k]] = q[k][ind[k]] + sample_values[ind[2]]*total_size;
			w[k] = w[k] + (x[k] - double(n*a[k]-1)/pow_n[k]/a[k]*q[k][ind[k]])/n/a[k]*(x[k]-double(n-1)/pow_n[k]/a[k]*q[k][ind[k]])/(n*a[k]-1);
			x[k] = x[k] + q[k][ind[k]]/pow_n[k]/a[k];
		}
}
void calc_ripple_ci_test(const size_t &K, std::vector<size_t> &ind, std::vector<double> &x, std::vector<double> &w, std::vector<double> *q,
		std::vector<double> &sample_values, size_t &n, std::vector<size_t> &a, size_t total_size, size_t col_id){
		std::vector<double> pow_n(K);
		double all = 1;
		for(size_t i = 0; i < K; ++i)
			all *= n *a[i];
		for(size_t i = 0; i < K; ++i)
			pow_n[i] = all / n/ a[i];
		for(size_t k = 0; k < K; ++k){
			x[k] = x[k] - q[k][ind[k]]/pow_n[k]/a[k];
			w[k] = w[k] - (x[k] - double(n*a[k]-1)/pow_n[k]/a[k]*q[k][ind[k]])/n/a[k]*(x[k]-double(n*a[k]-1)/pow_n[k]/a[k]*q[k][ind[k]])/(n*a[k]-1);
			q[k][ind[k]] = q[k][ind[k]] + sample_values[ind[col_id]]*total_size;
			w[k] = w[k] + (x[k] - double(n*a[k]-1)/pow_n[k]/a[k]*q[k][ind[k]])/n/a[k]*(x[k]-double(n*a[k]-1)/pow_n[k]/a[k]*q[k][ind[k]])/(n*a[k]-1);
			q[k][ind[k]] = q[k][ind[k]] + sample_values[ind[col_id]]*total_size;
			x[k] = x[k] + q[k][ind[k]]/pow_n[k]/a[k];
		}
}
void calc_ripple_ci_test(const size_t &K, std::vector<size_t> &ind, std::vector<double> &x, std::vector<double> &w, std::vector<double> *q,
		std::vector<double> &sample_values, size_t &n, std::vector<size_t> &a, double &total_size, size_t col_id){
		std::vector<double> pow_n(K);
		double all = 1;
		for(size_t i = 0; i < K; ++i)
			all *= n *a[i];
		for(size_t i = 0; i < K; ++i)
			pow_n[i] = all / n/ a[i];
		for(size_t k = 0; k < K; ++k){
			x[k] = x[k] - q[k][ind[k]]/pow_n[k]/a[k];
			w[k] = w[k] - (x[k] - double(n*a[k]-1)/pow_n[k]/a[k]*q[k][ind[k]])/n/a[k]*(x[k]-double(n*a[k]-1)/pow_n[k]/a[k]*q[k][ind[k]])/(n*a[k]-1);
			q[k][ind[k]] = q[k][ind[k]] + sample_values[ind[col_id]]*total_size;
			w[k] = w[k] + (x[k] - double(n*a[k]-1)/pow_n[k]/a[k]*q[k][ind[k]])/n/a[k]*(x[k]-double(n*a[k]-1)/pow_n[k]/a[k]*q[k][ind[k]])/(n*a[k]-1);
			q[k][ind[k]] = q[k][ind[k]] + sample_values[ind[col_id]]*total_size;
			x[k] = x[k] + q[k][ind[k]]/pow_n[k]/a[k];
		}
}
//! template hash join algorithm

//!3 table hash join with user defined aspect ratio
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
	for (size_t k = 0; k < K; ++k) q[k].resize(max_round);

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
		for (size_t k = 0; k < K; ++k) {
			q[k][0] = v;
			w[k] = 0;
			x[k] = y;
		}
	}
	size_t n = 1;
	double z = 0;
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
		const double &prob,\
		bool verbose)
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
	for (size_t k = 0; k < K; ++k) q[k].resize(max_round);

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
		for (size_t k = 0; k < K; ++k) {
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
		if(verbose && steptimer.get_elapsed()> STEP){
			ci = zp*std::sqrt(z) / sqrt(n);
			double tmp = y*table1.size() *table2.size() * table3.size() /n /n /n ;
			tool::report(STEP*report_round, n, n_rejected_join, n_rejected_cond, double(n_rejected_join +
						n_rejected_cond) /n /n /n, y/n, ci, prob);
			steptimer.restart();
			report_round++;
		}
	}
}
void ripple_join(std::vector<std::vector<base_raw *> > &input, cond_func_type
cond_func){
	std::cout<<"doing ripple join"<<std::endl;
}
//pointers with conditions

//! key func definition key_func_table#_query table
//ripple join for query3
void ripple_3(std::vector<std::vector<base_raw *> > &table,\
		std::vector<std::vector<key_func_type>> &key_func,\
		std::vector<cond_func_type> &cond_func,\
		std::function<double(base_raw *, base_raw *, base_raw*)> result_func,
		const double STEP,\
		const size_t &MAX,\
		const double &prob,\
		bool verbose,\
		double total_size){
	constexpr int K = 3;
	size_t max_round = 0;
	for(size_t i = 0; i < K; ++i){
		if(i == 0)
			max_round = table[i].size();
		else
			max_round = std::min(max_round, table[i].size());
	}
	std::vector<double> sample_values;

	tool::start_report();
	size_t test_count = 0;
	typedef unordered_multimap<uint32_t, size_t> index_type;
	std::vector<std::vector<index_type> > h_table;
	for(size_t i = 0; i < K; ++i){
		std::vector<index_type> tmp;
		h_table.push_back(tmp);
		for(size_t j = 0; j < K; ++j){
			index_type x;
			h_table[i].push_back(x);
		}
	}
	std::vector<double> q[K]; 
	std::vector<double> w(K, 0);
	std::vector<double> x(K, 0);
	for (size_t k = 0; k < K; ++k) q[k].resize(max_round);

	cout<<total_size<<endl;
	const double zp = tool::xql_erf_inv(prob);

	double y = 0.0;
	double ci = 1.0 / 0.0;
	size_t round = 0;
	/* round 1 */
	for(size_t i = 0; i < K; ++i){
		for(size_t j = 0; j < K; ++j){
			if(key_func[i][j] != nullptr)
				if(cond_func[i] == nullptr || (*cond_func[i])(table[i][round]))
					h_table[i][j].insert(std::make_pair((*key_func[i][j])(table[i][round]), round)); 
		}
	}
	double sum = 0;
	double v = result_func(table[0][round], table[1][round], table[2][round]);
	if(v > 0)
		test_count++;
	sample_values.push_back(v);
	if((*key_func[0][1])(table[0][round]) == (*key_func[1][0])(table[1][round]) && (*key_func[1][2])(table[1][round]) == (*key_func[2][1])(table[2][round])){
			if((cond_func[0] == nullptr || (*cond_func[0])(table[0][round])) &&\
					(cond_func[1] == nullptr || (*cond_func[1])(table[1][round])) &&\
					(cond_func[2] == nullptr || (*cond_func[2])(table[2][round]))){
				y = v *total_size;
				sum+=v;
				test_count++;
			}
	}
	for (size_t k = 0; k < K; ++k) {
		q[k][0] = v;
		w[k] = 0;
		x[k] = y;
	}
	size_t n = 1;
	double z = 0;
	size_t report_round = 1;
	timer maxtimer;
	timer steptimer;
	round = 1;
	while(maxtimer.get_elapsed() < MAX && round < max_round){
		n = round + 1;
		for(size_t i = 0; i < K; ++i){
			for(size_t j = 0; j < K; ++j){
				if(key_func[i][j] != nullptr)
					if(cond_func[i] == nullptr || (*cond_func[i])(table[i][round]))
						h_table[i][j].insert(std::make_pair((*key_func[i][j])(table[i][round]), round)); 
			}
		}
		sample_values.push_back(result_func(table[0][round], table[1][round], table[2][round]));	
		y = y * std::pow(double(round)/(round+1), K-1);
		v = sample_values.back();
		if((*key_func[0][1])(table[0][round]) == (*key_func[1][0])(table[1][round]) && (*key_func[1][2])(table[1][round]) == (*key_func[2][1])(table[2][round])){
			if((cond_func[0] == nullptr || (*cond_func[0])(table[0][round])) &&\
					(cond_func[1] == nullptr || (*cond_func[1])(table[1][round]))&&\
					(cond_func[2] == nullptr || (*cond_func[2])(table[2][round]))){
				y = v *total_size;
				sum+=v;
				test_count++;
			}
		}

		for(size_t k = 0; k < K; ++k){
			q[k][round] = 0;
			w[k] = w[k] * std::pow(double(round) / (round+1), 4) + y * y /(round+1)/round;
			x[k] = y;
		}
		//ik == n from table 0
		{
		std::vector<size_t> ind(K);
		ind[0] = round;
		auto tuple0 = table[0][ind[0]];
		if(cond_func[0] == nullptr || cond_func[0] != nullptr && (*cond_func[0])(tuple0)){
			auto range1 = h_table[1][0].equal_range((*key_func[0][1])(tuple0));
			for(auto iter1 = range1.first; iter1 != range1.second; iter1++){
				ind[1] = iter1->second;
				if(ind[1] >= round)
					continue;
				auto tuple1 = table[1][ind[1]];
				if(cond_func[1] != nullptr && (*cond_func[1])(tuple1) == false)
					continue;
				auto range2 = h_table[2][1].equal_range((*key_func[1][2])(tuple1));
				for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
					ind[2] = iter2->second;
					if(ind[2] >= round)
						continue;
					auto tuple2 = table[2][ind[2]];
					if(cond_func[2] != nullptr && (*cond_func[2])(tuple2) == false)
						continue;
					for(size_t k = 0; k < K; ++k){
						x[k] = x[k] - q[k][ind[k]]/n/n;
						w[k] = w[k] - (x[k] - double(n-1)/n/n*q[k][ind[k]])/n*(x[k]-double(n-1)/n/n*q[k][ind[k]])/(n-1);
						q[k][ind[k]] = q[k][ind[k]] + sample_values[ind[2]]*total_size;
						w[k] = w[k] + (x[k] - double(n-1)/n/n*q[k][ind[k]])/n*(x[k]-double(n-1)/n/n*q[k][ind[k]])/(n-1);
						x[k] = x[k] + q[k][ind[k]]/n/n;
					}
					sum += sample_values[ind[2]];
					test_count++;
				}
			}
			}
		}
		//ik == n from table 1
		{
		std::vector<size_t> ind(K);
		ind[1] = round;
		auto tuple1 = table[1][ind[1]];
		if(cond_func[1] == nullptr || (*cond_func[1])(tuple1)){
			auto range0 = h_table[0][1].equal_range((*key_func[1][0])(tuple1));
			for(auto iter0 = range0.first; iter0 != range0.second; iter0++){
				ind[0] = iter0->second;
				if(ind[0] >= round)
					continue;
				auto tuple0 = table[0][ind[0]];
				if(cond_func[0] != nullptr && !(*cond_func[0])(tuple0))
					continue;
				auto range2 = h_table[2][1].equal_range((*key_func[1][2])(tuple1));
				for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
					ind[2] = iter2->second;
					if(ind[2] >= round)
						continue;
					auto tuple2 = table[2][ind[2]];
					if(cond_func[2] != nullptr && !(*cond_func[2])(tuple2))
						continue;
					for(size_t k = 0; k < K; ++k){
						x[k] = x[k] - q[k][ind[k]]/n/n;
						w[k] = w[k] - (x[k] - double(n-1)/n/n*q[k][ind[k]])/n*(x[k]-double(n-1)/n/n*q[k][ind[k]])/(n-1);
						q[k][ind[k]] = q[k][ind[k]] + sample_values[ind[2]]*total_size;
						w[k] = w[k] + (x[k] - double(n-1)/n/n*q[k][ind[k]])/n*(x[k]-double(n-1)/n/n*q[k][ind[k]])/(n-1);
						x[k] = x[k] + q[k][ind[k]]/n/n;
					}
					sum += sample_values[ind[2]];
					test_count++;
				}
			}
		}
		}
		//ik == n from table 2
		{

		std::vector<size_t> ind(K);
		ind[2] = round;
		auto tuple2 = table[2][ind[2]];
		if(cond_func[2] == nullptr || (*cond_func[2])(tuple2)){
			auto range1 = h_table[1][2].equal_range((*key_func[2][1])(tuple2));
			for(auto iter1 = range1.first; iter1 != range1.second; iter1++){
				ind[1] = iter1->second;
				if(ind[1] >= round)
					continue;
				auto tuple1 = table[1][ind[1]];
				if(cond_func[1] != nullptr && !(*cond_func[1])(tuple1))
					continue;
				auto range0 = h_table[0][1].equal_range((*key_func[1][0])(tuple1));
				for(auto iter0 = range0.first; iter0 != range0.second; ++iter0){
					ind[0] = iter0->second;
					if(ind[0] >= round)
						continue;
					auto tuple0 = table[0][ind[0]];
					if(cond_func[0] != nullptr && !(*cond_func[0])(tuple0))
						continue;
					for(size_t k = 0; k < K; ++k){
						x[k] = x[k] - q[k][ind[k]]/n/n;
						w[k] = w[k] - (x[k] - double(n-1)/n/n*q[k][ind[k]])/n*(x[k]-double(n-1)/n/n*q[k][ind[k]])/(n-1);
						q[k][ind[k]] = q[k][ind[k]] + sample_values[ind[2]]*total_size;
						w[k] = w[k] + (x[k] - double(n-1)/n/n*q[k][ind[k]])/n*(x[k]-double(n-1)/n/n*q[k][ind[k]])/(n-1);
						x[k] = x[k] + q[k][ind[k]]/n/n;
					}
					sum += sample_values[ind[2]];
					test_count++;
				}
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
			tool::report(STEP*report_round, n, test_count, 0, 0, sum*total_size/nk, ci, prob);
			steptimer.restart();
			report_round++;
		}
		round++;
	}
}
void ripple_3_ratio(std::vector<std::vector<base_raw *> > &table,\
		std::vector<size_t> &a,\
		std::vector<std::vector<key_func_type>> &key_func,\
		std::vector<cond_func_type> &cond_func,\
		std::function<double(base_raw *, base_raw *, base_raw*)> result_func,\
		const double STEP,\
		const size_t &MAX,\
		const double &prob,\
		bool verbose,\
		double total_size){
	constexpr int K = 3;
	size_t max_round = 0;
	for(size_t i = 0; i < K; ++i){
		if(i == 0)
			max_round = table[i].size();
		else
			max_round = std::min(max_round, table[i].size());
	}
	std::vector<double> sample_values;

	tool::start_report();
	size_t test_count = 0;
	typedef unordered_multimap<uint32_t, size_t> index_type;
	std::vector<std::vector<index_type> > h_table;
	for(size_t i = 0; i < K; ++i){
		std::vector<index_type> tmp;
		h_table.push_back(tmp);
		for(size_t j = 0; j < K; ++j){
			index_type x;
			h_table[i].push_back(x);
		}
	}
	std::vector<double> q[K]; 
	std::vector<double> w(K, 0);
	std::vector<double> x(K, 0);
	for (size_t k = 0; k < K; ++k) q[k].resize(max_round * a[k]);
	for (size_t k = 0; k < K; ++k)
		for(size_t i = 0; i < q[k].size(); ++i)
			q[k][i] = 0;

	cout<<total_size<<endl;
	const double zp = tool::xql_erf_inv(prob);

	double y = 0.0;
	double ci = 1.0 / 0.0;
	size_t round = 0;
	/* round 1 */
	for(size_t i = 0; i < K; ++i){
		for(size_t k = 0; k < a[i]; ++k){
			for(size_t j = 0; j < K; ++j){
				if(key_func[i][j] != nullptr)
					if(cond_func[i] == nullptr || (*cond_func[i])(table[i][round*a[i] +k] ))
						h_table[i][j].insert(std::make_pair((*key_func[i][j])(table[i][round * a[i] + k]), round*a[i]+k)); 
			}
		}
	}
	double sum = 0;
	double v = 0;
	for(size_t i = 0; i < a[2]; ++i)
			sample_values.push_back(result_func(table[0][round], table[1][round], table[2][round*a[2] + i]));	
	for(size_t i0 = 0; i0 < a[0]; ++i0){
		for(size_t i1 = 0; i1 < a[1]; ++i1){
			for(size_t i2 = 0; i2 < a[2]; ++i2){
				std::vector<size_t> ind(K);
				ind[0] = round*a[0] + i0;
				ind[1] = round*a[1] + i1;
				ind[2] = round*a[2] + i2;
				if((*key_func[0][1])(table[0][ind[0]]) == (*key_func[1][0])(table[1][ind[1]]) &&
						(*key_func[1][2])(table[1][ind[1]]) == (*key_func[2][1])(table[2][ind[2]])){
						if((cond_func[0] == nullptr || (*cond_func[0])(table[0][ind[0]])) &&\
								(cond_func[1] == nullptr || (*cond_func[1])(table[1][ind[1]])) &&\
								(cond_func[2] == nullptr || (*cond_func[2])(table[2][ind[2]]))){
							v = sample_values[round*a[2] + i2];
							y += v *total_size / a[0] /a[1] / a[2];
							sum+=v;
							test_count++;
						}
				}
			}
		}
	}
	for (size_t k = 0; k < K; ++k) {
		q[k][0] = v;
		w[k] = 0;
		x[k] = y;
	}
	size_t n = 1;
	double z = 0;
	uint32_t n_rejected = 0;
	size_t report_round = 1;
	timer maxtimer;
	timer steptimer;
	round = 1;
	while(maxtimer.get_elapsed() < MAX && round < max_round){
		n = round + 1;
		for(size_t i = 0; i < K; ++i){
			for(size_t k = 0; k < a[i]; ++k){
				for(size_t j = 0; j < K; ++j){
					if(key_func[i][j] != nullptr)
						if(cond_func[i] == nullptr || (*cond_func[i])(table[i][round*a[i] +k] ))
							h_table[i][j].insert(std::make_pair((*key_func[i][j])(table[i][round * a[i] + k]), round*a[i]+k)); 
				}
			}
		}
		for(size_t i = 0; i < a[2]; ++i)
			sample_values.push_back(result_func(table[0][round], table[1][round], table[2][round*a[2] + i]));	
		y = y * std::pow(double(round)/(round+1), K-1);
		for(size_t i0 = 0; i0 < a[0]; ++i0){
			for(size_t i1 = 0; i1 < a[1]; ++i1){
				for(size_t i2 = 0; i2 < a[2]; ++i2){
					std::vector<size_t> ind(K);
					ind[0] = round*a[0] + i0;
					ind[1] = round*a[1] + i1;
					ind[2] = round*a[2] + i2;
					if((*key_func[0][1])(table[0][ind[0]]) == (*key_func[1][0])(table[1][ind[1]]) &&
							(*key_func[1][2])(table[1][ind[1]]) == (*key_func[2][1])(table[2][ind[2]])){
							if((cond_func[0] == nullptr || (*cond_func[0])(table[0][ind[0]])) &&\
									(cond_func[1] == nullptr || (*cond_func[1])(table[1][ind[1]])) &&\
									(cond_func[2] == nullptr || (*cond_func[2])(table[2][ind[2]]))){
								v = sample_values[round*a[2] + i2];
								y += v *total_size / a[0] /a[1] / a[2];
								sum+=v;
								test_count++;
							}
					}
				}
			}
		}

		for(size_t k = 0; k < K; ++k){
			q[k][round] = 0;
			w[k] = w[k] * std::pow(double(round) / (round+1), 4) + y * y /(round+1)/round;
			x[k] = y;
		}
		//ik == n from table 0
		for(size_t i = 0; i < a[0]; ++i)
		{
		std::vector<size_t> ind(K);
		ind[0] = round * a[0] + i;
		auto tuple0 = table[0][ind[0]];
		if(cond_func[0] == nullptr || cond_func[0] != nullptr && (*cond_func[0])(tuple0)){
			auto range1 = h_table[1][0].equal_range((*key_func[0][1])(tuple0));
			for(auto iter1 = range1.first; iter1 != range1.second; iter1++){
				ind[1] = iter1->second;
				if(ind[1]/a[1] >= round)
					continue;
				auto tuple1 = table[1][ind[1]];
				if(cond_func[1] != nullptr && (*cond_func[1])(tuple1) == false)
					continue;
				auto range2 = h_table[2][1].equal_range((*key_func[1][2])(tuple1));
				for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
					ind[2] = iter2->second;
					if(ind[2]/a[2] >= round)
						continue;
					auto tuple2 = table[2][ind[2]];
					if(cond_func[2] != nullptr && (*cond_func[2])(tuple2) == false)
						continue;
					calc_ripple_ci(K, ind, x, w, q, sample_values, n, a, total_size);
					sum += sample_values[ind[2]];
					test_count++;
				}
			}
			}
		}
		//ik == n from table 1
		for(size_t i = 0; i < a[1]; ++i)
		{
		std::vector<size_t> ind(K);
		ind[1] = round * a[1] + i;
		auto tuple1 = table[1][ind[1]];
		if(cond_func[1] == nullptr || (*cond_func[1])(tuple1)){
			auto range0 = h_table[0][1].equal_range((*key_func[1][0])(tuple1));
			for(auto iter0 = range0.first; iter0 != range0.second; iter0++){
				ind[0] = iter0->second;
				if(ind[0]/ a[0] >= round)
					continue;
				auto tuple0 = table[0][ind[0]];
				if(cond_func[0] != nullptr && !(*cond_func[0])(tuple0))
					continue;
				auto range2 = h_table[2][1].equal_range((*key_func[1][2])(tuple1));
				for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
					ind[2] = iter2->second;
					if(ind[2]/a[2] >= round)
						continue;
					auto tuple2 = table[2][ind[2]];
					if(cond_func[2] != nullptr && !(*cond_func[2])(tuple2))
						continue;
					calc_ripple_ci(K, ind, x, w, q, sample_values, n, a, total_size);
					sum += sample_values[ind[2]];
					test_count++;
				}
			}
		}
		}
		//ik == n from table 2
		for(size_t i = 0; i < a[2]; ++i)
		{
		std::vector<size_t> ind(K);
		ind[2] = round * a[2] + i;
		auto tuple2 = table[2][ind[2]];
		if(cond_func[2] == nullptr || (*cond_func[2])(tuple2)){
			auto range1 = h_table[1][2].equal_range((*key_func[2][1])(tuple2));
			for(auto iter1 = range1.first; iter1 != range1.second; iter1++){
				ind[1] = iter1->second;
				if(ind[1] /a[1] >= round)
					continue;
				auto tuple1 = table[1][ind[1]];
				if(cond_func[1] != nullptr && !(*cond_func[1])(tuple1))
					continue;
				auto range0 = h_table[0][1].equal_range((*key_func[1][0])(tuple1));
				for(auto iter0 = range0.first; iter0 != range0.second; ++iter0){
					ind[0] = iter0->second;
					if(ind[0]/a[0] >= round)
						continue;
					auto tuple0 = table[0][ind[0]];
					if(cond_func[0] != nullptr && !(*cond_func[0])(tuple0))
						continue;
					calc_ripple_ci(K, ind, x, w, q, sample_values, n, a, total_size);
					sum += sample_values[ind[2]];
					test_count++;
				}
			}
		}
		}
		y = x[0];
		z = 0;
		double nk=n * a[0] * n *a[1] * n *a[2];
		for(size_t i = 0; i < K; ++i)
			z += w[i] / (n*a[i]-1);
		if(steptimer.get_elapsed()> STEP){
			ci = zp*std::sqrt(z) / sqrt(n);
			tool::report(STEP*report_round, n, test_count, 0, 0, sum*total_size/nk, ci, prob);
			steptimer.restart();
			report_round++;
		}
		round++;
	}
}
std::pair<double, size_t> ripple_3_ratio(std::vector<std::vector<base_raw *> > &table,\
		std::vector<size_t> &a,\
		std::vector<std::vector<key_func_type>> &key_func,\
		std::vector<cond_func_type> &cond_func,\
		std::function<double(base_raw *, base_raw *, base_raw*)> result_func,\
		const double VALUE,\
		const double ERROR_LIMIT,\
		const double STEP,\
		const size_t &MAX,\
		const double &prob,\
		bool verbose){
	constexpr int K = 3;
	size_t max_round = 0;
	for(size_t i = 0; i < K; ++i){
		if(i == 0)
			max_round = table[i].size();
		else
			max_round = std::min(max_round, table[i].size());
	}
	std::vector<double> sample_values;
	double total_size = 1;
	for(size_t i = 0; i < K; ++i)
		total_size *= table[i].size();
	if(verbose)
		tool::start_report();
	size_t test_count = 0;
	typedef unordered_multimap<uint32_t, size_t> index_type;
	std::vector<std::vector<index_type> > h_table;
	for(size_t i = 0; i < K; ++i){
		std::vector<index_type> tmp;
		h_table.push_back(tmp);
		for(size_t j = 0; j < K; ++j){
			index_type x;
			h_table[i].push_back(x);
		}
	}
	std::vector<double> q[K]; 
	std::vector<double> w(K, 0);
	std::vector<double> x(K, 0);
	for (size_t k = 0; k < K; ++k) q[k].resize(max_round * a[k]);
	for (size_t k = 0; k < K; ++k)
		for(size_t i = 0; i < q[k].size(); ++i)
			q[k][i] = 0;

	const double zp = tool::xql_erf_inv(prob);

	double y = 0.0;
	double ci = 1.0 / 0.0;
	size_t round = 0;
	/* round 1 */
	for(size_t i = 0; i < K; ++i){
		for(size_t k = 0; k < a[i]; ++k){
			for(size_t j = 0; j < K; ++j){
				if(key_func[i][j] != nullptr)
					if(cond_func[i] == nullptr || (*cond_func[i])(table[i][round*a[i] +k] ))
						h_table[i][j].insert(std::make_pair((*key_func[i][j])(table[i][round * a[i] + k]), round*a[i]+k)); 
			}
		}
	}
	double sum = 0;
	double v = 0;
	for(size_t i = 0; i < a[2]; ++i)
			sample_values.push_back(result_func(table[0][round], table[1][round], table[2][round*a[2] + i]));	
	std::vector<size_t> ind(K);
	for(size_t i0 = 0; i0 < a[0]; ++i0){
		ind[0] = round*a[0] + i0;
		for(size_t i1 = 0; i1 < a[1]; ++i1){
			ind[1] = round*a[1] + i1;
			for(size_t i2 = 0; i2 < a[2]; ++i2){
				ind[2] = round*a[2] + i2;
				if((*key_func[0][1])(table[0][ind[0]]) == (*key_func[1][0])(table[1][ind[1]]) &&
						(*key_func[1][2])(table[1][ind[1]]) == (*key_func[2][1])(table[2][ind[2]])){
						if((cond_func[0] == nullptr || (*cond_func[0])(table[0][ind[0]])) &&\
								(cond_func[1] == nullptr || (*cond_func[1])(table[1][ind[1]])) &&\
								(cond_func[2] == nullptr || (*cond_func[2])(table[2][ind[2]]))){
							v = sample_values[round*a[2] + i2];
							y = v *total_size;
							q[0][ind[0]] += y;
							q[1][ind[1]] += y;
							q[2][ind[2]] += y;
							sum+=v;
							test_count++;
						}
				}
			}
		}
	}
	for(size_t i = 0; i < K; ++i){
		y = sum/a[i] * total_size;
	}
	double pow_n[K];
	double all = 1;
	for(size_t i = 0; i < K; ++i)
		all *= (round+1) * a[i];
	for(size_t i = 0; i < K; ++i)
		pow_n[i] = all/(round+1)/a[i];
	for(size_t i = 0; i < K; ++i){
		for(size_t j = 0; j < a[i]; ++j){
			x[i] += q[i][round + j];
		}
	}
	for(size_t i = 0; i < K; ++i){
		x[i] /= a[i];
	}
	for(size_t i = 0; i < K; ++i){
		w[i] = 0;
	}
	size_t n = 1;
	double z = 0;
	size_t report_round = 1;
	timer maxtimer;
	timer steptimer;
	round = 1;
	size_t old_time = 0;
	size_t cur_time;
	while(maxtimer.get_elapsed() < MAX && round < max_round && std::abs(y - VALUE) > ERROR_LIMIT){
		n = round + 1;
		for(size_t i = 0; i < K; ++i)
			all *= (round+1) * a[i];
		for(size_t i = 0; i < K; ++i)
			pow_n[i] = all/(round+1)/a[i];
		for(size_t i = 0; i < K; ++i){
			for(size_t k = 0; k < a[i]; ++k){
				for(size_t j = 0; j < K; ++j){
					if(key_func[i][j] != nullptr)
						if(cond_func[i] == nullptr || (*cond_func[i])(table[i][round*a[i] +k] ))
							h_table[i][j].insert(std::make_pair((*key_func[i][j])(table[i][round * a[i] + k]), round*a[i]+k)); 
				}
			}
		}
		for(size_t i = 0; i < a[2]; ++i)
			sample_values.push_back(result_func(table[0][round], table[1][round], table[2][round*a[2] + i]));	
		y = y * std::pow(double(round)/(round+1), K-1);
		std::vector<size_t> ind(K);
		for(size_t i0 = 0; i0 < a[0]; ++i0){
			ind[0] = round*a[0] + i0;
			for(size_t i1 = 0; i1 < a[1]; ++i1){
				ind[1] = round*a[1] + i1;
				for(size_t i2 = 0; i2 < a[2]; ++i2){
					ind[2] = round*a[2] + i2;
					if((*key_func[0][1])(table[0][ind[0]]) == (*key_func[1][0])(table[1][ind[1]]) &&
							(*key_func[1][2])(table[1][ind[1]]) == (*key_func[2][1])(table[2][ind[2]])){
							if((cond_func[0] == nullptr || (*cond_func[0])(table[0][ind[0]])) &&\
									(cond_func[1] == nullptr || (*cond_func[1])(table[1][ind[1]])) &&\
									(cond_func[2] == nullptr || (*cond_func[2])(table[2][ind[2]]))){
								v = sample_values[round*a[2] + i2];
								auto tmp = v *total_size;
								q[0][ind[0]] += tmp;
								q[1][ind[1]] += tmp;
								q[2][ind[2]] += tmp;
								sum+=v;
								test_count++;
							}
					}
				}
			}
		}
		for(size_t i = 0; i < K; ++i){
			auto tmp = 0;
			for(size_t j = 0; j < a[i]; ++j){
				tmp += q[i][round + j];
			}
			x[i] = y + tmp/a[i]/ pow_n[i];
		}

		for(size_t k = 0; k < K; ++k){
			q[k][round] = 0;
			w[k] = w[k] * std::pow(double(round) / (round+1), 4) + y * y /(round+1)/round;
		}
		//ik == n from table 0
		for(size_t i = 0; i < a[0]; ++i)
		{
		std::vector<size_t> ind(K);
		ind[0] = round * a[0] + i;
		auto tuple0 = table[0][ind[0]];
		if(cond_func[0] == nullptr || cond_func[0] != nullptr && (*cond_func[0])(tuple0)){
			auto range1 = h_table[1][0].equal_range((*key_func[0][1])(tuple0));
			for(auto iter1 = range1.first; iter1 != range1.second; iter1++){
				ind[1] = iter1->second;
				if(ind[1]/a[1] >= round)
					continue;
				auto tuple1 = table[1][ind[1]];
				auto range2 = h_table[2][1].equal_range((*key_func[1][2])(tuple1));
				for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
					ind[2] = iter2->second;
					if(ind[2]/a[2] >= round)
						continue;
					calc_ripple_ci_test(K, ind, x, w, q, sample_values, n, a, total_size);
					sum += sample_values[ind[2]];
					test_count++;
				}
			}
			}
		}
		//ik == n from table 1
		for(size_t i = 0; i < a[1]; ++i)
		{
		std::vector<size_t> ind(K);
		ind[1] = round * a[1] + i;
		auto tuple1 = table[1][ind[1]];
		if(cond_func[1] == nullptr || (*cond_func[1])(tuple1)){
			auto range0 = h_table[0][1].equal_range((*key_func[1][0])(tuple1));
			for(auto iter0 = range0.first; iter0 != range0.second; iter0++){
				ind[0] = iter0->second;
				if(ind[0]/ a[0] >= round)
					continue;
				auto range2 = h_table[2][1].equal_range((*key_func[1][2])(tuple1));
				for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
					ind[2] = iter2->second;
					if(ind[2]/a[2] >= round)
						continue;
					calc_ripple_ci_test(K, ind, x, w, q, sample_values, n, a, total_size);
					sum += sample_values[ind[2]];
					test_count++;
				}
			}
		}
		}
		//ik == n from table 2
		for(size_t i = 0; i < a[2]; ++i)
		{
		std::vector<size_t> ind(K);
		ind[2] = round * a[2] + i;
		auto tuple2 = table[2][ind[2]];
		if(cond_func[2] == nullptr || (*cond_func[2])(tuple2)){
			auto range1 = h_table[1][2].equal_range((*key_func[2][1])(tuple2));
			for(auto iter1 = range1.first; iter1 != range1.second; iter1++){
				ind[1] = iter1->second;
				if(ind[1] /a[1] >= round)
					continue;
				auto tuple1 = table[1][ind[1]];
				auto range0 = h_table[0][1].equal_range((*key_func[1][0])(tuple1));
				for(auto iter0 = range0.first; iter0 != range0.second; ++iter0){
					ind[0] = iter0->second;
					if(ind[0]/a[0] >= round)
						continue;
					auto tuple0 = table[0][ind[0]];
					calc_ripple_ci_test(K, ind, x, w, q, sample_values, n, a, total_size);
					sum += sample_values[ind[2]];
					test_count++;
				}
			}
		}
		}
		z = 0;
		double nk=n * a[0] * n *a[1] * n *a[2];
		y = sum*total_size/nk;
		for(size_t i = 0; i < K; ++i)
			z += w[i] / (n*a[i]-1);
		cur_time = maxtimer.get_elapsed();
		if(verbose && cur_time - old_time> STEP){
			ci = zp*std::sqrt(z) / sqrt(n);
			tool::report(cur_time, n, test_count, 0, 0, sum*total_size/nk, ci, prob);
			old_time = cur_time;
			report_round++;
		}
		round++;
	}
	return std::make_pair(y, cur_time);
}
std::pair<double, size_t> ripple_4_ratio(std::vector<std::vector<base_raw *> > &table,\
		std::vector<size_t> &a,\
		std::vector<std::vector<key_func_type>> &key_func,\
		std::vector<cond_func_type> &cond_func,\
		std::function<double(base_raw*)> result_func,\
		const double VALUE,\
		const double ERROR_LIMIT,\
		const double STEP,\
		const size_t &MAX,\
		const double &prob,\
		bool verbose,\
		size_t col_id){
	constexpr int K = 4;
	size_t max_round = 0;
	for(size_t i = 0; i < K; ++i){
		if(i == 0)
			max_round = table[i].size();
		else
			max_round = std::min(max_round, table[i].size());
	}
	std::vector<double> sample_values;
	double total_size = 1;
	for(size_t i = 0; i < K; ++i)
		total_size *= table[i].size();
	if(verbose)
		tool::start_report();
	size_t test_count = 0;
	typedef unordered_multimap<uint32_t, size_t> index_type;
	std::vector<std::vector<index_type> > h_table;
	for(size_t i = 0; i < K; ++i){
		std::vector<index_type> tmp;
		h_table.push_back(tmp);
		for(size_t j = 0; j < K; ++j){
			index_type x;
			h_table[i].push_back(x);
		}
	}
	std::vector<double> q[K]; 
	std::vector<double> w(K, 0);
	std::vector<double> x(K, 0);
	for (size_t k = 0; k < K; ++k) q[k].resize(max_round * a[k]);
	for (size_t k = 0; k < K; ++k)
		for(size_t i = 0; i < q[k].size(); ++i)
			q[k][i] = 0;


	const double zp = tool::xql_erf_inv(prob);

	double y = 0.0;
	double ci = 1.0 / 0.0;
	size_t round = 0;
	/* round 0 */
	for(size_t i = 0; i < K; ++i){
		for(size_t k = 0; k < a[i]; ++k){
			for(size_t j = 0; j < K; ++j){
					if(key_func[i][j] != nullptr)
						if(cond_func[i] == nullptr || (*cond_func[i])(table[i][round*a[i] +k]))
							h_table[i][j].insert(std::make_pair((*key_func[i][j])(table[i][round * a[i] + k]), round*a[i]+k)); 
			}
		}
	}
	double sum = 0;
	double v = 0;
	for(size_t i = 0; i < a[col_id]; ++i)
			sample_values.push_back(result_func(table[col_id][round*a[col_id] + i]));	
	size_t i[4];
	for(i[0] = 0; i[0] < a[0]; ++i[0]){
		for(i[1] = 0; i[1] < a[1]; ++i[1]){
			for(i[2] = 0; i[2] < a[2]; ++i[2]){
				for(i[3] = 0; i[3] < a[3]; ++i[3]){
					std::vector<size_t> ind(K);
					ind[0] = round*a[0] + i[0];
					ind[1] = round*a[1] + i[1];
					ind[2] = round*a[2] + i[2];
					ind[3] = round*a[3] + i[3];
				if((*key_func[0][1])(table[0][ind[0]]) == (*key_func[1][0])(table[1][ind[1]]) \
						&& (*key_func[1][2])(table[1][ind[1]]) == (*key_func[2][1])(table[2][ind[2]]) \
						&& (*key_func[2][3])(table[2][ind[2]])== (*key_func[3][2])(table[3][ind[3]])){
						if((cond_func[0] == nullptr || (*cond_func[0])(table[0][ind[0]])) &&\
								(cond_func[1] == nullptr || (*cond_func[1])(table[1][ind[1]])) &&\
								(cond_func[2] == nullptr || (*cond_func[2])(table[2][ind[2]])) &&\
								 (cond_func[3]==nullptr || (*cond_func[3])(table[3][ind[3]]))){
							v = sample_values[round*a[col_id] + i[col_id]];
							y = v * total_size;
							q[0][ind[0]] += y;
							q[1][ind[1]] += y;
							q[2][ind[2]] += y;
							q[3][ind[3]] += y;
							sum+=v;
							test_count++;
						}
				}
				}
			}
		}
	}
	for(size_t i = 0; i < K; ++i){
		y = sum/a[i] * total_size;
	}
	double pow_n[K];
	double all = 1;
	for(size_t i = 0; i < K; ++i)
		all *= (round+1) * a[i];
	for(size_t i = 0; i < K; ++i)
		pow_n[i] = all/(round+1)/a[i];
	for(size_t i = 0; i < K; ++i){
		for(size_t j = 0; j < a[i]; ++j){
			x[i] += q[i][round*a[i]+ j]/a[i];
		}
	}
	for(size_t i = 0; i < K; ++i){
		w[i] = 0;
	}
	size_t n = 1;
	double z = 0;
	uint32_t n_rejected = 0;
	size_t report_round = 1;
	timer maxtimer;
	timer steptimer;
	size_t old_time = 0;
	round = 1;
	while(maxtimer.get_elapsed() < MAX && round < max_round && std::abs(y - VALUE) > ERROR_LIMIT){
		n = round + 1;
		for(size_t i = 0; i < K; ++i)
			all *= (round+1) * a[i];
		for(size_t i = 0; i < K; ++i)
			pow_n[i] = all/(round+1)/a[i];
		for(size_t i = 0; i < K; ++i){
			for(size_t k = 0; k < a[i]; ++k){
				for(size_t j = 0; j < K; ++j){
					if(key_func[i][j] != nullptr){
						if(cond_func[i] == nullptr || (*cond_func[i])(table[i][round*a[i] +k] ))
							h_table[i][j].insert(std::make_pair((*key_func[i][j])(table[i][round * a[i] + k]), round*a[i]+k)); 
					}
				}
			}
		}
		for(size_t i = 0; i < a[col_id]; ++i)
			sample_values.push_back(result_func(table[col_id][round*a[col_id] + i]));	
		y = y * std::pow(double(round)/(round+1), K-1);
		for(size_t k = 0; k < K; ++k){
			for(size_t i = 0; i < a[k]; ++i)
				q[k][round*a[k] + i ] = 0;
		}
		for(i[0] = 0; i[0] < a[0]; ++i[0]){
			for(i[1] = 0; i[1] < a[1]; ++i[1]){
				for(i[2] = 0; i[2] < a[2]; ++i[2]){
					for(i[3] = 0; i[3] < a[3]; ++i[3]){
						std::vector<size_t> ind(K);
						ind[0] = round*a[0] + i[0];
						ind[1] = round*a[1] + i[1];
						ind[2] = round*a[2] + i[2];
						ind[3] = round*a[3] + i[3];
					if((*key_func[0][1])(table[0][ind[0]]) == (*key_func[1][0])(table[1][ind[1]]) \
							&& (*key_func[1][2])(table[1][ind[1]]) == (*key_func[2][1])(table[2][ind[2]]) \
							&& (*key_func[2][3])(table[2][ind[2]])== (*key_func[3][2])(table[3][ind[3]])){
							if((cond_func[0] == nullptr || (*cond_func[0])(table[0][ind[0]])) &&\
									(cond_func[1] == nullptr || (*cond_func[1])(table[1][ind[1]])) &&\
									(cond_func[2] == nullptr || (*cond_func[2])(table[2][ind[2]])) &&\
									 (cond_func[3]==nullptr || (*cond_func[3])(table[3][ind[3]]))){
								v = sample_values[round*a[col_id] + i[col_id]];
								auto tmp = v *total_size;
								q[0][ind[0]] += tmp;
								q[1][ind[1]] += tmp;
								q[2][ind[2]] += tmp;
								q[3][ind[3]] += tmp;
								sum += v;
								test_count++;
							}
					}
					}
				}
			}
		}

		for(size_t i = 0; i < K; ++i){
			auto tmp = 0;
			for(size_t j = 0; j < a[i]; ++j){
				tmp += q[i][round*a[i] + j];
			}
			x[i] = y + tmp/a[i]/ pow_n[i];
		}

		for(size_t k = 0; k < K; ++k){
			for(size_t i = 0; i < a[k]; ++i)
				q[k][round*a[k] + i ] = 0;
			w[k] = w[k] * std::pow(double(round) / (round+1), 6) + y * y /(round+1)/round;
		}
		//ik == n from table 0
		for(size_t i = 0; i < a[0]; ++i)
		{
		std::vector<size_t> ind(K);
		ind[0] = round*a[0] + i;
		auto tuple0 = table[0][ind[0]];
		if(cond_func[0] == nullptr || cond_func[0] != nullptr && (*cond_func[0])(tuple0)){
			auto range1 = h_table[1][0].equal_range((*key_func[0][1])(tuple0));
			for(auto iter1 = range1.first; iter1 != range1.second; iter1++){
				ind[1] = iter1->second;
				if(ind[1]/a[1] >= round)
					continue;
				auto tuple1 = table[1][ind[1]];
				auto range2 = h_table[2][1].equal_range((*key_func[1][2])(tuple1));
				for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
					ind[2] = iter2->second;
					if(ind[2]/a[2] >= round)
						continue;
					auto tuple2 = table[2][ind[2]];
					auto range3 = h_table[3][2].equal_range((*key_func[2][3])(tuple2));
					for(auto iter3 = range3.first; iter3 != range3.second; ++iter3){
						ind[3] = iter3->second;
						if(ind[3]/a[3] >= round)
							continue;
						auto tuple3 = table[3][ind[3]];
						calc_ripple_ci_test(K, ind, x, w, q, sample_values, n, a, total_size, col_id);
						sum += sample_values[ind[col_id]];
						test_count++;
					}
				}
			}
			}
		}
		//ik == n from table 1
		for(size_t i = 0; i < a[1]; ++i)
		{
		std::vector<size_t> ind(K);
		ind[1] = round * a[1] + i;
		auto tuple1 = table[1][ind[1]];
		if(cond_func[1] == nullptr || (*cond_func[1])(tuple1)){
			auto range0 = h_table[0][1].equal_range((*key_func[1][0])(tuple1));
			for(auto iter0 = range0.first; iter0 != range0.second; iter0++){
				ind[0] = iter0->second;
				if(ind[0]/a[0]>= round)
					continue;
				auto tuple0 = table[0][ind[0]];
				auto range2 = h_table[2][1].equal_range((*key_func[1][2])(tuple1));
				for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
					ind[2] = iter2->second;
					if(ind[2]/a[2]>= round)
						continue;
					auto tuple2 = table[2][ind[2]];
					auto range3 = h_table[3][2].equal_range((*key_func[2][3])(tuple2));
					for(auto iter3 = range3.first; iter3 != range3.second; ++iter3){
						ind[3] = iter3->second;
						if(ind[3]/a[3] >= round)
							continue;
						auto tuple3 = table[3][ind[3]];
						calc_ripple_ci_test(K, ind, x, w, q, sample_values, n, a, total_size, col_id);
						sum += sample_values[ind[col_id]];
						test_count++;
					}
				}
			}
		}
		}
		//ik == n from table 2
		for(size_t i = 0; i < a[2]; ++i)
		{
		std::vector<size_t> ind(K);
		ind[2] = a[2]*round + i;
		auto tuple2 = table[2][ind[2]];
		if(cond_func[2] == nullptr || (*cond_func[2])(tuple2)){
			auto range3 = h_table[3][2].equal_range((*key_func[2][3])(tuple2));
			for(auto iter3 = range3.first; iter3 != range3.second; ++iter3){
				ind[3] = iter3->second;
				if(ind[3]/a[3]>= round)
					continue;
				auto tuple3 = table[3][ind[3]];
				auto range1 = h_table[1][2].equal_range((*key_func[2][1])(tuple2));
				for(auto iter1 = range1.first; iter1 != range1.second; iter1++){
					ind[1] = iter1->second;
					if(ind[1]/a[1] >= round)
						continue;
					auto tuple1 = table[1][ind[1]];
					auto range0 = h_table[0][1].equal_range((*key_func[1][0])(tuple1));
					for(auto iter0 = range0.first; iter0 != range0.second; ++iter0){
						ind[0] = iter0->second;
						if(ind[0]/a[0] >= round)
							continue;
						auto tuple0 = table[0][ind[0]];
						calc_ripple_ci_test(K, ind, x, w, q, sample_values, n, a, total_size, col_id);
						sum += sample_values[ind[col_id]];
						test_count++;
					}
				}
			}
		}
		}
		//! i_n from table 3
		for(size_t i = 0; i < a[3]; ++i)
		{
		std::vector<size_t> ind(K);
		ind[3] = round*a[3] + i;
		auto tuple3 = table[3][ind[3]];
		if(cond_func[3] == nullptr || (*cond_func[3])(tuple3)){
			auto range2 = h_table[2][3].equal_range((*key_func[3][2])(tuple3));
			for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
				ind[2] = iter2->second;
				if(ind[2]/a[2]>= round)
					continue;
				auto tuple2 = table[2][ind[2]];
				auto range1 = h_table[1][2].equal_range((*key_func[2][1])(tuple2));
				for(auto iter1 = range1.first; iter1 != range1.second; iter1++){
					ind[1] = iter1->second;
					if(ind[1]/a[1]>= round)
						continue;
					auto tuple1 = table[1][ind[1]];
					auto range0 = h_table[0][1].equal_range((*key_func[1][0])(tuple1));
					for(auto iter0 = range0.first; iter0 != range0.second; ++iter0){
						ind[0] = iter0->second;
						if(ind[0]/a[0] >= round)
							continue;
						auto tuple0 = table[0][ind[0]];
						calc_ripple_ci_test(K, ind, x, w, q, sample_values, n, a, total_size, col_id);
						sum += sample_values[ind[col_id]];
						test_count++;
					}
				}
			}
		}
		}
		z = 0;
		double nk = 1;
		for(size_t i = 0; i < K; ++i)
			nk *= a[i] * n;
		y = sum/nk*total_size;
		for(size_t i = 0; i < K; ++i)
			z += w[i] / (n * a[i]-1);
		auto cur_time = maxtimer.get_elapsed();
		if(verbose && cur_time - old_time> STEP){
			ci = zp*std::sqrt(z) / sqrt(n);
			tool::report(cur_time, n, test_count, 0, 0, sum/nk*total_size, ci, prob);
			old_time = cur_time;
			report_round++;
		}
		round++;
	}
	return std::make_pair(y,maxtimer.get_elapsed());
}
void ripple_3_ratio_test(std::vector<std::vector<base_raw *> > &table,\
		std::vector<size_t> &a,\
		std::vector<std::vector<key_func_type>> &key_func,\
		std::vector<cond_func_type> &cond_func,\
		std::function<double(base_raw *, base_raw *, base_raw*)> result_func,\
		const double STEP,\
		const size_t &MAX,\
		const double &prob,\
		bool verbose,\
		const double &VALUE,\
		double total_size){
	constexpr int K = 3;
	size_t max_round = 0;
	for(size_t i = 0; i < K; ++i){
		if(i == 0)
			max_round = table[i].size();
		else
			max_round = std::min(max_round, table[i].size());
	}
	std::vector<double> sample_values;

	size_t test_count = 0;
	typedef unordered_multimap<uint32_t, size_t> index_type;
	std::vector<std::vector<index_type> > h_table;
	for(size_t i = 0; i < K; ++i){
		std::vector<index_type> tmp;
		h_table.push_back(tmp);
		for(size_t j = 0; j < K; ++j){
			index_type x;
			h_table[i].push_back(x);
		}
	}
	std::vector<double> q[K]; 
	std::vector<double> w(K, 0);
	std::vector<double> x(K, 0);
	for (size_t k = 0; k < K; ++k) q[k].resize(max_round * a[k]);
	for (size_t k = 0; k < K; ++k)
		for(size_t i = 0; i < q[k].size(); ++i)
			q[k][i] = 0;

	cout<<total_size<<endl;
	const double zp = tool::xql_erf_inv(prob);

	double y = 0.0;
	double ci = 1.0 / 0.0;
	size_t round = 0;
	/* round 1 */
	for(size_t i = 0; i < K; ++i){
		for(size_t k = 0; k < a[i]; ++k){
			for(size_t j = 0; j < K; ++j){
				if(key_func[i][j] != nullptr)
					if(cond_func[i] == nullptr || (*cond_func[i])(table[i][round*a[i] +k] ))
						h_table[i][j].insert(std::make_pair((*key_func[i][j])(table[i][round * a[i] + k]), round*a[i]+k)); 
			}
		}
	}
	double sum = 0;
	double v = 0;
	for(size_t i = 0; i < a[2]; ++i)
			sample_values.push_back(result_func(table[0][round], table[1][round], table[2][round*a[2] + i]));	
	std::vector<size_t> ind(K);
	for(size_t i0 = 0; i0 < a[0]; ++i0){
		ind[0] = round*a[0] + i0;
		for(size_t i1 = 0; i1 < a[1]; ++i1){
			ind[1] = round*a[1] + i1;
			for(size_t i2 = 0; i2 < a[2]; ++i2){
				ind[2] = round*a[2] + i2;
				if((*key_func[0][1])(table[0][ind[0]]) == (*key_func[1][0])(table[1][ind[1]]) &&
						(*key_func[1][2])(table[1][ind[1]]) == (*key_func[2][1])(table[2][ind[2]])){
						if((cond_func[0] == nullptr || (*cond_func[0])(table[0][ind[0]])) &&\
								(cond_func[1] == nullptr || (*cond_func[1])(table[1][ind[1]])) &&\
								(cond_func[2] == nullptr || (*cond_func[2])(table[2][ind[2]]))){
							v = sample_values[round*a[2] + i2];
							y = v *total_size;
							q[0][ind[0]] += y;
							q[1][ind[1]] += y;
							q[2][ind[2]] += y;
							sum+=v;
							test_count++;
						}
				}
			}
		}
	}
	for(size_t i = 0; i < K; ++i){
		y = sum/a[i] * total_size;
	}
	double pow_n[K];
	double all = 1;
	for(size_t i = 0; i < K; ++i)
		all *= (round+1) * a[i];
	for(size_t i = 0; i < K; ++i)
		pow_n[i] = all/(round+1)/a[i];
	for(size_t i = 0; i < K; ++i){
		for(size_t j = 0; j < a[i]; ++j){
			x[i] += q[i][round + j];
		}
	}
	for(size_t i = 0; i < K; ++i){
		x[i] /= a[i];
	}
	for(size_t i = 0; i < K; ++i){
		w[i] = 0;
	}
	size_t n = 1;
	double z = 0;
	uint32_t n_rejected = 0;
	size_t report_round = 1;
	timer maxtimer;
	timer steptimer;
	round = 1;
	size_t old_time = 0;
	size_t cur_time;
	while(maxtimer.get_elapsed() < MAX && round < max_round){
		n = round + 1;
		for(size_t i = 0; i < K; ++i)
			all *= (round+1) * a[i];
		for(size_t i = 0; i < K; ++i)
			pow_n[i] = all/(round+1)/a[i];
		for(size_t i = 0; i < K; ++i){
			for(size_t k = 0; k < a[i]; ++k){
				for(size_t j = 0; j < K; ++j){
					if(key_func[i][j] != nullptr)
						if(cond_func[i] == nullptr || (*cond_func[i])(table[i][round*a[i] +k] ))
							h_table[i][j].insert(std::make_pair((*key_func[i][j])(table[i][round * a[i] + k]), round*a[i]+k)); 
				}
			}
		}
		for(size_t i = 0; i < a[2]; ++i)
			sample_values.push_back(result_func(table[0][round], table[1][round], table[2][round*a[2] + i]));	
		y = y * std::pow(double(round)/(round+1), K-1);
		std::vector<size_t> ind(K);
		for(size_t i0 = 0; i0 < a[0]; ++i0){
			ind[0] = round*a[0] + i0;
			for(size_t i1 = 0; i1 < a[1]; ++i1){
				ind[1] = round*a[1] + i1;
				for(size_t i2 = 0; i2 < a[2]; ++i2){
					ind[2] = round*a[2] + i2;
					if((*key_func[0][1])(table[0][ind[0]]) == (*key_func[1][0])(table[1][ind[1]]) &&
							(*key_func[1][2])(table[1][ind[1]]) == (*key_func[2][1])(table[2][ind[2]])){
							if((cond_func[0] == nullptr || (*cond_func[0])(table[0][ind[0]])) &&\
									(cond_func[1] == nullptr || (*cond_func[1])(table[1][ind[1]])) &&\
									(cond_func[2] == nullptr || (*cond_func[2])(table[2][ind[2]]))){
								v = sample_values[round*a[2] + i2];
								auto tmp = v *total_size;
								q[0][ind[0]] += tmp;
								q[1][ind[1]] += tmp;
								q[2][ind[2]] += tmp;
								sum+=v;
								test_count++;
							}
					}
				}
			}
		}
		for(size_t i = 0; i < K; ++i){
			auto tmp = 0;
			for(size_t j = 0; j < a[i]; ++j){
				tmp += q[i][round + j];
			}
			x[i] = y + tmp/a[i]/ pow_n[i];
		}

		for(size_t k = 0; k < K; ++k){
			q[k][round] = 0;
			w[k] = w[k] * std::pow(double(round) / (round+1), 4) + y * y /(round+1)/round;
		}
		//ik == n from table 0
		for(size_t i = 0; i < a[0]; ++i)
		{
		std::vector<size_t> ind(K);
		ind[0] = round * a[0] + i;
		auto tuple0 = table[0][ind[0]];
		if(cond_func[0] == nullptr || cond_func[0] != nullptr && (*cond_func[0])(tuple0)){
			auto range1 = h_table[1][0].equal_range((*key_func[0][1])(tuple0));
			for(auto iter1 = range1.first; iter1 != range1.second; iter1++){
				ind[1] = iter1->second;
				if(ind[1]/a[1] >= round)
					continue;
				auto tuple1 = table[1][ind[1]];
				auto range2 = h_table[2][1].equal_range((*key_func[1][2])(tuple1));
				for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
					ind[2] = iter2->second;
					if(ind[2]/a[2] >= round)
						continue;
					auto tuple2 = table[2][ind[2]];
					calc_ripple_ci_test(K, ind, x, w, q, sample_values, n, a, total_size);
					sum += sample_values[ind[2]];
					test_count++;
				}
			}
			}
		}
		//ik == n from table 1
		for(size_t i = 0; i < a[1]; ++i)
		{
		std::vector<size_t> ind(K);
		ind[1] = round * a[1] + i;
		auto tuple1 = table[1][ind[1]];
		if(cond_func[1] == nullptr || (*cond_func[1])(tuple1)){
			auto range0 = h_table[0][1].equal_range((*key_func[1][0])(tuple1));
			for(auto iter0 = range0.first; iter0 != range0.second; iter0++){
				ind[0] = iter0->second;
				if(ind[0]/ a[0] >= round)
					continue;
				auto tuple0 = table[0][ind[0]];
				auto range2 = h_table[2][1].equal_range((*key_func[1][2])(tuple1));
				for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
					ind[2] = iter2->second;
					if(ind[2]/a[2] >= round)
						continue;
					auto tuple2 = table[2][ind[2]];
					calc_ripple_ci_test(K, ind, x, w, q, sample_values, n, a, total_size);
					sum += sample_values[ind[2]];
					test_count++;
				}
			}
		}
		}
		//ik == n from table 2
		for(size_t i = 0; i < a[2]; ++i)
		{
		std::vector<size_t> ind(K);
		ind[2] = round * a[2] + i;
		auto tuple2 = table[2][ind[2]];
		if(cond_func[2] == nullptr || (*cond_func[2])(tuple2)){
			auto range1 = h_table[1][2].equal_range((*key_func[2][1])(tuple2));
			for(auto iter1 = range1.first; iter1 != range1.second; iter1++){
				ind[1] = iter1->second;
				if(ind[1] /a[1] >= round)
					continue;
				auto tuple1 = table[1][ind[1]];
				auto range0 = h_table[0][1].equal_range((*key_func[1][0])(tuple1));
				for(auto iter0 = range0.first; iter0 != range0.second; ++iter0){
					ind[0] = iter0->second;
					if(ind[0]/a[0] >= round)
						continue;
					auto tuple0 = table[0][ind[0]];
					calc_ripple_ci_test(K, ind, x, w, q, sample_values, n, a, total_size);
					sum += sample_values[ind[2]];
					test_count++;
				}
			}
		}
		}
		y = x[0];
		z = 0;
		double nk=n * a[0] * n *a[1] * n *a[2];
		for(size_t i = 0; i < K; ++i)
			z += w[i] / (n*a[i]-1);
		cur_time = maxtimer.get_elapsed();
		ci = zp*std::sqrt(z) / sqrt(n);
		if(ci > 0 && ci < VALUE)
			return;
		if(cur_time - old_time> STEP){
			tool::report(cur_time, n, test_count, 0, 0, sum*total_size/nk, ci, prob);
			old_time = cur_time;
			report_round++;
		}
		round++;
	}
}
void ripple_3_ratio_test(std::vector<std::vector<base_raw *> > &table,\
		std::vector<size_t> &a,\
		std::vector<std::vector<key_func_type>> &key_func,\
		std::vector<cond_func_type> &cond_func,\
		std::function<double(base_raw *, base_raw *, base_raw*)> result_func,\
		const double STEP,\
		const size_t &MAX,\
		const double &prob,\
		bool verbose,\
		double total_size){
	ripple_3_ratio_test(table, a, key_func, cond_func, result_func, STEP, MAX, prob, verbose, -1, total_size);
}
void ripple_4(std::vector<std::vector<base_raw *> > &table,\
		std::vector<std::vector<key_func_type>> &key_func,\
		std::vector<cond_func_type> &cond_func,\
		std::function<double(base_raw *)> result_func,
		const double STEP,\
		const size_t &MAX,\
		const double &prob,\
		bool verbose,\
		double total_size){
	constexpr int K = 4;
	size_t max_round = 0;
	for(size_t i = 0; i < K; ++i){
		if(i == 0)
			max_round = table[i].size();
		else
			max_round = std::min(max_round, table[i].size());
	}
	std::vector<double> sample_values;

	tool::start_report();
	size_t test_count = 0;
	typedef unordered_multimap<uint32_t, size_t> index_type;
	std::vector<std::vector<index_type> > h_table;
	for(size_t i = 0; i < K; ++i){
		std::vector<index_type> tmp;
		h_table.push_back(tmp);
		for(size_t j = 0; j < K; ++j){
			index_type x;
			h_table[i].push_back(x);
		}
	}
	std::vector<double> q[K]; 
	std::vector<double> w(K, 0);
	std::vector<double> x(K, 0);
	for (size_t k = 0; k < K; ++k) q[k].resize(max_round);

	cout<<total_size<<endl;
	const double zp = tool::xql_erf_inv(prob);

	double y = 0.0;
	double ci = 1.0 / 0.0;
	size_t round = 0;
	/* round 1 */
	for(size_t i = 0; i < K; ++i){
		for(size_t j = 0; j < K; ++j){
			if(key_func[i][j] != nullptr)
					if(cond_func[i] == nullptr || (*cond_func[i])(table[i][round]))
				h_table[i][j].insert(std::make_pair((*key_func[i][j])(table[i][round]), round)); 
		}
	}
	double sum = 0;
	double v = result_func(table[3][round]);
	sample_values.push_back(v);
	if((*key_func[0][1])(table[0][round]) == (*key_func[1][0])(table[1][round]) \
			&& (*key_func[1][2])(table[1][round]) == (*key_func[2][1])(table[2][round]) \
			&& (*key_func[2][3])(table[2][round])== (*key_func[3][2])(table[3][round])){
			if((cond_func[0] == nullptr || (*cond_func[0])(table[0][round])) &&\
					(cond_func[1] == nullptr || (*cond_func[1])(table[1][round])) &&\
					(cond_func[2] == nullptr || (*cond_func[2])(table[2][round])) &&\
					 (cond_func[3]==nullptr || (*cond_func[3])(table[3][round]))){
				y = v *total_size;
				sum+=v;
				test_count++;
			}
	}
	for (size_t k = 0; k < K; ++k) {
		q[k][0] = v;
		w[k] = 0;
		x[k] = y;
	}
	size_t n = 1;
	double z = 0;
	uint32_t n_rejected = 0;
	size_t report_round = 1;
	timer maxtimer;
	timer steptimer;
	round = 1;
	while(maxtimer.get_elapsed() < MAX && round < max_round){
		n = round + 1;
		for(size_t i = 0; i < K; ++i){
			for(size_t j = 0; j < K; ++j){
				if(key_func[i][j] != nullptr){
					if(cond_func[i] == nullptr || (*cond_func[i])(table[i][round]))
					h_table[i][j].insert(std::make_pair((*key_func[i][j])(table[i][round]), round)); 
				}
			}
		}
		sample_values.push_back(result_func(table[3][round]));	
		y = y * std::pow(double(round)/(round+1), K-1);
		base_raw *e[K];
		v = sample_values.back();
		if((*key_func[0][1])(table[0][round]) == (*key_func[1][0])(table[1][round]) \
				&& (*key_func[1][2])(table[1][round]) == (*key_func[2][1])(table[2][round])\
				&& (*key_func[2][3])(table[2][round]) == (*key_func[3][2])(table[3][round])\
				){
			if((cond_func[0] == nullptr || (*cond_func[0])(table[0][round])) &&\
					(cond_func[1] == nullptr || (*cond_func[1])(table[1][round])) &&\
					(cond_func[2] == nullptr || (*cond_func[2])(table[2][round]))&&\
					 (cond_func[3]==nullptr || (*cond_func[3])(table[3][round]))){
				y = v *total_size;
				sum+=v;
				test_count++;
			}
		}

		for(size_t k = 0; k < K; ++k){
			q[k][round] = 0;
			w[k] = w[k] * std::pow(double(round) / (round+1), 6) + y * y /(round+1)/round;
			x[k] = y;
		}
		//ik == n from table 0
		{
		std::vector<size_t> ind(K);
		ind[0] = round;
		auto tuple0 = table[0][ind[0]];
		if(cond_func[0] == nullptr || cond_func[0] != nullptr && (*cond_func[0])(tuple0)){
			auto range1 = h_table[1][0].equal_range((*key_func[0][1])(tuple0));
			for(auto iter1 = range1.first; iter1 != range1.second; iter1++){
				ind[1] = iter1->second;
				if(ind[1] >= round)
					continue;
				auto tuple1 = table[1][ind[1]];
				if(cond_func[1] != nullptr && (*cond_func[1])(tuple1) == false)
					continue;
				auto range2 = h_table[2][1].equal_range((*key_func[1][2])(tuple1));
				for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
					ind[2] = iter2->second;
					if(ind[2] >= round)
						continue;
					auto tuple2 = table[2][ind[2]];
					if(cond_func[2] != nullptr && (*cond_func[2])(tuple2) == false)
						continue;
					auto range3 = h_table[3][2].equal_range((*key_func[2][3])(tuple2));
					for(auto iter3 = range3.first; iter3 != range3.second; ++iter3){
						ind[3] = iter3->second;
						if(ind[3] >= round)
							continue;
						auto tuple3 = table[3][ind[3]];
						if(cond_func[3] != nullptr && (*cond_func[3])(tuple3) == false)
							continue;
						for(size_t k = 0; k < K; ++k){
							x[k] = x[k] - q[k][ind[k]]/n/n;
							w[k] = w[k] - (x[k] - double(n-1)/n/n*q[k][ind[k]])/n*(x[k]-double(n-1)/n/n*q[k][ind[k]])/(n-1);
							q[k][ind[k]] = q[k][ind[k]] + sample_values[ind[3]]*total_size;
							w[k] = w[k] + (x[k] - double(n-1)/n/n*q[k][ind[k]])/n*(x[k]-double(n-1)/n/n*q[k][ind[k]])/(n-1);
							x[k] = x[k] + q[k][ind[k]]/n/n;
						}
					sum += sample_values[ind[3]];
					test_count++;
					}
				}
			}
			}
		}
		//ik == n from table 1
		{
		std::vector<size_t> ind(K);
		ind[1] = round;
		auto tuple1 = table[1][ind[1]];
		if(cond_func[1] == nullptr || (*cond_func[1])(tuple1)){
			auto range0 = h_table[0][1].equal_range((*key_func[1][0])(tuple1));
			for(auto iter0 = range0.first; iter0 != range0.second; iter0++){
				ind[0] = iter0->second;
				if(ind[0] >= round)
					continue;
				auto tuple0 = table[0][ind[0]];
				if(cond_func[0] != nullptr && !(*cond_func[0])(tuple0))
					continue;
				auto range2 = h_table[2][1].equal_range((*key_func[1][2])(tuple1));
				for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
					ind[2] = iter2->second;
					if(ind[2] >= round)
						continue;
					auto tuple2 = table[2][ind[2]];
					if(cond_func[2] != nullptr && (*cond_func[2])(tuple2) == false)
						continue;
					auto range3 = h_table[3][2].equal_range((*key_func[2][3])(tuple2));
					for(auto iter3 = range3.first; iter3 != range3.second; ++iter3){
						ind[3] = iter3->second;
						if(ind[3] >= round)
							continue;
						auto tuple3 = table[3][ind[3]];
						if(cond_func[3] != nullptr && (*cond_func[3])(tuple3) == false)
							continue;
						for(size_t k = 0; k < K; ++k){
							x[k] = x[k] - q[k][ind[k]]/n/n;
							w[k] = w[k] - (x[k] - double(n-1)/n/n*q[k][ind[k]])/n*(x[k]-double(n-1)/n/n*q[k][ind[k]])/(n-1);
							q[k][ind[k]] = q[k][ind[k]] + sample_values[ind[3]]*total_size;
							w[k] = w[k] + (x[k] - double(n-1)/n/n*q[k][ind[k]])/n*(x[k]-double(n-1)/n/n*q[k][ind[k]])/(n-1);
							x[k] = x[k] + q[k][ind[k]]/n/n;
						}
					sum += sample_values[ind[3]];
					test_count++;
					}
				}
			}
		}
		}
		//ik == n from table 2
		{
		std::vector<size_t> ind(K);
		ind[2] = round;
		auto tuple2 = table[2][ind[2]];
		if(cond_func[2] == nullptr || (*cond_func[2])(tuple2)){
			auto range3 = h_table[3][2].equal_range((*key_func[2][3])(tuple2));
			for(auto iter3 = range3.first; iter3 != range3.second; ++iter3){
				ind[3] = iter3->second;
				if(ind[3] >= round)
					continue;
				auto tuple3 = table[3][ind[3]];
				if(cond_func[3] != nullptr && (*cond_func[3])(tuple3) == false)
					continue;
				auto range1 = h_table[1][2].equal_range((*key_func[2][1])(tuple2));
				for(auto iter1 = range1.first; iter1 != range1.second; iter1++){
					ind[1] = iter1->second;
					if(ind[1] >= round)
						continue;
					auto tuple1 = table[1][ind[1]];
					if(cond_func[1] != nullptr && !(*cond_func[1])(tuple1))
						continue;
					auto range0 = h_table[0][1].equal_range((*key_func[1][0])(tuple1));
					for(auto iter0 = range0.first; iter0 != range0.second; ++iter0){
						ind[0] = iter0->second;
						if(ind[0] >= round)
							continue;
						auto tuple0 = table[0][ind[0]];
						if(cond_func[0] != nullptr && !(*cond_func[0])(tuple0))
							continue;
						for(size_t k = 0; k < K; ++k){
							x[k] = x[k] - q[k][ind[k]]/n/n;
							w[k] = w[k] - (x[k] - double(n-1)/n/n*q[k][ind[k]])/n*(x[k]-double(n-1)/n/n*q[k][ind[k]])/(n-1);
							q[k][ind[k]] = q[k][ind[k]] + sample_values[ind[3]]*total_size;
							w[k] = w[k] + (x[k] - double(n-1)/n/n*q[k][ind[k]])/n*(x[k]-double(n-1)/n/n*q[k][ind[k]])/(n-1);
							x[k] = x[k] + q[k][ind[k]]/n/n;
						}
						sum += sample_values[ind[3]];
						test_count++;
					}
				}
			}
		}
		}
		//! i_n from table 3
		{
		std::vector<size_t> ind(K);
		ind[3] = round;
		auto tuple3 = table[3][ind[3]];
		if(cond_func[3] == nullptr || (*cond_func[3])(tuple3)){
			auto range2 = h_table[2][3].equal_range((*key_func[3][2])(tuple3));
			for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
				ind[2] = iter2->second;
				if(ind[2] >= round)
					continue;
				auto tuple2 = table[2][ind[2]];
				if(cond_func[2] != nullptr && (*cond_func[2])(tuple2) == false)
					continue;
				auto range1 = h_table[1][2].equal_range((*key_func[2][1])(tuple2));
				for(auto iter1 = range1.first; iter1 != range1.second; iter1++){
					ind[1] = iter1->second;
					if(ind[1] >= round)
						continue;
					auto tuple1 = table[1][ind[1]];
					if(cond_func[1] != nullptr && !(*cond_func[1])(tuple1))
						continue;
					auto range0 = h_table[0][1].equal_range((*key_func[1][0])(tuple1));
					for(auto iter0 = range0.first; iter0 != range0.second; ++iter0){
						ind[0] = iter0->second;
						if(ind[0] >= round)
							continue;
						auto tuple0 = table[0][ind[0]];
						if(cond_func[0] != nullptr && !(*cond_func[0])(tuple0))
							continue;
						for(size_t k = 0; k < K; ++k){
							x[k] = x[k] - q[k][ind[k]]/n/n;
							w[k] = w[k] - (x[k] - double(n-1)/n/n*q[k][ind[k]])/n*(x[k]-double(n-1)/n/n*q[k][ind[k]])/(n-1);
							q[k][ind[k]] = q[k][ind[k]] + sample_values[ind[3]]*total_size;
							w[k] = w[k] + (x[k] - double(n-1)/n/n*q[k][ind[k]])/n*(x[k]-double(n-1)/n/n*q[k][ind[k]])/(n-1);
							x[k] = x[k] + q[k][ind[k]]/n/n;
						}
						sum += sample_values[ind[3]];
						test_count++;
					}
				}
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
			tool::report(STEP*report_round, n, test_count, 0, 0, sum*total_size/nk, ci, prob);
			steptimer.restart();
			report_round++;
		}
		round++;
	}
}
void ripple_6(std::vector<std::vector<base_raw *> > &table,\
		std::vector<std::vector<key_func_type>> &key_func,\
		std::vector<cond_func_type> &cond_func,\
		std::function<double(base_raw *)> result_func,
		const double STEP,\
		const size_t &MAX,\
		const double &prob,\
		bool verbose,\
		double total_size){
	constexpr int K = 6;
	size_t max_round = 0;
	for(size_t i = 0; i < K; ++i){
		if(i == 0)
			max_round = table[i].size();
		else
			max_round = std::min(max_round, table[i].size());
	}
	std::vector<double> sample_values;

	tool::start_report();
	size_t test_count = 0;
	typedef unordered_multimap<uint32_t, size_t> index_type;
	std::vector<std::vector<index_type> > h_table;
	for(size_t i = 0; i < K; ++i){
		std::vector<index_type> tmp;
		h_table.push_back(tmp);
		for(size_t j = 0; j < K; ++j){
			index_type x;
			h_table[i].push_back(x);
		}
	}
	std::vector<double> q[K]; 
	std::vector<double> w(K, 0);
	std::vector<double> x(K, 0);
	for (size_t k = 0; k < K; ++k) q[k].resize(max_round);

	cout<<total_size<<endl;
	const double zp = tool::xql_erf_inv(prob);

	double y = 0.0;
	double ci = 1.0 / 0.0;
	size_t round = 0;
	/* round 1 */
	for(size_t i = 0; i < K; ++i){
		for(size_t j = 0; j < K; ++j){
			if(key_func[i][j] != nullptr)
				h_table[i][j].insert(std::make_pair((*key_func[i][j])(table[i][round]), round)); 
		}
	}
	double sum = 0;
	double v = result_func(table[2][round]);
	if(v > 0)
		test_count++;
	sample_values.push_back(v);
	if((*key_func[0][1])(table[0][round]) == (*key_func[1][0])(table[1][round]) \
			&& (*key_func[1][2])(table[1][round]) == (*key_func[2][1])(table[2][round]) \
			&& (*key_func[2][3])(table[2][round])== (*key_func[3][2])(table[3][round]) \
			&& (*key_func[4][3])(table[4][round]) == (*key_func[3][4])(table[3][round]) \
			&& (*key_func[5][4])(table[5][round]) == (*key_func[4][5])(table[4][round]))

	{
		y = v *total_size;
		sum+=v;
	}
	for (size_t k = 0; k < K; ++k) {
		q[k][0] = v;
		w[k] = 0;
		x[k] = y;
	}
	size_t n = 1;
	double z = 0;
	uint32_t n_rejected = 0;
	size_t report_round = 1;
	timer maxtimer;
	timer steptimer;
	round = 1;
	while(maxtimer.get_elapsed() < MAX && round < max_round){
		n = round + 1;
		for(size_t i = 0; i < K; ++i){
			for(size_t j = 0; j < K; ++j){
				if(key_func[i][j] != nullptr)
						h_table[i][j].insert(std::make_pair((*key_func[i][j])(table[i][round]), round)); 
			}
		}
		sample_values.push_back(result_func(table[2][round]));	
		y = y * std::pow(double(round)/(round+1), K-1);
		v = sample_values.back();
		if((*key_func[0][1])(table[0][round]) == (*key_func[1][0])(table[1][round]) \
				&& (*key_func[1][2])(table[1][round]) == (*key_func[2][1])(table[2][round]) \
				&& (*key_func[2][3])(table[2][round])== (*key_func[3][2])(table[3][round]) \
				&& (*key_func[4][3])(table[4][round]) == (*key_func[3][4])(table[3][round]) \
				&& (*key_func[5][4])(table[5][round]) == (*key_func[4][5])(table[4][round]))

		{
			y = v *total_size;
			sum+=v;
		}

		for(size_t k = 0; k < K; ++k){
			q[k][round] = 0;
			w[k] = w[k] * std::pow(double(round) / (round+1), 6) + y * y /(round+1)/round;
			x[k] = y;
		}
		//ik == n from table 0
		{
		std::vector<size_t> ind(K);
		ind[0] = round;
		auto tuple0 = table[0][ind[0]];
		if(cond_func[0] == nullptr || cond_func[0] != nullptr && (*cond_func[0])(tuple0)){
			auto range1 = h_table[1][0].equal_range((*key_func[0][1])(tuple0));
			for(auto iter1 = range1.first; iter1 != range1.second; iter1++){
				ind[1] = iter1->second;
				if(ind[1] >= round)
					continue;
				auto tuple1 = table[1][ind[1]];
				if(cond_func[1] != nullptr && (*cond_func[1])(tuple1) == false)
					continue;
				auto range2 = h_table[2][1].equal_range((*key_func[1][2])(tuple1));
				for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
					ind[2] = iter2->second;
					if(ind[2] >= round)
						continue;
					auto tuple2 = table[2][ind[2]];
					if(cond_func[2] != nullptr && (*cond_func[2])(tuple2) == false)
						continue;
					auto range3 = h_table[3][2].equal_range((*key_func[2][3])(tuple2));
					for(auto iter3 = range3.first; iter3 != range3.second; ++iter3){
						ind[3] = iter3->second;
						if(ind[3] >= round)
							continue;
						auto tuple3 = table[3][ind[3]];
						if(cond_func[3] != nullptr && (*cond_func[3])(tuple3) == false)
							continue;
						auto range4 = h_table[4][3].equal_range((*key_func[3][4])(tuple3));
						for(auto iter4 = range4.first; iter4 != range4.second; ++iter4){
							ind[4] = iter4->second;
							if(ind[4] >= round)
								continue;
							auto tuple4 = table[4][ind[4]];
							if(cond_func[4] != nullptr && (*cond_func[4])(tuple4) == false)
								continue;
							auto range5 = h_table[5][4].equal_range((*key_func[4][5])(tuple4));
							for(auto iter5 = range5.first; iter5 != range5.second; ++iter5){
								ind[5] = iter5->second;
								if(ind[5] >= round)
									continue;
								auto tuple5 = table[5][ind[5]];
								if(cond_func[5] != nullptr && (*cond_func[5])(tuple5) == false)
									continue;
								for(size_t k = 0; k < K; ++k){
									x[k] = x[k] - q[k][ind[k]]/n/n;
									w[k] = w[k] - (x[k] - double(n-1)/n/n*q[k][ind[k]])/n*(x[k]-double(n-1)/n/n*q[k][ind[k]])/(n-1);
									q[k][ind[k]] = q[k][ind[k]] + sample_values[ind[2]]*total_size;
									w[k] = w[k] + (x[k] - double(n-1)/n/n*q[k][ind[k]])/n*(x[k]-double(n-1)/n/n*q[k][ind[k]])/(n-1);
									x[k] = x[k] + q[k][ind[k]]/n/n;
								}
								sum += sample_values[ind[2]];
								test_count++;
							}
						}
					}
				}
			}
		}
		}
		//ik == n from table 1
		{
		std::vector<size_t> ind(K);
		ind[1] = round;
		auto tuple1 = table[1][ind[1]];
		if(cond_func[1] == nullptr || (*cond_func[1])(tuple1)){
			auto range0 = h_table[0][1].equal_range((*key_func[1][0])(tuple1));
			for(auto iter0 = range0.first; iter0 != range0.second; iter0++){
				ind[0] = iter0->second;
				if(ind[0] >= round)
					continue;
				auto tuple0 = table[0][ind[0]];
				if(cond_func[0] != nullptr && !(*cond_func[0])(tuple0))
					continue;
				auto range2 = h_table[2][1].equal_range((*key_func[1][2])(tuple1));
				for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
					ind[2] = iter2->second;
					if(ind[2] >= round)
						continue;
					auto tuple2 = table[2][ind[2]];
					if(cond_func[2] != nullptr && (*cond_func[2])(tuple2) == false)
						continue;
					auto range3 = h_table[3][2].equal_range((*key_func[2][3])(tuple2));
					for(auto iter3 = range3.first; iter3 != range3.second; ++iter3){
						ind[3] = iter3->second;
						if(ind[3] >= round)
							continue;
						auto tuple3 = table[3][ind[3]];
						if(cond_func[3] != nullptr && (*cond_func[3])(tuple3) == false)
							continue;
						auto range4 = h_table[4][3].equal_range((*key_func[3][4])(tuple3));
						for(auto iter4 = range4.first; iter4 != range4.second; ++iter4){
							ind[4] = iter4->second;
							if(ind[4] >= round)
								continue;
							auto tuple4 = table[4][ind[4]];
							if(cond_func[4] != nullptr && (*cond_func[4])(tuple4) == false)
								continue;
							auto range5 = h_table[5][4].equal_range((*key_func[4][5])(tuple4));
							for(auto iter5 = range5.first; iter5 != range5.second; ++iter5){
								ind[5] = iter5->second;
								if(ind[5] >= round)
									continue;
								auto tuple5 = table[5][ind[5]];
								if(cond_func[5] != nullptr && (*cond_func[5])(tuple5) == false)
									continue;
								for(size_t k = 0; k < K; ++k){
									x[k] = x[k] - q[k][ind[k]]/n/n;
									w[k] = w[k] - (x[k] - double(n-1)/n/n*q[k][ind[k]])/n*(x[k]-double(n-1)/n/n*q[k][ind[k]])/(n-1);
									q[k][ind[k]] = q[k][ind[k]] + sample_values[ind[2]]*total_size;
									w[k] = w[k] + (x[k] - double(n-1)/n/n*q[k][ind[k]])/n*(x[k]-double(n-1)/n/n*q[k][ind[k]])/(n-1);
									x[k] = x[k] + q[k][ind[k]]/n/n;
								}
								sum += sample_values[ind[2]];
								test_count++;
							}
						}
					}
				}
			}
		}
		}
		//ik == n from table 2
		{
		std::vector<size_t> ind(K);
		ind[2] = round;
		auto tuple2 = table[2][ind[2]];
		if(cond_func[2] == nullptr || (*cond_func[2])(tuple2)){
			auto range1 = h_table[1][2].equal_range((*key_func[2][1])(tuple2));
			for(auto iter1 = range1.first; iter1 != range1.second; iter1++){
				ind[1] = iter1->second;
				if(ind[1] >= round)
					continue;
				auto tuple1 = table[1][ind[1]];
				if(cond_func[1] != nullptr && !(*cond_func[1])(tuple1))
					continue;
				auto range0 = h_table[0][1].equal_range((*key_func[1][0])(tuple1));
				for(auto iter0 = range0.first; iter0 != range0.second; ++iter0){
					ind[0] = iter0->second;
					if(ind[0] >= round)
						continue;
					auto tuple0 = table[0][ind[0]];
					if(cond_func[0] != nullptr && (*cond_func[0])(tuple0) == false)
						continue;
					auto range3 = h_table[3][2].equal_range((*key_func[2][3])(tuple2));
					for(auto iter3 = range3.first; iter3 != range3.second; ++iter3){
						ind[3] = iter3->second;
						if(ind[3] >= round)
							continue;
						auto tuple3 = table[3][ind[3]];
						if(cond_func[3] != nullptr && (*cond_func[3])(tuple3) == false)
							continue;
						auto range4 = h_table[4][3].equal_range((*key_func[3][4])(tuple3));
						for(auto iter4 = range4.first; iter4 != range4.second; ++iter4){
							ind[4] = iter4->second;
							if(ind[4] >= round)
								continue;
							auto tuple4 = table[4][ind[4]];
							if(cond_func[4] != nullptr && (*cond_func[4])(tuple4) == false)
								continue;
							auto range5 = h_table[5][4].equal_range((*key_func[4][5])(tuple4));
							for(auto iter5 = range5.first; iter5 != range5.second; ++iter5){
								ind[5] = iter5->second;
								if(ind[5] >= round)
									continue;
								auto tuple5 = table[5][ind[5]];
								if(cond_func[5] != nullptr && (*cond_func[5])(tuple5) == false)
									continue;
								for(size_t k = 0; k < K; ++k){
									x[k] = x[k] - q[k][ind[k]]/n/n;
									w[k] = w[k] - (x[k] - double(n-1)/n/n*q[k][ind[k]])/n*(x[k]-double(n-1)/n/n*q[k][ind[k]])/(n-1);
									q[k][ind[k]] = q[k][ind[k]] + sample_values[ind[2]]*total_size;
									w[k] = w[k] + (x[k] - double(n-1)/n/n*q[k][ind[k]])/n*(x[k]-double(n-1)/n/n*q[k][ind[k]])/(n-1);
									x[k] = x[k] + q[k][ind[k]]/n/n;
								}
								sum += sample_values[ind[3]];
								test_count++;
							}
						}
					}
				}
			}
		}
		}
		//! i_n from table 3
		{
		std::vector<size_t> ind(K);
		ind[3] = round;
		auto tuple3 = table[3][ind[3]];
		if(cond_func[3] == nullptr || (*cond_func[3])(tuple3)){
			auto range2 = h_table[2][3].equal_range((*key_func[3][2])(tuple3));
			for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
				ind[2] = iter2->second;
				if(ind[2] >= round)
					continue;
				auto tuple2 = table[2][ind[2]];
				if(cond_func[2] != nullptr && (*cond_func[2])(tuple2) == false)
					continue;
				auto range1 = h_table[1][2].equal_range((*key_func[2][1])(tuple2));
				for(auto iter1 = range1.first; iter1 != range1.second; iter1++){
					ind[1] = iter1->second;
					if(ind[1] >= round)
						continue;
					auto tuple1 = table[1][ind[1]];
					if(cond_func[1] != nullptr && !(*cond_func[1])(tuple1))
						continue;
					auto range0 = h_table[0][1].equal_range((*key_func[1][0])(tuple1));
					for(auto iter0 = range0.first; iter0 != range0.second; ++iter0){
						ind[0] = iter0->second;
						if(ind[0] >= round)
							continue;
						auto tuple0 = table[0][ind[0]];
						if(cond_func[0] != nullptr && !(*cond_func[0])(tuple0))
							continue;
						auto range4 = h_table[4][3].equal_range((*key_func[3][4])(tuple3));
						for(auto iter4 = range4.first; iter4 != range4.second; ++iter4){
							ind[4] = iter4->second;
							if(ind[4] >= round)
								continue;
							auto tuple4 = table[4][ind[4]];
							if(cond_func[4] != nullptr && (*cond_func[4])(tuple4) == false)
								continue;
							auto range5 = h_table[5][4].equal_range((*key_func[4][5])(tuple4));
							for(auto iter5 = range5.first; iter5 != range5.second; ++iter5){
								ind[5] = iter5->second;
								if(ind[5] >= round)
									continue;
								auto tuple5 = table[5][ind[5]];
								if(cond_func[5] != nullptr && (*cond_func[5])(tuple5) == false)
									continue;
								for(size_t k = 0; k < K; ++k){
									x[k] = x[k] - q[k][ind[k]]/n/n;
									w[k] = w[k] - (x[k] - double(n-1)/n/n*q[k][ind[k]])/n*(x[k]-double(n-1)/n/n*q[k][ind[k]])/(n-1);
									q[k][ind[k]] = q[k][ind[k]] + sample_values[ind[2]]*total_size;
									w[k] = w[k] + (x[k] - double(n-1)/n/n*q[k][ind[k]])/n*(x[k]-double(n-1)/n/n*q[k][ind[k]])/(n-1);
									x[k] = x[k] + q[k][ind[k]]/n/n;
								}
								sum += sample_values[ind[2]];
								test_count++;
							}
						}
					}
				}
			}
		}
		}
		//ik == n from table 4
		{
		std::vector<size_t> ind(K);
		ind[4] = round;
		auto tuple4 = table[4][ind[4]];
		if(cond_func[4] == nullptr || cond_func[4] != nullptr && (*cond_func[4])(tuple4)){
			auto range3 = h_table[3][4].equal_range((*key_func[4][3])(tuple4));
			for(auto iter3 = range3.first; iter3 != range3.second; iter3++){
				ind[3] = iter3->second;
				if(ind[3] >= round)
					continue;
				auto tuple3 = table[3][ind[3]];
				if(cond_func[3] != nullptr && (*cond_func[3])(tuple3) == false)
					continue;
				auto range2 = h_table[2][3].equal_range((*key_func[3][2])(tuple3));
				for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
					ind[2] = iter2->second;
					if(ind[2] >= round)
						continue;
					auto tuple2 = table[2][ind[2]];
					if(cond_func[2] != nullptr && (*cond_func[2])(tuple2) == false)
						continue;
					auto range1 = h_table[1][2].equal_range((*key_func[2][1])(tuple2));
					for(auto iter1 = range1.first; iter1 != range1.second; ++iter1){
						ind[1] = iter1->second;
						if(ind[1] >= round)
							continue;
						auto tuple1 = table[1][ind[1]];
						if(cond_func[1] != nullptr && (*cond_func[1])(tuple1) == false)
							continue;
						auto range0 = h_table[0][1].equal_range((*key_func[1][0])(tuple1));
						for(auto iter0 = range0.first; iter0 != range0.second; ++iter0){
							ind[0] = iter0->second;
							if(ind[0] >= round)
								continue;
							auto tuple0 = table[0][ind[0]];
							if(cond_func[0] != nullptr && (*cond_func[0])(tuple0) == false)
								continue;
							auto range5 = h_table[5][4].equal_range((*key_func[4][5])(tuple4));
							for(auto iter5 = range5.first; iter5 != range5.second; ++iter5){
								ind[5] = iter5->second;
								if(ind[5] >= round)
									continue;
								auto tuple5 = table[5][ind[5]];
								if(cond_func[5] != nullptr && (*cond_func[5])(tuple5) == false)
									continue;
								for(size_t k = 0; k < K; ++k){
									x[k] = x[k] - q[k][ind[k]]/n/n;
									w[k] = w[k] - (x[k] - double(n-1)/n/n*q[k][ind[k]])/n*(x[k]-double(n-1)/n/n*q[k][ind[k]])/(n-1);
									q[k][ind[k]] = q[k][ind[k]] + sample_values[ind[2]]*total_size;
									w[k] = w[k] + (x[k] - double(n-1)/n/n*q[k][ind[k]])/n*(x[k]-double(n-1)/n/n*q[k][ind[k]])/(n-1);
									x[k] = x[k] + q[k][ind[k]]/n/n;
								}
								sum += sample_values[ind[2]];
								test_count++;
							}
						}
					}
				}
			}
		}
		}
		//ik == n from table 5
		{
		std::vector<size_t> ind(K);
		ind[5] = round;
		auto tuple5 = table[5][ind[5]];
		if(cond_func[5] == nullptr || cond_func[5] != nullptr && (*cond_func[5])(tuple5)){
			auto range4 = h_table[4][5].equal_range((*key_func[5][4])(tuple5));
			for(auto iter4 = range4.first; iter4 != range4.second; iter4++){
				ind[4] = iter4->second;
				if(ind[4] >= round)
					continue;
				auto tuple4 = table[4][ind[4]];
				auto range3 = h_table[3][4].equal_range((*key_func[4][3])(tuple4));
				for(auto iter3 = range3.first; iter3 != range3.second; iter3++){
					ind[3] = iter3->second;
					if(ind[3] >= round)
						continue;
					auto tuple3 = table[3][ind[3]];
					if(cond_func[3] != nullptr && (*cond_func[3])(tuple3) == false)
						continue;
					auto range2 = h_table[2][3].equal_range((*key_func[3][2])(tuple3));
					for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
						ind[2] = iter2->second;
						if(ind[2] >= round)
							continue;
						auto tuple2 = table[2][ind[2]];
						if(cond_func[2] != nullptr && (*cond_func[2])(tuple2) == false)
							continue;
						auto range1 = h_table[1][2].equal_range((*key_func[2][1])(tuple2));
						for(auto iter1 = range1.first; iter1 != range1.second; ++iter1){
							ind[1] = iter1->second;
							if(ind[1] >= round)
								continue;
							auto tuple1 = table[1][ind[1]];
							if(cond_func[1] != nullptr && (*cond_func[1])(tuple1) == false)
								continue;
							auto range0 = h_table[0][1].equal_range((*key_func[1][0])(tuple1));
							for(auto iter0 = range0.first; iter0 != range0.second; ++iter0){
								ind[0] = iter0->second;
								if(ind[0] >= round)
									continue;
								auto tuple0 = table[0][ind[0]];
								if(cond_func[0] != nullptr && (*cond_func[0])(tuple0) == false)
									continue;
									for(size_t k = 0; k < K; ++k){
										x[k] = x[k] - q[k][ind[k]]/n/n;
										w[k] = w[k] - (x[k] - double(n-1)/n/n*q[k][ind[k]])/n*(x[k]-double(n-1)/n/n*q[k][ind[k]])/(n-1);
										q[k][ind[k]] = q[k][ind[k]] + sample_values[ind[2]]*total_size;
										w[k] = w[k] + (x[k] - double(n-1)/n/n*q[k][ind[k]])/n*(x[k]-double(n-1)/n/n*q[k][ind[k]])/(n-1);
										x[k] = x[k] + q[k][ind[k]]/n/n;
									}
									sum += sample_values[ind[2]];
									test_count++;
							}
						}
					}
				}
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
			tool::report(STEP*report_round, n, test_count, 0, 0, sum*total_size/nk, ci, prob);
			steptimer.restart();
			report_round++;
		}
		round++;
		}
		
}
void ripple_4_ratio(std::vector<std::vector<base_raw *> > &table,\
		std::vector<size_t> &a,\
		std::vector<std::vector<key_func_type>> &key_func,\
		std::vector<cond_func_type> &cond_func,\
		std::function<double(base_raw*)> result_func,\
		const double STEP,\
		const size_t &MAX,\
		const double &prob,\
		bool verbose,\
		double total_size){
	constexpr int K = 4;
	size_t max_round = 0;
	for(size_t i = 0; i < K; ++i){
		if(i == 0)
			max_round = table[i].size();
		else
			max_round = std::min(max_round, table[i].size());
	}
	std::vector<double> sample_values;

	tool::start_report();
	size_t test_count = 0;
	typedef unordered_multimap<uint32_t, size_t> index_type;
	std::vector<std::vector<index_type> > h_table;
	for(size_t i = 0; i < K; ++i){
		std::vector<index_type> tmp;
		h_table.push_back(tmp);
		for(size_t j = 0; j < K; ++j){
			index_type x;
			h_table[i].push_back(x);
		}
	}
	std::vector<double> q[K]; 
	std::vector<double> w(K, 0);
	std::vector<double> x(K, 0);
	for (size_t k = 0; k < K; ++k) q[k].resize(max_round);

	cout<<total_size<<endl;
	const double zp = tool::xql_erf_inv(prob);

	double y = 0.0;
	double ci = 1.0 / 0.0;
	size_t round = 0;
	/* round 1 */
	for(size_t i = 0; i < K; ++i){
		for(size_t k = 0; k < a[i]; ++k){
			for(size_t j = 0; j < K; ++j){
					if(key_func[i][j] != nullptr)
						if(cond_func[i] == nullptr || (*cond_func[i])(table[i][round*a[i] +k]))
							h_table[i][j].insert(std::make_pair((*key_func[i][j])(table[i][round * a[i] + k]), round*a[i]+k)); 
			}
		}
	}
	double sum = 0;
	double v = 0;
	for(size_t i = 0; i < a[3]; ++i)
			sample_values.push_back(result_func(table[3][round*a[3] + i]));	
	for(size_t i0 = 0; i0 < a[0]; ++i0){
		for(size_t i1 = 0; i1 < a[1]; ++i1){
			for(size_t i2 = 0; i2 < a[2]; ++i2){
				for(size_t i3 = 0; i3 < a[3]; ++i3){
					std::vector<size_t> ind(K);
					ind[0] = round*a[0] + i0;
					ind[1] = round*a[1] + i1;
					ind[2] = round*a[2] + i2;
					ind[3] = round*a[3] + i3;
				if((*key_func[0][1])(table[0][ind[0]]) == (*key_func[1][0])(table[1][ind[1]]) \
						&& (*key_func[1][2])(table[1][ind[1]]) == (*key_func[2][1])(table[2][ind[2]]) \
						&& (*key_func[2][3])(table[2][ind[2]])== (*key_func[3][2])(table[3][ind[3]])){
						if((cond_func[0] == nullptr || (*cond_func[0])(table[0][ind[0]])) &&\
								(cond_func[1] == nullptr || (*cond_func[1])(table[1][ind[1]])) &&\
								(cond_func[2] == nullptr || (*cond_func[2])(table[2][ind[2]])) &&\
								 (cond_func[3]==nullptr || (*cond_func[3])(table[3][ind[3]]))){
							v = sample_values[round*a[3] + i3];
							y += v *total_size/a[0]/a[1]/a[2]/a[3];
							sum+=v;
							test_count++;
						}
				}
				}
			}
		}
	}
	for (size_t k = 0; k < K; ++k) {
		q[k][0] = 0;
		w[k] = 0;
		x[k] = y;
	}
	size_t n = 1;
	double z = 0;
	uint32_t n_rejected = 0;
	size_t report_round = 1;
	timer maxtimer;
	timer steptimer;
	round = 1;
	while(maxtimer.get_elapsed() < MAX && round < max_round){
		n = round + 1;
		for(size_t i = 0; i < K; ++i){
			for(size_t k = 0; k < a[i]; ++k){
				for(size_t j = 0; j < K; ++j){
					if(key_func[i][j] != nullptr){
						if(cond_func[i] == nullptr || (*cond_func[i])(table[i][round*a[i] +k] ))
							h_table[i][j].insert(std::make_pair((*key_func[i][j])(table[i][round * a[i] + k]), round*a[i]+k)); 
					}
				}
			}
		}
		for(size_t i = 0; i < a[3]; ++i)
			sample_values.push_back(result_func(table[3][round*a[3] + i]));	
		y = y * std::pow(double(round)/(round+1), K-1);
		for(size_t i0 = 0; i0 < a[0]; ++i0){
			for(size_t i1 = 0; i1 < a[1]; ++i1){
				for(size_t i2 = 0; i2 < a[2]; ++i2){
					for(size_t i3 = 0; i3 < a[3]; ++i3){
						std::vector<size_t> ind(K);
						ind[0] = round*a[0] + i0;
						ind[1] = round*a[1] + i1;
						ind[2] = round*a[2] + i2;
						ind[3] = round*a[3] + i3;
					if((*key_func[0][1])(table[0][ind[0]]) == (*key_func[1][0])(table[1][ind[1]]) \
							&& (*key_func[1][2])(table[1][ind[1]]) == (*key_func[2][1])(table[2][ind[2]]) \
							&& (*key_func[2][3])(table[2][ind[2]])== (*key_func[3][2])(table[3][ind[3]])){
							if((cond_func[0] == nullptr || (*cond_func[0])(table[0][ind[0]])) &&\
									(cond_func[1] == nullptr || (*cond_func[1])(table[1][ind[1]])) &&\
									(cond_func[2] == nullptr || (*cond_func[2])(table[2][ind[2]])) &&\
									 (cond_func[3]==nullptr || (*cond_func[3])(table[3][ind[3]]))){
								v = sample_values[round*a[3] + i3];
								y += v *total_size/a[0]/a[1]/a[2]/a[3];
								sum+=v;
								test_count++;
							}
					}
					}
				}
			}
		}

		for(size_t k = 0; k < K; ++k){
			q[k][0] = 0;
			w[k] = w[k] * std::pow(double(round) / (round+1), 2*(K-1)) + y * y /(round+1)/round;
			x[k] = y;
		}
		//ik == n from table 0
		for(size_t i = 0; i < a[0]; ++i)
		{
		std::vector<size_t> ind(K);
		ind[0] = round*a[0] + i;
		auto tuple0 = table[0][ind[0]];
		if(cond_func[0] == nullptr || cond_func[0] != nullptr && (*cond_func[0])(tuple0)){
			auto range1 = h_table[1][0].equal_range((*key_func[0][1])(tuple0));
			for(auto iter1 = range1.first; iter1 != range1.second; iter1++){
				ind[1] = iter1->second;
				if(ind[1]/a[1] >= round)
					continue;
				auto tuple1 = table[1][ind[1]];
				if(cond_func[1] != nullptr && (*cond_func[1])(tuple1) == false)
					continue;
				auto range2 = h_table[2][1].equal_range((*key_func[1][2])(tuple1));
				for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
					ind[2] = iter2->second;
					if(ind[2]/a[2] >= round)
						continue;
					auto tuple2 = table[2][ind[2]];
					if(cond_func[2] != nullptr && (*cond_func[2])(tuple2) == false)
						continue;
					auto range3 = h_table[3][2].equal_range((*key_func[2][3])(tuple2));
					for(auto iter3 = range3.first; iter3 != range3.second; ++iter3){
						ind[3] = iter3->second;
						if(ind[3]/a[3] >= round)
							continue;
						auto tuple3 = table[3][ind[3]];
						if(cond_func[3] != nullptr && (*cond_func[3])(tuple3) == false)
							continue;
						calc_ripple_ci(K, ind, x, w, q, sample_values, n, a, total_size);
						sum += sample_values[ind[3]];
						test_count++;
					}
				}
			}
			}
		}
		//ik == n from table 1
		for(size_t i = 0; i < a[1]; ++i)
		{
		std::vector<size_t> ind(K);
		ind[1] = round * a[1] + i;
		auto tuple1 = table[1][ind[1]];
		if(cond_func[1] == nullptr || (*cond_func[1])(tuple1)){
			auto range0 = h_table[0][1].equal_range((*key_func[1][0])(tuple1));
			for(auto iter0 = range0.first; iter0 != range0.second; iter0++){
				ind[0] = iter0->second;
				if(ind[0]/a[0]>= round)
					continue;
				auto tuple0 = table[0][ind[0]];
				if(cond_func[0] != nullptr && !(*cond_func[0])(tuple0))
					continue;
				auto range2 = h_table[2][1].equal_range((*key_func[1][2])(tuple1));
				for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
					ind[2] = iter2->second;
					if(ind[2]/a[2]>= round)
						continue;
					auto tuple2 = table[2][ind[2]];
					if(cond_func[2] != nullptr && (*cond_func[2])(tuple2) == false)
						continue;
					auto range3 = h_table[3][2].equal_range((*key_func[2][3])(tuple2));
					for(auto iter3 = range3.first; iter3 != range3.second; ++iter3){
						ind[3] = iter3->second;
						if(ind[3]/a[3] >= round)
							continue;
						auto tuple3 = table[3][ind[3]];
						if(cond_func[3] != nullptr && (*cond_func[3])(tuple3) == false)
							continue;
						calc_ripple_ci(K, ind, x, w, q, sample_values, n, a, total_size);
						sum += sample_values[ind[3]];
						test_count++;
					}
				}
			}
		}
		}
		//ik == n from table 2
		for(size_t i = 0; i < a[2]; ++i)
		{
		std::vector<size_t> ind(K);
		ind[2] = a[2]*round + i;
		auto tuple2 = table[2][ind[2]];
		if(cond_func[2] == nullptr || (*cond_func[2])(tuple2)){
			auto range3 = h_table[3][2].equal_range((*key_func[2][3])(tuple2));
			for(auto iter3 = range3.first; iter3 != range3.second; ++iter3){
				ind[3] = iter3->second;
				if(ind[3]/a[3]>= round)
					continue;
				auto tuple3 = table[3][ind[3]];
				if(cond_func[3] != nullptr && (*cond_func[3])(tuple3) == false)
					continue;
				auto range1 = h_table[1][2].equal_range((*key_func[2][1])(tuple2));
				for(auto iter1 = range1.first; iter1 != range1.second; iter1++){
					ind[1] = iter1->second;
					if(ind[1]/a[1] >= round)
						continue;
					auto tuple1 = table[1][ind[1]];
					if(cond_func[1] != nullptr && !(*cond_func[1])(tuple1))
						continue;
					auto range0 = h_table[0][1].equal_range((*key_func[1][0])(tuple1));
					for(auto iter0 = range0.first; iter0 != range0.second; ++iter0){
						ind[0] = iter0->second;
						if(ind[0]/a[0] >= round)
							continue;
						auto tuple0 = table[0][ind[0]];
						if(cond_func[0] != nullptr && !(*cond_func[0])(tuple0))
							continue;
						calc_ripple_ci(K, ind, x, w, q, sample_values, n, a, total_size);
						sum += sample_values[ind[3]];
						test_count++;
					}
				}
			}
		}
		}
		//! i_n from table 3
		for(size_t i = 0; i < a[3]; ++i)
		{
		std::vector<size_t> ind(K);
		ind[3] = round*a[3] + i;
		auto tuple3 = table[3][ind[3]];
		if(cond_func[3] == nullptr || (*cond_func[3])(tuple3)){
			auto range2 = h_table[2][3].equal_range((*key_func[3][2])(tuple3));
			for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
				ind[2] = iter2->second;
				if(ind[2]/a[2]>= round)
					continue;
				auto tuple2 = table[2][ind[2]];
				if(cond_func[2] != nullptr && (*cond_func[2])(tuple2) == false)
					continue;
				auto range1 = h_table[1][2].equal_range((*key_func[2][1])(tuple2));
				for(auto iter1 = range1.first; iter1 != range1.second; iter1++){
					ind[1] = iter1->second;
					if(ind[1]/a[1]>= round)
						continue;
					auto tuple1 = table[1][ind[1]];
					if(cond_func[1] != nullptr && !(*cond_func[1])(tuple1))
						continue;
					auto range0 = h_table[0][1].equal_range((*key_func[1][0])(tuple1));
					for(auto iter0 = range0.first; iter0 != range0.second; ++iter0){
						ind[0] = iter0->second;
						if(ind[0]/a[0] >= round)
							continue;
						auto tuple0 = table[0][ind[0]];
						if(cond_func[0] != nullptr && !(*cond_func[0])(tuple0))
							continue;
						calc_ripple_ci(K, ind, x, w, q, sample_values, n, a, total_size);
						sum += sample_values[ind[3]];
						test_count++;
					}
				}
			}
		}
		}
		y = x[0];
		z = 0;
		double nk=std::pow(n, K);
		for(size_t i = 0; i < K; ++i)
			z += w[i] / (n * a[i]-1);
		if(steptimer.get_elapsed()> STEP){
			ci = zp*std::sqrt(z) / sqrt(n);
			tool::report(STEP*report_round, n, test_count, 0, 0, sum*total_size/nk, ci, prob);
			steptimer.restart();
			report_round++;
		}
		round++;
	}
}
void ripple_4_ratio(std::vector<std::vector<base_raw *> > &table,\
		std::vector<size_t> &a,\
		std::vector<std::vector<key_func_type>> &key_func,\
		std::vector<cond_func_type> &cond_func,\
		std::function<double(base_raw*)> result_func,\
		const double STEP,\
		const size_t &MAX,\
		const double &prob,\
		bool verbose,\
		double total_size,\
		size_t col_id){
	constexpr int K = 4;
	size_t max_round = 0;
	for(size_t i = 0; i < K; ++i){
		if(i == 0)
			max_round = table[i].size();
		else
			max_round = std::min(max_round, table[i].size());
	}
	std::vector<double> sample_values;

	tool::start_report();
	size_t test_count = 0;
	typedef unordered_multimap<uint32_t, size_t> index_type;
	std::vector<std::vector<index_type> > h_table;
	for(size_t i = 0; i < K; ++i){
		std::vector<index_type> tmp;
		h_table.push_back(tmp);
		for(size_t j = 0; j < K; ++j){
			index_type x;
			h_table[i].push_back(x);
		}
	}
	std::vector<double> q[K]; 
	std::vector<double> w(K, 0);
	std::vector<double> x(K, 0);
	for (size_t k = 0; k < K; ++k) q[k].resize(max_round);

	cout<<total_size<<endl;
	const double zp = tool::xql_erf_inv(prob);

	double y = 0.0;
	double ci = 1.0 / 0.0;
	size_t round = 0;
	/* round 1 */
	for(size_t i = 0; i < K; ++i){
		for(size_t k = 0; k < a[i]; ++k){
			for(size_t j = 0; j < K; ++j){
					if(key_func[i][j] != nullptr)
						if(cond_func[i] == nullptr || (*cond_func[i])(table[i][round*a[i] +k]))
							h_table[i][j].insert(std::make_pair((*key_func[i][j])(table[i][round * a[i] + k]), round*a[i]+k)); 
			}
		}
	}
	double sum = 0;
	double v = 0;
	for(size_t i = 0; i < a[3]; ++i)
			sample_values.push_back(result_func(table[2][round*a[2] + i]));	
	for(size_t i0 = 0; i0 < a[0]; ++i0){
		for(size_t i1 = 0; i1 < a[1]; ++i1){
			for(size_t i2 = 0; i2 < a[2]; ++i2){
				for(size_t i3 = 0; i3 < a[3]; ++i3){
					std::vector<size_t> ind(K);
					ind[0] = round*a[0] + i0;
					ind[1] = round*a[1] + i1;
					ind[2] = round*a[2] + i2;
					ind[3] = round*a[3] + i3;
				if((*key_func[0][1])(table[0][ind[0]]) == (*key_func[1][0])(table[1][ind[1]]) \
						&& (*key_func[1][2])(table[1][ind[1]]) == (*key_func[2][1])(table[2][ind[2]]) \
						&& (*key_func[2][3])(table[2][ind[2]])== (*key_func[3][2])(table[3][ind[3]])){
						if((cond_func[0] == nullptr || (*cond_func[0])(table[0][ind[0]])) &&\
								(cond_func[1] == nullptr || (*cond_func[1])(table[1][ind[1]])) &&\
								(cond_func[2] == nullptr || (*cond_func[2])(table[2][ind[2]])) &&\
								 (cond_func[3]==nullptr || (*cond_func[3])(table[3][ind[3]]))){
							v = sample_values[round*a[2] + i2];
							y += v *total_size/a[0]/a[1]/a[2]/a[3];
							sum+=v;
							test_count++;
						}
				}
				}
			}
		}
	}
	for (size_t k = 0; k < K; ++k) {
		q[k][0] = 0;
		w[k] = 0;
		x[k] = y;
	}
	size_t n = 1;
	double z = 0;
	uint32_t n_rejected = 0;
	size_t report_round = 1;
	timer maxtimer;
	timer steptimer;
	round = 1;
	while(maxtimer.get_elapsed() < MAX && round < max_round){
		n = round + 1;
		for(size_t i = 0; i < K; ++i){
			for(size_t k = 0; k < a[i]; ++k){
				for(size_t j = 0; j < K; ++j){
					if(key_func[i][j] != nullptr){
						if(cond_func[i] == nullptr || (*cond_func[i])(table[i][round*a[i] +k] ))
							h_table[i][j].insert(std::make_pair((*key_func[i][j])(table[i][round * a[i] + k]), round*a[i]+k)); 
					}
				}
			}
		}
		for(size_t i = 0; i < a[col_id]; ++i)
			sample_values.push_back(result_func(table[col_id][round*a[col_id] + i]));	
		y = y * std::pow(double(round)/(round+1), K-1);
		for(size_t i0 = 0; i0 < a[0]; ++i0){
			for(size_t i1 = 0; i1 < a[1]; ++i1){
				for(size_t i2 = 0; i2 < a[2]; ++i2){
					for(size_t i3 = 0; i3 < a[3]; ++i3){
						std::vector<size_t> ind(K);
						ind[0] = round*a[0] + i0;
						ind[1] = round*a[1] + i1;
						ind[2] = round*a[2] + i2;
						ind[3] = round*a[3] + i3;
					if((*key_func[0][1])(table[0][ind[0]]) == (*key_func[1][0])(table[1][ind[1]]) \
							&& (*key_func[1][2])(table[1][ind[1]]) == (*key_func[2][1])(table[2][ind[2]]) \
							&& (*key_func[2][3])(table[2][ind[2]])== (*key_func[3][2])(table[3][ind[3]])){
							if((cond_func[0] == nullptr || (*cond_func[0])(table[0][ind[0]])) &&\
									(cond_func[1] == nullptr || (*cond_func[1])(table[1][ind[1]])) &&\
									(cond_func[2] == nullptr || (*cond_func[2])(table[2][ind[2]])) &&\
									 (cond_func[3]==nullptr || (*cond_func[3])(table[3][ind[3]]))){
								v = sample_values[round*a[2] + i2];
								y += v *total_size/a[0]/a[1]/a[2]/a[3];
								sum+=v;
								test_count++;
							}
					}
					}
				}
			}
		}

		for(size_t k = 0; k < K; ++k){
			q[k][round] = 0;
			w[k] = w[k] * std::pow(double(round) / (round+1), 2*(K-1)) + y * y /(round+1)/round;
			x[k] = y;
		}
		//ik == n from table 0
		for(size_t i = 0; i < a[0]; ++i)
		{
		std::vector<size_t> ind(K);
		ind[0] = round*a[0] + i;
		auto tuple0 = table[0][ind[0]];
		if(cond_func[0] == nullptr || cond_func[0] != nullptr && (*cond_func[0])(tuple0)){
			auto range1 = h_table[1][0].equal_range((*key_func[0][1])(tuple0));
			for(auto iter1 = range1.first; iter1 != range1.second; iter1++){
				ind[1] = iter1->second;
				if(ind[1]/a[1] >= round)
					continue;
				auto tuple1 = table[1][ind[1]];
				if(cond_func[1] != nullptr && (*cond_func[1])(tuple1) == false)
					continue;
				auto range2 = h_table[2][1].equal_range((*key_func[1][2])(tuple1));
				for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
					ind[2] = iter2->second;
					if(ind[2]/a[2] >= round)
						continue;
					auto tuple2 = table[2][ind[2]];
					if(cond_func[2] != nullptr && (*cond_func[2])(tuple2) == false)
						continue;
					auto range3 = h_table[3][2].equal_range((*key_func[2][3])(tuple2));
					for(auto iter3 = range3.first; iter3 != range3.second; ++iter3){
						ind[3] = iter3->second;
						if(ind[3]/a[3] >= round)
							continue;
						auto tuple3 = table[3][ind[3]];
						if(cond_func[3] != nullptr && (*cond_func[3])(tuple3) == false)
							continue;
						calc_ripple_ci(K, ind, x, w, q, sample_values, n, a, total_size);
						sum += sample_values[ind[2]];
						test_count++;
					}
				}
			}
			}
		}
		//ik == n from table 1
		for(size_t i = 0; i < a[1]; ++i)
		{
		std::vector<size_t> ind(K);
		ind[1] = round * a[1] + i;
		auto tuple1 = table[1][ind[1]];
		if(cond_func[1] == nullptr || (*cond_func[1])(tuple1)){
			auto range0 = h_table[0][1].equal_range((*key_func[1][0])(tuple1));
			for(auto iter0 = range0.first; iter0 != range0.second; iter0++){
				ind[0] = iter0->second;
				if(ind[0]/a[0]>= round)
					continue;
				auto tuple0 = table[0][ind[0]];
				if(cond_func[0] != nullptr && !(*cond_func[0])(tuple0))
					continue;
				auto range2 = h_table[2][1].equal_range((*key_func[1][2])(tuple1));
				for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
					ind[2] = iter2->second;
					if(ind[2]/a[2]>= round)
						continue;
					auto tuple2 = table[2][ind[2]];
					if(cond_func[2] != nullptr && (*cond_func[2])(tuple2) == false)
						continue;
					auto range3 = h_table[3][2].equal_range((*key_func[2][3])(tuple2));
					for(auto iter3 = range3.first; iter3 != range3.second; ++iter3){
						ind[3] = iter3->second;
						if(ind[3]/a[3] >= round)
							continue;
						auto tuple3 = table[3][ind[3]];
						if(cond_func[3] != nullptr && (*cond_func[3])(tuple3) == false)
							continue;
						calc_ripple_ci(K, ind, x, w, q, sample_values, n, a, total_size);
						sum += sample_values[ind[2]];
						test_count++;
					}
				}
			}
		}
		}
		//ik == n from table 2
		for(size_t i = 0; i < a[2]; ++i)
		{
		std::vector<size_t> ind(K);
		ind[2] = a[2]*round + i;
		auto tuple2 = table[2][ind[2]];
		if(cond_func[2] == nullptr || (*cond_func[2])(tuple2)){
			auto range3 = h_table[3][2].equal_range((*key_func[2][3])(tuple2));
			for(auto iter3 = range3.first; iter3 != range3.second; ++iter3){
				ind[3] = iter3->second;
				if(ind[3]/a[3]>= round)
					continue;
				auto tuple3 = table[3][ind[3]];
				if(cond_func[3] != nullptr && (*cond_func[3])(tuple3) == false)
					continue;
				auto range1 = h_table[1][2].equal_range((*key_func[2][1])(tuple2));
				for(auto iter1 = range1.first; iter1 != range1.second; iter1++){
					ind[1] = iter1->second;
					if(ind[1]/a[1] >= round)
						continue;
					auto tuple1 = table[1][ind[1]];
					if(cond_func[1] != nullptr && !(*cond_func[1])(tuple1))
						continue;
					auto range0 = h_table[0][1].equal_range((*key_func[1][0])(tuple1));
					for(auto iter0 = range0.first; iter0 != range0.second; ++iter0){
						ind[0] = iter0->second;
						if(ind[0]/a[0] >= round)
							continue;
						auto tuple0 = table[0][ind[0]];
						if(cond_func[0] != nullptr && !(*cond_func[0])(tuple0))
							continue;
						calc_ripple_ci(K, ind, x, w, q, sample_values, n, a, total_size);
						sum += sample_values[ind[2]];
						test_count++;
					}
				}
			}
		}
		}
		//! i_n from table 3
		for(size_t i = 0; i < a[3]; ++i)
		{
		std::vector<size_t> ind(K);
		ind[3] = round*a[3] + i;
		auto tuple3 = table[3][ind[3]];
		if(cond_func[3] == nullptr || (*cond_func[3])(tuple3)){
			auto range2 = h_table[2][3].equal_range((*key_func[3][2])(tuple3));
			for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
				ind[2] = iter2->second;
				if(ind[2]/a[2]>= round)
					continue;
				auto tuple2 = table[2][ind[2]];
				if(cond_func[2] != nullptr && (*cond_func[2])(tuple2) == false)
					continue;
				auto range1 = h_table[1][2].equal_range((*key_func[2][1])(tuple2));
				for(auto iter1 = range1.first; iter1 != range1.second; iter1++){
					ind[1] = iter1->second;
					if(ind[1]/a[1]>= round)
						continue;
					auto tuple1 = table[1][ind[1]];
					if(cond_func[1] != nullptr && !(*cond_func[1])(tuple1))
						continue;
					auto range0 = h_table[0][1].equal_range((*key_func[1][0])(tuple1));
					for(auto iter0 = range0.first; iter0 != range0.second; ++iter0){
						ind[0] = iter0->second;
						if(ind[0]/a[0] >= round)
							continue;
						auto tuple0 = table[0][ind[0]];
						if(cond_func[0] != nullptr && !(*cond_func[0])(tuple0))
							continue;
						calc_ripple_ci(K, ind, x, w, q, sample_values, n, a, total_size);
						sum += sample_values[ind[2]];
						test_count++;
					}
				}
			}
		}
		}
		z = 0;
		double nk=std::pow(n, K);
		y = sum*total_size/nk;
		for(size_t i = 0; i < K; ++i)
			z += w[i] / (n * a[i]-1);
		if(steptimer.get_elapsed()> STEP){
			ci = zp*std::sqrt(z) / sqrt(n);
			tool::report(STEP*report_round, n, test_count, 0, 0, sum*total_size/nk, ci, prob);
			steptimer.restart();
			report_round++;
		}
		round++;
	}
}
void ripple_4_ratio_test(std::vector<std::vector<base_raw *> > &table,\
		std::vector<size_t> &a,\
		std::vector<std::vector<key_func_type>> &key_func,\
		std::vector<cond_func_type> &cond_func,\
		std::function<double(base_raw*)> result_func,\
		const double STEP,\
		const size_t &MAX,\
		const double &prob,\
		bool verbose,\
		double &total_size,\
		const double &VALUE,\
		size_t col_id){
	constexpr int K = 4;
	size_t max_round = 0;
	for(size_t i = 0; i < K; ++i){
		if(i == 0)
			max_round = table[i].size();
		else
			max_round = std::min(max_round, table[i].size());
	}
	std::vector<double> sample_values;
	total_size = 1;
	for(size_t i = 0; i < K; ++i)
		total_size *= table[i].size();

	tool::start_report();
	size_t test_count = 0;
	typedef unordered_multimap<uint32_t, size_t> index_type;
	std::vector<std::vector<index_type> > h_table;
	for(size_t i = 0; i < K; ++i){
		std::vector<index_type> tmp;
		h_table.push_back(tmp);
		for(size_t j = 0; j < K; ++j){
			index_type x;
			h_table[i].push_back(x);
		}
	}
	std::vector<double> q[K]; 
	std::vector<double> w(K, 0);
	std::vector<double> x(K, 0);
	for (size_t k = 0; k < K; ++k) q[k].resize(max_round * a[k]);
	for (size_t k = 0; k < K; ++k)
		for(size_t i = 0; i < q[k].size(); ++i)
			q[k][i] = 0;

	cout<<total_size<<endl;

	const double zp = tool::xql_erf_inv(prob);

	double y = 0.0;
	double ci = 1.0 / 0.0;
	size_t round = 0;
	/* round 0 */
	for(size_t i = 0; i < K; ++i){
		for(size_t k = 0; k < a[i]; ++k){
			for(size_t j = 0; j < K; ++j){
					if(key_func[i][j] != nullptr)
						if(cond_func[i] == nullptr || (*cond_func[i])(table[i][round*a[i] +k]))
							h_table[i][j].insert(std::make_pair((*key_func[i][j])(table[i][round * a[i] + k]), round*a[i]+k)); 
			}
		}
	}
	double sum = 0;
	double v = 0;
	for(size_t i = 0; i < a[col_id]; ++i)
			sample_values.push_back(result_func(table[col_id][round*a[col_id] + i]));	
	size_t i[4];
	for(i[0] = 0; i[0] < a[0]; ++i[0]){
		for(i[1] = 0; i[1] < a[1]; ++i[1]){
			for(i[2] = 0; i[2] < a[2]; ++i[2]){
				for(i[3] = 0; i[3] < a[3]; ++i[3]){
					std::vector<size_t> ind(K);
					ind[0] = round*a[0] + i[0];
					ind[1] = round*a[1] + i[1];
					ind[2] = round*a[2] + i[2];
					ind[3] = round*a[3] + i[3];
				if((*key_func[0][1])(table[0][ind[0]]) == (*key_func[1][0])(table[1][ind[1]]) \
						&& (*key_func[1][2])(table[1][ind[1]]) == (*key_func[2][1])(table[2][ind[2]]) \
						&& (*key_func[2][3])(table[2][ind[2]])== (*key_func[3][2])(table[3][ind[3]])){
						if((cond_func[0] == nullptr || (*cond_func[0])(table[0][ind[0]])) &&\
								(cond_func[1] == nullptr || (*cond_func[1])(table[1][ind[1]])) &&\
								(cond_func[2] == nullptr || (*cond_func[2])(table[2][ind[2]])) &&\
								 (cond_func[3]==nullptr || (*cond_func[3])(table[3][ind[3]]))){
							v = sample_values[round*a[col_id] + i[col_id]];
							y = v * total_size;
							q[0][ind[0]] += y;
							q[1][ind[1]] += y;
							q[2][ind[2]] += y;
							q[3][ind[3]] += y;
							sum+=v;
							test_count++;
						}
				}
				}
			}
		}
	}
	for(size_t i = 0; i < K; ++i){
		y = sum/a[i] * total_size;
	}
	double pow_n[K];
	double all = 1;
	for(size_t i = 0; i < K; ++i)
		all *= (round+1) * a[i];
	for(size_t i = 0; i < K; ++i)
		pow_n[i] = all/(round+1)/a[i];
	for(size_t i = 0; i < K; ++i){
		for(size_t j = 0; j < a[i]; ++j){
			x[i] += q[i][round*a[i]+ j]/a[i];
		}
	}
	for(size_t i = 0; i < K; ++i){
		w[i] = 0;
	}
	size_t n = 1;
	double z = 0;
	uint32_t n_rejected = 0;
	size_t report_round = 1;
	timer maxtimer;
	timer steptimer;
	size_t old_time = 0;
	round = 1;
	while(maxtimer.get_elapsed() < MAX && round < max_round){
		n = round + 1;
		for(size_t i = 0; i < K; ++i)
			all *= (round+1) * a[i];
		for(size_t i = 0; i < K; ++i)
			pow_n[i] = all/(round+1)/a[i];
		for(size_t i = 0; i < K; ++i){
			for(size_t k = 0; k < a[i]; ++k){
				for(size_t j = 0; j < K; ++j){
					if(key_func[i][j] != nullptr){
						if(cond_func[i] == nullptr || (*cond_func[i])(table[i][round*a[i] +k] ))
							h_table[i][j].insert(std::make_pair((*key_func[i][j])(table[i][round * a[i] + k]), round*a[i]+k)); 
					}
				}
			}
		}
		for(size_t i = 0; i < a[col_id]; ++i)
			sample_values.push_back(result_func(table[col_id][round*a[col_id] + i]));	
		y = y * std::pow(double(round)/(round+1), K-1);
		for(size_t k = 0; k < K; ++k){
			for(size_t i = 0; i < a[k]; ++i)
				q[k][round*a[k] + i ] = 0;
		}
		for(i[0] = 0; i[0] < a[0]; ++i[0]){
			for(i[1] = 0; i[1] < a[1]; ++i[1]){
				for(i[2] = 0; i[2] < a[2]; ++i[2]){
					for(i[3] = 0; i[3] < a[3]; ++i[3]){
						std::vector<size_t> ind(K);
						ind[0] = round*a[0] + i[0];
						ind[1] = round*a[1] + i[1];
						ind[2] = round*a[2] + i[2];
						ind[3] = round*a[3] + i[3];
					if((*key_func[0][1])(table[0][ind[0]]) == (*key_func[1][0])(table[1][ind[1]]) \
							&& (*key_func[1][2])(table[1][ind[1]]) == (*key_func[2][1])(table[2][ind[2]]) \
							&& (*key_func[2][3])(table[2][ind[2]])== (*key_func[3][2])(table[3][ind[3]])){
							if((cond_func[0] == nullptr || (*cond_func[0])(table[0][ind[0]])) &&\
									(cond_func[1] == nullptr || (*cond_func[1])(table[1][ind[1]])) &&\
									(cond_func[2] == nullptr || (*cond_func[2])(table[2][ind[2]])) &&\
									 (cond_func[3]==nullptr || (*cond_func[3])(table[3][ind[3]]))){
								v = sample_values[round*a[col_id] + i[col_id]];
								auto tmp = v *total_size;
								q[0][ind[0]] += tmp;
								q[1][ind[1]] += tmp;
								q[2][ind[2]] += tmp;
								q[3][ind[3]] += tmp;
								sum += v;
								test_count++;
							}
					}
					}
				}
			}
		}

		for(size_t i = 0; i < K; ++i){
			auto tmp = 0;
			for(size_t j = 0; j < a[i]; ++j){
				tmp += q[i][round*a[i] + j];
			}
			x[i] = y + tmp/a[i]/ pow_n[i];
		}

		for(size_t k = 0; k < K; ++k){
			for(size_t i = 0; i < a[k]; ++i)
				q[k][round*a[k] + i ] = 0;
			w[k] = w[k] * std::pow(double(round) / (round+1), 6) + y * y /(round+1)/round;
		}
		//ik == n from table 0
		for(size_t i = 0; i < a[0]; ++i)
		{
		std::vector<size_t> ind(K);
		ind[0] = round*a[0] + i;
		auto tuple0 = table[0][ind[0]];
		if(cond_func[0] == nullptr || cond_func[0] != nullptr && (*cond_func[0])(tuple0)){
			auto range1 = h_table[1][0].equal_range((*key_func[0][1])(tuple0));
			for(auto iter1 = range1.first; iter1 != range1.second; iter1++){
				ind[1] = iter1->second;
				if(ind[1]/a[1] >= round)
					continue;
				auto tuple1 = table[1][ind[1]];
				auto range2 = h_table[2][1].equal_range((*key_func[1][2])(tuple1));
				for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
					ind[2] = iter2->second;
					if(ind[2]/a[2] >= round)
						continue;
					auto tuple2 = table[2][ind[2]];
					auto range3 = h_table[3][2].equal_range((*key_func[2][3])(tuple2));
					for(auto iter3 = range3.first; iter3 != range3.second; ++iter3){
						ind[3] = iter3->second;
						if(ind[3]/a[3] >= round)
							continue;
						auto tuple3 = table[3][ind[3]];
						calc_ripple_ci_test(K, ind, x, w, q, sample_values, n, a, total_size, col_id);
						sum += sample_values[ind[col_id]];
						test_count++;
					}
				}
			}
			}
		}
		//ik == n from table 1
		for(size_t i = 0; i < a[1]; ++i)
		{
		std::vector<size_t> ind(K);
		ind[1] = round * a[1] + i;
		auto tuple1 = table[1][ind[1]];
		if(cond_func[1] == nullptr || (*cond_func[1])(tuple1)){
			auto range0 = h_table[0][1].equal_range((*key_func[1][0])(tuple1));
			for(auto iter0 = range0.first; iter0 != range0.second; iter0++){
				ind[0] = iter0->second;
				if(ind[0]/a[0]>= round)
					continue;
				auto tuple0 = table[0][ind[0]];
				auto range2 = h_table[2][1].equal_range((*key_func[1][2])(tuple1));
				for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
					ind[2] = iter2->second;
					if(ind[2]/a[2]>= round)
						continue;
					auto tuple2 = table[2][ind[2]];
					auto range3 = h_table[3][2].equal_range((*key_func[2][3])(tuple2));
					for(auto iter3 = range3.first; iter3 != range3.second; ++iter3){
						ind[3] = iter3->second;
						if(ind[3]/a[3] >= round)
							continue;
						auto tuple3 = table[3][ind[3]];
						calc_ripple_ci_test(K, ind, x, w, q, sample_values, n, a, total_size, col_id);
						sum += sample_values[ind[col_id]];
						test_count++;
					}
				}
			}
		}
		}
		//ik == n from table 2
		for(size_t i = 0; i < a[2]; ++i)
		{
		std::vector<size_t> ind(K);
		ind[2] = a[2]*round + i;
		auto tuple2 = table[2][ind[2]];
		if(cond_func[2] == nullptr || (*cond_func[2])(tuple2)){
			auto range3 = h_table[3][2].equal_range((*key_func[2][3])(tuple2));
			for(auto iter3 = range3.first; iter3 != range3.second; ++iter3){
				ind[3] = iter3->second;
				if(ind[3]/a[3]>= round)
					continue;
				auto tuple3 = table[3][ind[3]];
				auto range1 = h_table[1][2].equal_range((*key_func[2][1])(tuple2));
				for(auto iter1 = range1.first; iter1 != range1.second; iter1++){
					ind[1] = iter1->second;
					if(ind[1]/a[1] >= round)
						continue;
					auto tuple1 = table[1][ind[1]];
					auto range0 = h_table[0][1].equal_range((*key_func[1][0])(tuple1));
					for(auto iter0 = range0.first; iter0 != range0.second; ++iter0){
						ind[0] = iter0->second;
						if(ind[0]/a[0] >= round)
							continue;
						auto tuple0 = table[0][ind[0]];
						calc_ripple_ci_test(K, ind, x, w, q, sample_values, n, a, total_size, col_id);
						sum += sample_values[ind[col_id]];
						test_count++;
					}
				}
			}
		}
		}
		//! i_n from table 3
		for(size_t i = 0; i < a[3]; ++i)
		{
		std::vector<size_t> ind(K);
		ind[3] = round*a[3] + i;
		auto tuple3 = table[3][ind[3]];
		if(cond_func[3] == nullptr || (*cond_func[3])(tuple3)){
			auto range2 = h_table[2][3].equal_range((*key_func[3][2])(tuple3));
			for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
				ind[2] = iter2->second;
				if(ind[2]/a[2]>= round)
					continue;
				auto tuple2 = table[2][ind[2]];
				auto range1 = h_table[1][2].equal_range((*key_func[2][1])(tuple2));
				for(auto iter1 = range1.first; iter1 != range1.second; iter1++){
					ind[1] = iter1->second;
					if(ind[1]/a[1]>= round)
						continue;
					auto tuple1 = table[1][ind[1]];
					auto range0 = h_table[0][1].equal_range((*key_func[1][0])(tuple1));
					for(auto iter0 = range0.first; iter0 != range0.second; ++iter0){
						ind[0] = iter0->second;
						if(ind[0]/a[0] >= round)
							continue;
						auto tuple0 = table[0][ind[0]];
						calc_ripple_ci_test(K, ind, x, w, q, sample_values, n, a, total_size, col_id);
						sum += sample_values[ind[col_id]];
						test_count++;
					}
				}
			}
		}
		}
		y = x[0];
		z = 0;
		double nk = 1;
		for(size_t i = 0; i < K; ++i)
			nk *= a[i] * n;
		for(size_t i = 0; i < K; ++i)
			z += w[i] / (n * a[i]-1);
		auto cur_time = maxtimer.get_elapsed();
		ci = zp*std::sqrt(z) / sqrt(n);
		if(ci > 0 && ci < VALUE)
			return;
		if(verbose && cur_time - old_time> STEP){
			tool::report(cur_time, n, test_count, 0, 0, sum/nk*total_size, ci, prob);
			old_time = cur_time;
			report_round++;
		}
		round++;
	}
	return;
}
void ripple_4_ratio_test(std::vector<std::vector<base_raw *> > &table,\
		std::vector<size_t> &a,\
		std::vector<std::vector<key_func_type>> &key_func,\
		std::vector<cond_func_type> &cond_func,\
		std::function<double(base_raw*)> result_func,\
		const double STEP,\
		const size_t &MAX,\
		const double &prob,\
		bool verbose,\
		double &total_size,\
		size_t col_id){
	ripple_4_ratio_test(table, a, key_func, cond_func, result_func, STEP, MAX, prob, verbose,\
			total_size, -1, col_id);
}
void ripple_6_ratio(std::vector<std::vector<base_raw *> > &table,\
		std::vector<size_t> &a,\
		std::vector<std::vector<key_func_type>> &key_func,\
		std::vector<cond_func_type> &cond_func,\
		std::function<double(base_raw*)> result_func,\
		const double STEP,\
		const size_t &MAX,\
		const double &prob,\
		bool verbose,\
		double total_size){
	constexpr int K = 6;
	size_t max_round = 0;
	for(size_t i = 0; i < K; ++i){
		if(i == 0)
			max_round = table[i].size();
		else
			max_round = std::min(max_round, table[i].size());
	}
	std::vector<double> sample_values;

	tool::start_report();
	size_t test_count = 0;
	typedef unordered_multimap<uint32_t, size_t> index_type;
	std::vector<std::vector<index_type> > h_table;
	for(size_t i = 0; i < K; ++i){
		std::vector<index_type> tmp;
		h_table.push_back(tmp);
		for(size_t j = 0; j < K; ++j){
			index_type x;
			h_table[i].push_back(x);
		}
	}
	std::vector<double> q[K]; 
	std::vector<double> w(K, 0);
	std::vector<double> x(K, 0);
	for (size_t k = 0; k < K; ++k) q[k].resize(max_round);

	cout<<total_size<<endl;
	const double zp = tool::xql_erf_inv(prob);

	double y = 0.0;
	double ci = 1.0 / 0.0;
	size_t round = 0;
	/* round 1 */
	for(size_t i = 0; i < K; ++i){
		for(size_t k = 0; k < a[i]; ++k){
			for(size_t j = 0; j < K; ++j){
					if(key_func[i][j] != nullptr)
						if(cond_func[i] == nullptr || (*cond_func[i])(table[i][round*a[i] +k]))
							h_table[i][j].insert(std::make_pair((*key_func[i][j])(table[i][round * a[i] + k]), round*a[i]+k)); 
			}
		}
	}
	double sum = 0;
	double v = 0;
	for(size_t i = 0; i < a[2]; ++i)
			sample_values.push_back(result_func(table[2][round*a[2] + i]));	
	for(size_t i0 = 0; i0 < a[0]; ++i0){
		for(size_t i1 = 0; i1 < a[1]; ++i1){
			for(size_t i2 = 0; i2 < a[2]; ++i2){
				for(size_t i3 = 0; i3 < a[3]; ++i3){
					for(size_t i4 = 0; i4< a[4]; ++i4){
						for(size_t i5 = 0; i5 < a[5]; ++i5){
					std::vector<size_t> ind(K);
					ind[0] = round*a[0] + i0;
					ind[1] = round*a[1] + i1;
					ind[2] = round*a[2] + i2;
					ind[3] = round*a[3] + i3;
					ind[4] = round*a[4] + i4;
					ind[5] = round*a[5] + i5;
						if((*key_func[0][1])(table[0][ind[0]]) == (*key_func[1][0])(table[1][ind[1]]) \
								&& (*key_func[1][2])(table[1][ind[1]]) == (*key_func[2][1])(table[2][ind[2]]) \
								&& (*key_func[2][3])(table[2][ind[2]])== (*key_func[3][2])(table[3][ind[3]])\
								&& (*key_func[3][4])(table[3][ind[3]]) == (*key_func[4][3])(table[4][ind[4]])\
								&& (*key_func[4][5])(table[4][ind[4]]) == (*key_func[5][4])(table[5][ind[5]])){
								if((cond_func[0] == nullptr || (*cond_func[0])(table[0][ind[0]])) &&\
										(cond_func[1] == nullptr || (*cond_func[1])(table[1][ind[1]])) &&\
										(cond_func[2] == nullptr || (*cond_func[2])(table[2][ind[2]])) &&\
										 (cond_func[3]==nullptr || (*cond_func[3])(table[3][ind[3]])) &&\
										 (cond_func[4] == nullptr || (*cond_func[4])(table[4][ind[4]])) &&\
										 (cond_func[5] == nullptr || (*cond_func[5])(table[5][ind[5]]))){
									v = sample_values[round*a[3] + i3];
									y += v *total_size/a[0]/a[1]/a[2]/a[3]/a[4]/a[5];
									sum+=v;
									test_count++;
								}
						}
						}
					}
				}
			}
		}
	}
	for (size_t k = 0; k < K; ++k) {
		q[k][0] = v;
		w[k] = 0;
		x[k] = y;
	}
	size_t n = 1;
	double z = 0;
	uint32_t n_rejected = 0;
	size_t report_round = 1;
	timer maxtimer;
	timer steptimer;
	round = 1;
	while(maxtimer.get_elapsed() < MAX && round < max_round){
		n = round + 1;
		for(size_t i = 0; i < K; ++i){
			for(size_t k = 0; k < a[i]; ++k){
				for(size_t j = 0; j < K; ++j){
						if(key_func[i][j] != nullptr)
							if(cond_func[i] == nullptr || (*cond_func[i])(table[i][round*a[i] +k]))
								h_table[i][j].insert(std::make_pair((*key_func[i][j])(table[i][round * a[i] + k]), round*a[i]+k)); 
				}
			}
		}
		for(size_t i = 0; i < a[2]; ++i)
				sample_values.push_back(result_func(table[2][round*a[2] + i]));	
		y = y * std::pow(double(round)/(round+1), K-1);
		if((*key_func[0][1])(table[0][round]) == (*key_func[1][0])(table[1][round]) \
				&& (*key_func[1][2])(table[1][round]) == (*key_func[2][1])(table[2][round]) \
				&& (*key_func[2][3])(table[2][round])== (*key_func[3][2])(table[3][round]) \
				&& (*key_func[4][3])(table[4][round]) == (*key_func[3][4])(table[3][round]) \
				&& (*key_func[5][4])(table[5][round]) == (*key_func[4][5])(table[4][round]))

		{
			y = v *total_size;
			sum+=v;
		}

		for(size_t k = 0; k < K; ++k){
			q[k][round] = 0;
			w[k] = w[k] * std::pow(double(round) / (round+1), 6) + y * y /(round+1)/round;
			x[k] = y;
		}
		//ik == n from table 0
		for(size_t i = 0; i < a[0]; ++i)
		{
		std::vector<size_t> ind(K);
		ind[0] = round*a[0] + i;
		auto tuple0 = table[0][ind[0]];
		if(cond_func[0] == nullptr || cond_func[0] != nullptr && (*cond_func[0])(tuple0)){
			auto range1 = h_table[1][0].equal_range((*key_func[0][1])(tuple0));
			for(auto iter1 = range1.first; iter1 != range1.second; iter1++){
				ind[1] = iter1->second;
				if(ind[1]/a[1]>= round)
					continue;
				auto tuple1 = table[1][ind[1]];
				auto range2 = h_table[2][1].equal_range((*key_func[1][2])(tuple1));
				for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
					ind[2] = iter2->second;
					if(ind[2]/a[2] >= round)
						continue;
					auto tuple2 = table[2][ind[2]];
					auto range3 = h_table[3][2].equal_range((*key_func[2][3])(tuple2));
					for(auto iter3 = range3.first; iter3 != range3.second; ++iter3){
						ind[3] = iter3->second;
						if(ind[3]/a[3] >= round)
							continue;
						auto tuple3 = table[3][ind[3]];
						auto range4 = h_table[4][3].equal_range((*key_func[3][4])(tuple3));
						for(auto iter4 = range4.first; iter4 != range4.second; ++iter4){
							ind[4] = iter4->second;
							if(ind[4]/a[4] >= round)
								continue;
							auto tuple4 = table[4][ind[4]];
							auto range5 = h_table[5][4].equal_range((*key_func[4][5])(tuple4));
							for(auto iter5 = range5.first; iter5 != range5.second; ++iter5){
								ind[5] = iter5->second;
								if(ind[5] >= round)
									continue;
								auto tuple5 = table[5][ind[5]];
								calc_ripple_ci(K, ind, x, w, q, sample_values, n, a, total_size);
								sum += sample_values[ind[2]];
								test_count++;
							}
						}
					}
				}
			}
		}
		}
		//ik == n from table 1
		for(size_t i = 0; i < a[1]; ++i)
		{
		std::vector<size_t> ind(K);
		ind[1] = round * a[1] + i;
		auto tuple1 = table[1][ind[1]];
		if(cond_func[1] == nullptr || (*cond_func[1])(tuple1)){
			auto range0 = h_table[0][1].equal_range((*key_func[1][0])(tuple1));
			for(auto iter0 = range0.first; iter0 != range0.second; iter0++){
				ind[0] = iter0->second;
				if(ind[0]/a[0]>= round)
					continue;
				auto tuple0 = table[0][ind[0]];
				auto range2 = h_table[2][1].equal_range((*key_func[1][2])(tuple1));
				for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
					ind[2] = iter2->second;
					if(ind[2] / a[2]>= round)
						continue;
					auto tuple2 = table[2][ind[2]];
					auto range3 = h_table[3][2].equal_range((*key_func[2][3])(tuple2));
					for(auto iter3 = range3.first; iter3 != range3.second; ++iter3){
						ind[3] = iter3->second;
						if(ind[3] /a[3]>= round)
							continue;
						auto tuple3 = table[3][ind[3]];
						auto range4 = h_table[4][3].equal_range((*key_func[3][4])(tuple3));
						for(auto iter4 = range4.first; iter4 != range4.second; ++iter4){
							ind[4] = iter4->second;
							if(ind[4] / a[4]>= round)
								continue;
							auto tuple4 = table[4][ind[4]];
							auto range5 = h_table[5][4].equal_range((*key_func[4][5])(tuple4));
							for(auto iter5 = range5.first; iter5 != range5.second; ++iter5){
								ind[5] = iter5->second;
								if(ind[5]/a[5]>= round)
									continue;
								auto tuple5 = table[5][ind[5]];
								calc_ripple_ci(K, ind, x, w, q, sample_values, n, a, total_size);
								sum += sample_values[ind[2]];
								test_count++;
							}
						}
					}
				}
			}
		}
		}
		//ik == n from table 2
		for(size_t i = 0; i < a[2]; ++i)
		{
		std::vector<size_t> ind(K);
		ind[2] = round * a[2] + i;
		auto tuple2 = table[2][ind[2]];
		if(cond_func[2] == nullptr || (*cond_func[2])(tuple2)){
			auto range1 = h_table[1][2].equal_range((*key_func[2][1])(tuple2));
			for(auto iter1 = range1.first; iter1 != range1.second; iter1++){
				ind[1] = iter1->second;
				if(ind[1] /a[1]>= round)
					continue;
				auto tuple1 = table[1][ind[1]];
				auto range0 = h_table[0][1].equal_range((*key_func[1][0])(tuple1));
				for(auto iter0 = range0.first; iter0 != range0.second; ++iter0){
					ind[0] = iter0->second;
					if(ind[0]/a[0]>= round)
						continue;
					auto tuple0 = table[0][ind[0]];
					auto range3 = h_table[3][2].equal_range((*key_func[2][3])(tuple2));
					for(auto iter3 = range3.first; iter3 != range3.second; ++iter3){
						ind[3] = iter3->second;
						if(ind[3]/a[3] >= round)
							continue;
						auto tuple3 = table[3][ind[3]];
						auto range4 = h_table[4][3].equal_range((*key_func[3][4])(tuple3));
						for(auto iter4 = range4.first; iter4 != range4.second; ++iter4){
							ind[4] = iter4->second;
							if(ind[4]/a[4]>= round)
								continue;
							auto tuple4 = table[4][ind[4]];
							auto range5 = h_table[5][4].equal_range((*key_func[4][5])(tuple4));
							for(auto iter5 = range5.first; iter5 != range5.second; ++iter5){
								ind[5] = iter5->second;
								if(ind[5]/a[5]>= round)
									continue;
								auto tuple5 = table[5][ind[5]];
								calc_ripple_ci(K, ind, x, w, q, sample_values, n, a, total_size);
								sum += sample_values[ind[2]];
								test_count++;
							}
						}
					}
				}
			}
		}
		}
		//! i_n from table 3
		for(size_t i = 0; i < a[3]; ++i)
		{
		std::vector<size_t> ind(K);
		ind[3] = round*a[3] + i;
		auto tuple3 = table[3][ind[3]];
		if(cond_func[3] == nullptr || (*cond_func[3])(tuple3)){
			auto range2 = h_table[2][3].equal_range((*key_func[3][2])(tuple3));
			for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
				ind[2] = iter2->second;
				if(ind[2] / a[2]>= round)
					continue;
				auto tuple2 = table[2][ind[2]];
				auto range1 = h_table[1][2].equal_range((*key_func[2][1])(tuple2));
				for(auto iter1 = range1.first; iter1 != range1.second; iter1++){
					ind[1] = iter1->second;
					if(ind[1] / a[1]>= round)
						continue;
					auto tuple1 = table[1][ind[1]];
					auto range0 = h_table[0][1].equal_range((*key_func[1][0])(tuple1));
					for(auto iter0 = range0.first; iter0 != range0.second; ++iter0){
						ind[0] = iter0->second;
						if(ind[0] / a[0]>= round)
							continue;
						auto tuple0 = table[0][ind[0]];
						auto range4 = h_table[4][3].equal_range((*key_func[3][4])(tuple3));
						for(auto iter4 = range4.first; iter4 != range4.second; ++iter4){
							ind[4] = iter4->second;
							if(ind[4]/a[4]>= round)
								continue;
							auto tuple4 = table[4][ind[4]];
							auto range5 = h_table[5][4].equal_range((*key_func[4][5])(tuple4));
							for(auto iter5 = range5.first; iter5 != range5.second; ++iter5){
								ind[5] = iter5->second;
								if(ind[5]/a[5]>= round)
									continue;
								auto tuple5 = table[5][ind[5]];
								calc_ripple_ci(K, ind, x, w, q, sample_values, n, a, total_size);
								sum += sample_values[ind[2]];
								test_count++;
							}
						}
					}
				}
			}
		}
		}
		//ik == n from table 4
		for(size_t i = 0; i < a[4]; ++i)
		{
		std::vector<size_t> ind(K);
		ind[4] = round*a[4] + i;
		auto tuple4 = table[4][ind[4]];
		if(cond_func[4] == nullptr || cond_func[4] != nullptr && (*cond_func[4])(tuple4)){
			auto range3 = h_table[3][4].equal_range((*key_func[4][3])(tuple4));
			for(auto iter3 = range3.first; iter3 != range3.second; iter3++){
				ind[3] = iter3->second;
				if(ind[3] / a[3]>= round)
					continue;
				auto tuple3 = table[3][ind[3]];
				auto range2 = h_table[2][3].equal_range((*key_func[3][2])(tuple3));
				for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
					ind[2] = iter2->second;
					if(ind[2] / a[2]>= round)
						continue;
					auto tuple2 = table[2][ind[2]];
					auto range1 = h_table[1][2].equal_range((*key_func[2][1])(tuple2));
					for(auto iter1 = range1.first; iter1 != range1.second; ++iter1){
						ind[1] = iter1->second;
						if(ind[1] / a[1]>= round)
							continue;
						auto tuple1 = table[1][ind[1]];
						auto range0 = h_table[0][1].equal_range((*key_func[1][0])(tuple1));
						for(auto iter0 = range0.first; iter0 != range0.second; ++iter0){
							ind[0] = iter0->second;
							if(ind[0] / a[0]>= round)
								continue;
							auto tuple0 = table[0][ind[0]];
							auto range5 = h_table[5][4].equal_range((*key_func[4][5])(tuple4));
							for(auto iter5 = range5.first; iter5 != range5.second; ++iter5){
								ind[5] = iter5->second;
								if(ind[5] / a[5]>= round)
									continue;
								auto tuple5 = table[5][ind[5]];
								calc_ripple_ci(K, ind, x, w, q, sample_values, n, a, total_size);
								sum += sample_values[ind[2]];
								test_count++;
							}
						}
					}
				}
			}
		}
		}
		//ik == n from table 5
		for(size_t i = 0; i < a[5]; ++i)
		{
		std::vector<size_t> ind(K);
		ind[5] = round * a[5] + i;
		auto tuple5 = table[5][ind[5]];
		if(cond_func[5] == nullptr || cond_func[5] != nullptr && (*cond_func[5])(tuple5)){
			auto range4 = h_table[4][5].equal_range((*key_func[5][4])(tuple5));
			for(auto iter4 = range4.first; iter4 != range4.second; iter4++){
				ind[4] = iter4->second;
				if(ind[4] / a[4]>= round)
					continue;
				auto tuple4 = table[4][ind[4]];
				auto range3 = h_table[3][4].equal_range((*key_func[4][3])(tuple4));
				for(auto iter3 = range3.first; iter3 != range3.second; iter3++){
					ind[3] = iter3->second;
					if(ind[3] / a[3]>= round)
						continue;
					auto tuple3 = table[3][ind[3]];
					auto range2 = h_table[2][3].equal_range((*key_func[3][2])(tuple3));
					for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
						ind[2] = iter2->second;
						if(ind[2] / a[2]>= round)
							continue;
						auto tuple2 = table[2][ind[2]];
						auto range1 = h_table[1][2].equal_range((*key_func[2][1])(tuple2));
						for(auto iter1 = range1.first; iter1 != range1.second; ++iter1){
							ind[1] = iter1->second;
							if(ind[1] / a[1]>= round)
								continue;
							auto tuple1 = table[1][ind[1]];
							auto range0 = h_table[0][1].equal_range((*key_func[1][0])(tuple1));
							for(auto iter0 = range0.first; iter0 != range0.second; ++iter0){
								ind[0] = iter0->second;
								if(ind[0] / a[0]>= round)
									continue;
								auto tuple0 = table[0][ind[0]];
								calc_ripple_ci(K, ind, x, w, q, sample_values, n, a, total_size);
									sum += sample_values[ind[2]];
									test_count++;
							}
						}
					}
				}
			}
		}
		}
		y = x[0];
		z = 0;
		double nk=std::pow(n, K);
		for(size_t i = 0; i < K; ++i)
			z += w[i] / (n * a[i]-1);
		if(steptimer.get_elapsed()> STEP){
			ci = zp*std::sqrt(z) / sqrt(n);
			tool::report(STEP*report_round, n, test_count, 0, 0, sum*total_size/nk, ci, prob);
			steptimer.restart();
			report_round++;
		}
		round++;
		}
		
}
void index_ripple_query3(
		std::vector<raw_customer> &vec_c,\
		std::vector<raw_orders> &vec_o,\
		std::vector<raw_lineitem> &vec_l,\
		std::function<double(base_raw *, base_raw *, base_raw*)> result_func,
		std::string *p1,\
		date_t *p2_start,\
		date_t *p2_end,\
		date_t *p3_start,\
		date_t *p3_end,\
		std::vector<std::vector<key_func_type>> &key_func,\
		const double STEP,\
		const size_t &MAX,\
		const double &prob,\
		const double &VALUE,\
		bool verbose){
	/* build data structure */
	std::vector<size_t> a(3, 1);
	const size_t K = 3;
	double total_size = vec_c.size() * 1.0 * vec_o.size() * vec_l.size();
	auto h_customer_mktsegment = new hash_tree<raw_customer, std::string>(vec_c, key_func_customer_mktsegment);
	auto h_orders_orderdate = new range_tree<raw_orders, date_t>(vec_o, key_func_orders_orderdate);
	auto h_lineitem_shipdate = new range_tree<raw_lineitem, date_t>(vec_l, key_func_lineitem_shipdate);
	std::vector<std::vector<base_raw *> > table;
	std::mt19937 gen;
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	gen.seed(seed);

	for(size_t i = 0; i < K; ++i){
		std::vector<base_raw *> tmp;
		table.emplace_back(tmp);
	}
	if(p1 == nullptr){
		std::uniform_int_distribution<uint64_t> dis(0, vec_c.size() - 1);
		size_t ind = dis(gen);
		table[0].push_back(&vec_c[ind]);
	}
	else{
		size_t d = 0;
		auto iter = h_customer_mktsegment->get_sample(*p1, d);
		total_size = total_size / vec_c.size() * d;
		table[0].push_back(iter);
	}
	if(p2_start == nullptr){
		std::uniform_int_distribution<uint64_t> dis(0, vec_o.size() - 1);
		size_t ind = dis(gen);
		table[1].push_back(&vec_o[ind]);
	}
	else{
		size_t d = 0;
		auto iter = h_orders_orderdate->get_sample(*p2_start, *p2_end, d);
		total_size = total_size / vec_o.size() * d;
		table[1].push_back(iter);
	}
	if(p3_start == nullptr){
		std::uniform_int_distribution<uint64_t> dis(0, vec_l.size() - 1);
		size_t ind = dis(gen);
		table[2].push_back(&vec_l[ind]);
	}
	else{
		size_t d = 0;
		auto iter = h_lineitem_shipdate->get_sample(*p3_start, *p3_end, d);
		total_size = total_size / vec_l.size() * d;
		table[2].push_back(iter);
	}

	size_t max_round = vec_l.size();;
	std::vector<double> sample_values;

	size_t test_count = 0;
	typedef unordered_multimap<uint32_t, size_t> index_type;
	std::vector<std::vector<index_type> > h_table;
	for(size_t i = 0; i < K; ++i){
		std::vector<index_type> tmp;
		h_table.push_back(tmp);
		for(size_t j = 0; j < K; ++j){
			index_type x;
			h_table[i].push_back(x);
		}
	}
	std::vector<double> q[K]; 
	std::vector<double> w(K, 0);
	std::vector<double> x(K, 0);
	for (size_t k = 0; k < K; ++k) q[k].resize(max_round);
	for (size_t k = 0; k < K; ++k)
		for(size_t i = 0; i < q[k].size(); ++i)
			q[k][i] = 0;

	cout<<total_size<<endl;
	const double zp = tool::xql_erf_inv(prob);

	double y = 0.0;
	double ci = 1.0 / 0.0;
	size_t round = 0;
	/* round 1 */
	for(size_t i = 0; i < K; ++i){
		for(size_t j = 0; j < K; ++j){
			if(key_func[i][j] != nullptr){
					h_table[i][j].insert(std::make_pair((*key_func[i][j])(table[i][round]), round)); 
			}
		}
	}
	double sum = 0;
	double v = 0;
	sample_values.push_back(result_func(table[0][round], table[1][round], table[2][round]));	
	std::vector<size_t> ind(K);
	for(size_t i = 0; i < K; ++i)
		ind[i] = round;
	if((*key_func[0][1])(table[0][ind[0]]) == (*key_func[1][0])(table[1][ind[1]]) &&
			(*key_func[1][2])(table[1][ind[1]]) == (*key_func[2][1])(table[2][ind[2]])){
				v = sample_values[round];
				y = v *total_size;
				q[0][ind[0]] += y;
				q[1][ind[1]] += y;
				q[2][ind[2]] += y;
				sum+=v;
				test_count++;
	}
	for(size_t i = 0; i < K; ++i){
		y = sum/a[i] * total_size;
	}
	for(size_t i = 0; i < K; ++i){
		x[i] = q[i][round];
	}
	for(size_t i = 0; i < K; ++i){
		w[i] = 0;
	}
	size_t n = 1;
	double z = 0;
	uint32_t n_rejected = 0;
	size_t report_round = 1;
	timer maxtimer;
	timer steptimer;
	round = 1;
	size_t old_time = 0;
	size_t cur_time;
	while(maxtimer.get_elapsed() < MAX && round < max_round){
		if(p1 == nullptr){
			std::uniform_int_distribution<uint64_t> dis(0, vec_c.size() - 1);
			size_t ind = dis(gen);
			table[0].push_back(&vec_c[ind]);
		}
		else{
			size_t d = 0;
			auto iter = h_customer_mktsegment->get_sample(*p1, d);
			table[0].push_back(iter);
		}
		if(p2_start == nullptr){
			std::uniform_int_distribution<uint64_t> dis(0, vec_o.size() - 1);
			size_t ind = dis(gen);
			table[1].push_back(&vec_o[ind]);
		}
		else{
			size_t d = 0;
			auto iter = h_orders_orderdate->get_sample(*p2_start, *p2_end, d);
			table[1].push_back(iter);
		}
		if(p3_start == nullptr){
			std::uniform_int_distribution<uint64_t> dis(0, vec_l.size() - 1);
			size_t ind = dis(gen);
			table[2].push_back(&vec_l[ind]);
		}
		else{
			size_t d = 0;
			auto iter = h_lineitem_shipdate->get_sample(*p3_start, *p3_end, d);
			table[2].push_back(iter);
		}
		n = round + 1;
		for(size_t i = 0; i < K; ++i){
			for(size_t j = 0; j < K; ++j){
				if(key_func[i][j] != nullptr)
						h_table[i][j].insert(std::make_pair((*key_func[i][j])(table[i][round]), round)); 
			}
		}
		sample_values.push_back(result_func(table[0][round], table[1][round], table[2][round]));	
		y = y * std::pow(double(round)/(round+1), K-1);
		std::vector<size_t> ind(K, round);
		if((*key_func[0][1])(table[0][ind[0]]) == (*key_func[1][0])(table[1][ind[1]]) &&
				(*key_func[1][2])(table[1][ind[1]]) == (*key_func[2][1])(table[2][ind[2]])){
					v = sample_values[round];
					auto tmp = v *total_size;
					q[0][ind[0]] += tmp;
					q[1][ind[1]] += tmp;
					q[2][ind[2]] += tmp;
					sum+=v;
					test_count++;
		}
		for(size_t i = 0; i < K; ++i){
			auto tmp = 0;
			for(size_t j = 0; j < a[i]; ++j){
				tmp += q[i][round + j];
			}
			x[i] = y + tmp/n/n;
		}

		for(size_t k = 0; k < K; ++k){
			q[k][round] = 0;
			w[k] = w[k] * std::pow(double(round) / (round+1), 4) + y * y /(round+1)/round;
		}
		//ik == n from table 0
		for(size_t i = 0; i < a[0]; ++i)
		{
		std::vector<size_t> ind(K, round);
		auto tuple0 = table[0][ind[0]];
		{
			auto range1 = h_table[1][0].equal_range((*key_func[0][1])(tuple0));
			for(auto iter1 = range1.first; iter1 != range1.second; iter1++){
				ind[1] = iter1->second;
				if(ind[1] >= round)
					continue;
				auto tuple1 = table[1][ind[1]];
				auto range2 = h_table[2][1].equal_range((*key_func[1][2])(tuple1));
				for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
					ind[2] = iter2->second;
					if(ind[2] >= round)
						continue;
					auto tuple2 = table[2][ind[2]];
					calc_ripple_ci_test(K, ind, x, w, q, sample_values, n, a, total_size);
					sum += sample_values[ind[2]];
					test_count++;
				}
			}
			}
		}
		//ik == n from table 1
		for(size_t i = 0; i < a[1]; ++i)
		{
		std::vector<size_t> ind(K, round);
		auto tuple1 = table[1][ind[1]];
		{
			auto range0 = h_table[0][1].equal_range((*key_func[1][0])(tuple1));
			for(auto iter0 = range0.first; iter0 != range0.second; iter0++){
				ind[0] = iter0->second;
				if(ind[0]>= round)
					continue;
				auto tuple0 = table[0][ind[0]];
				auto range2 = h_table[2][1].equal_range((*key_func[1][2])(tuple1));
				for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
					ind[2] = iter2->second;
					if(ind[2] >= round)
						continue;
					auto tuple2 = table[2][ind[2]];
					calc_ripple_ci_test(K, ind, x, w, q, sample_values, n, a, total_size);
					sum += sample_values[ind[2]];
					test_count++;
				}
			}
		}
		}
		//ik == n from table 2
		for(size_t i = 0; i < a[2]; ++i)
		{
		std::vector<size_t> ind(K, round);
		auto tuple2 = table[2][ind[2]];
		{
			auto range1 = h_table[1][2].equal_range((*key_func[2][1])(tuple2));
			for(auto iter1 = range1.first; iter1 != range1.second; iter1++){
				ind[1] = iter1->second;
				if(ind[1] >= round)
					continue;
				auto tuple1 = table[1][ind[1]];
				auto range0 = h_table[0][1].equal_range((*key_func[1][0])(tuple1));
				for(auto iter0 = range0.first; iter0 != range0.second; ++iter0){
					ind[0] = iter0->second;
					if(ind[0] >= round)
						continue;
					auto tuple0 = table[0][ind[0]];
					calc_ripple_ci_test(K, ind, x, w, q, sample_values, n, a, total_size);
					sum += sample_values[ind[2]];
					test_count++;
				}
			}
		}
		}
		z = 0;
		double nk=n * a[0] * n *a[1] * n *a[2];
		y = sum*total_size/nk;
		for(size_t i = 0; i < K; ++i)
			z += w[i] / (n*a[i]-1);
		ci = zp*std::sqrt(z) / sqrt(n);
		if(ci > 0 && ci < VALUE)
			return;
		cur_time = maxtimer.get_elapsed();
		if(verbose && cur_time - old_time> STEP){
			tool::report(cur_time, n, test_count, 0, 0, sum*total_size/nk, ci, prob);
			old_time = cur_time;
			report_round++;
		}
		round++;
	}
}
void index_ripple_query7(\
		std::vector<raw_supplier> &vec_s,\
		std::vector<raw_lineitem> &vec_l,\
		std::vector<raw_orders> &vec_o,\
		std::vector<raw_customer> &vec_c,\
		uint32_t *nationkey_1,\
		date_t *p2_start,\
		date_t *p2_end,\
		uint32_t *nationkey_2,\
		std::vector<std::vector<key_func_type>> &key_func,\
		std::function<double(base_raw*)> result_func,\
		const double STEP,\
		const size_t &MAX,\
		const double &prob,\
		const double &VALUE,\
		bool verbose,\
		size_t col_id){
	constexpr int K = 4;
	std::vector<size_t> a(K, 1);
	size_t max_round = vec_l.size();
	double total_size = vec_s.size() * 1.0 * vec_l.size() * vec_o.size() * vec_c.size();
	auto h_supplier_nationkey = new hash_tree<raw_supplier, uint32_t>(vec_s, key_func_supplier_nationkey);
	auto h_lineitem_shipdate = new range_tree<raw_lineitem, date_t>(vec_l, key_func_lineitem_shipdate);
	auto h_customer_nationkey = new hash_tree<raw_customer, uint32_t>(vec_c, key_func_customer_nationkey);
	std::vector<double> sample_values;

	std::mt19937 gen;
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	gen.seed(seed);
	std::vector<std::vector<base_raw *> > table;
	for(size_t i = 0; i < K; ++i){
		std::vector<base_raw *> tmp;
		table.emplace_back(tmp);
	}
	if(nationkey_1 == nullptr){
		std::uniform_int_distribution<uint64_t> dis(0, vec_s.size() - 1);
		size_t ind = dis(gen);
		table[0].push_back(&vec_s[ind]);
	}
	else{
		size_t d = 0;
		auto iter = h_supplier_nationkey->get_sample(*nationkey_1, d);
		total_size = total_size / vec_s.size() * d;
		table[0].push_back(iter);
	}
	if(p2_start == nullptr){
		std::uniform_int_distribution<uint64_t> dis(0, vec_l.size() - 1);
		size_t ind = dis(gen);
		table[1].push_back(&vec_l[ind]);
	}
	else{
		size_t d = 0;
		auto iter = h_lineitem_shipdate->get_sample(*p2_start, *p2_end, d);
		total_size = total_size / vec_l.size() * d;
		table[1].push_back(iter);
	}
	{
		std::uniform_int_distribution<uint64_t> dis(0, vec_o.size() - 1);
		size_t ind = dis(gen);
		table[2].push_back(&vec_o[ind]);
	}
	if(nationkey_2 == nullptr){
		std::uniform_int_distribution<uint64_t> dis(0, vec_c.size() - 1);
		size_t ind = dis(gen);
		table[3].push_back(&vec_c[ind]);
	}
	else{
		size_t d = 0;
		auto iter = h_customer_nationkey->get_sample(*nationkey_2, d);
		total_size = total_size /vec_c.size() * d;
		table[3].push_back(iter);
	}
	tool::start_report();

	size_t test_count = 0;
	typedef unordered_multimap<uint32_t, size_t> index_type;
	std::vector<std::vector<index_type> > h_table;
	for(size_t i = 0; i < K; ++i){
		std::vector<index_type> tmp;
		h_table.push_back(tmp);
		for(size_t j = 0; j < K; ++j){
			index_type x;
			h_table[i].push_back(x);
		}
	}
	std::vector<double> q[K]; 
	std::vector<double> w(K, 0);
	std::vector<double> x(K, 0);
	for (size_t k = 0; k < K; ++k) q[k].resize(max_round * a[k]);
	for (size_t k = 0; k < K; ++k)
		for(size_t i = 0; i < q[k].size(); ++i)
			q[k][i] = 0;

	cout<<total_size<<endl;

	const double zp = tool::xql_erf_inv(prob);

	double y = 0.0;
	double ci = 1.0 / 0.0;
	size_t round = 0;
	/* round 0 */
	for(size_t i = 0; i < K; ++i){
		for(size_t k = 0; k < a[i]; ++k){
			for(size_t j = 0; j < K; ++j){
					if(key_func[i][j] != nullptr)
							h_table[i][j].insert(std::make_pair((*key_func[i][j])(table[i][round * a[i] + k]), round*a[i]+k)); 
			}
		}
	}
	double sum = 0;
	double v = 0;
	for(size_t i = 0; i < a[col_id]; ++i)
			sample_values.push_back(result_func(table[col_id][round*a[col_id] + i]));	
	size_t i[4];
	for(i[0] = 0; i[0] < a[0]; ++i[0]){
		for(i[1] = 0; i[1] < a[1]; ++i[1]){
			for(i[2] = 0; i[2] < a[2]; ++i[2]){
				for(i[3] = 0; i[3] < a[3]; ++i[3]){
					std::vector<size_t> ind(K);
					ind[0] = round*a[0] + i[0];
					ind[1] = round*a[1] + i[1];
					ind[2] = round*a[2] + i[2];
					ind[3] = round*a[3] + i[3];
				if((*key_func[0][1])(table[0][ind[0]]) == (*key_func[1][0])(table[1][ind[1]]) \
						&& (*key_func[1][2])(table[1][ind[1]]) == (*key_func[2][1])(table[2][ind[2]]) \
						&& (*key_func[2][3])(table[2][ind[2]])== (*key_func[3][2])(table[3][ind[3]])){
							v = sample_values[round*a[col_id] + i[col_id]];
							y = v * total_size;
							q[0][ind[0]] += y;
							q[1][ind[1]] += y;
							q[2][ind[2]] += y;
							q[3][ind[3]] += y;
							sum+=v;
							test_count++;
				}
				}
			}
		}
	}
	for(size_t i = 0; i < K; ++i){
		y = sum/a[i] * total_size;
	}
	double pow_n[K];
	double all = 1;
	for(size_t i = 0; i < K; ++i)
		all *= (round+1) * a[i];
	for(size_t i = 0; i < K; ++i)
		pow_n[i] = all/(round+1)/a[i];
	for(size_t i = 0; i < K; ++i){
		for(size_t j = 0; j < a[i]; ++j){
			x[i] += q[i][round*a[i]+ j]/a[i];
		}
	}
	for(size_t i = 0; i < K; ++i){
		w[i] = 0;
	}
	size_t n = 1;
	double z = 0;
	uint32_t n_rejected = 0;
	size_t report_round = 1;
	timer maxtimer;
	timer steptimer;
	size_t old_time = 0;
	round = 1;
	while(maxtimer.get_elapsed() < MAX && round < max_round){
		n = round + 1;
		if(nationkey_1 == nullptr){
			std::uniform_int_distribution<uint64_t> dis(0, vec_s.size() - 1);
			size_t ind = dis(gen);
			table[0].push_back(&vec_s[ind]);
		}
		else{
			size_t d = 0;
			auto iter = h_supplier_nationkey->get_sample(*nationkey_1, d);
			table[0].push_back(iter);
		}
		if(p2_start == nullptr){
			std::uniform_int_distribution<uint64_t> dis(0, vec_l.size() - 1);
			size_t ind = dis(gen);
			table[1].push_back(&vec_l[ind]);
		}
		else{
			size_t d = 0;
			auto iter = h_lineitem_shipdate->get_sample(*p2_start, *p2_end, d);
			table[1].push_back(iter);
		}
		{
			std::uniform_int_distribution<uint64_t> dis(0, vec_o.size() - 1);
			size_t ind = dis(gen);
			table[2].push_back(&vec_o[ind]);
		}
		if(nationkey_2 == nullptr){
			std::uniform_int_distribution<uint64_t> dis(0, vec_c.size() - 1);
			size_t ind = dis(gen);
			table[3].push_back(&vec_c[ind]);
		}
		else{
			size_t d = 0;
			auto iter = h_customer_nationkey->get_sample(*nationkey_2, d);
			table[3].push_back(iter);
		}
		for(size_t i = 0; i < K; ++i)
			all *= (round+1) * a[i];
		for(size_t i = 0; i < K; ++i)
			pow_n[i] = all/(round+1)/a[i];
		for(size_t i = 0; i < K; ++i){
			for(size_t k = 0; k < a[i]; ++k){
				for(size_t j = 0; j < K; ++j){
					if(key_func[i][j] != nullptr){
							h_table[i][j].insert(std::make_pair((*key_func[i][j])(table[i][round * a[i] + k]), round*a[i]+k)); 
					}
				}
			}
		}
		for(size_t i = 0; i < a[col_id]; ++i)
			sample_values.push_back(result_func(table[col_id][round*a[col_id] + i]));	
		y = y * std::pow(double(round)/(round+1), K-1);
		for(size_t k = 0; k < K; ++k){
			for(size_t i = 0; i < a[k]; ++i)
				q[k][round*a[k] + i ] = 0;
		}
		for(i[0] = 0; i[0] < a[0]; ++i[0]){
			for(i[1] = 0; i[1] < a[1]; ++i[1]){
				for(i[2] = 0; i[2] < a[2]; ++i[2]){
					for(i[3] = 0; i[3] < a[3]; ++i[3]){
						std::vector<size_t> ind(K);
						ind[0] = round*a[0] + i[0];
						ind[1] = round*a[1] + i[1];
						ind[2] = round*a[2] + i[2];
						ind[3] = round*a[3] + i[3];
					if((*key_func[0][1])(table[0][ind[0]]) == (*key_func[1][0])(table[1][ind[1]]) \
							&& (*key_func[1][2])(table[1][ind[1]]) == (*key_func[2][1])(table[2][ind[2]]) \
							&& (*key_func[2][3])(table[2][ind[2]])== (*key_func[3][2])(table[3][ind[3]])){
								v = sample_values[round*a[col_id] + i[col_id]];
								auto tmp = v *total_size;
								q[0][ind[0]] += tmp;
								q[1][ind[1]] += tmp;
								q[2][ind[2]] += tmp;
								q[3][ind[3]] += tmp;
								sum += v;
								test_count++;
					}
					}
				}
			}
		}

		for(size_t i = 0; i < K; ++i){
			auto tmp = 0;
			for(size_t j = 0; j < a[i]; ++j){
				tmp += q[i][round*a[i] + j];
			}
			x[i] = y + tmp/a[i]/ pow_n[i];
		}

		for(size_t k = 0; k < K; ++k){
			for(size_t i = 0; i < a[k]; ++i)
				q[k][round*a[k] + i ] = 0;
			w[k] = w[k] * std::pow(double(round) / (round+1), 6) + y * y /(round+1)/round;
		}
		//ik == n from table 0
		for(size_t i = 0; i < a[0]; ++i)
		{
		std::vector<size_t> ind(K);
		ind[0] = round*a[0] + i;
		auto tuple0 = table[0][ind[0]];
		{
			auto range1 = h_table[1][0].equal_range((*key_func[0][1])(tuple0));
			for(auto iter1 = range1.first; iter1 != range1.second; iter1++){
				ind[1] = iter1->second;
				if(ind[1]/a[1] >= round)
					continue;
				auto tuple1 = table[1][ind[1]];
				auto range2 = h_table[2][1].equal_range((*key_func[1][2])(tuple1));
				for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
					ind[2] = iter2->second;
					if(ind[2]/a[2] >= round)
						continue;
					auto tuple2 = table[2][ind[2]];
					auto range3 = h_table[3][2].equal_range((*key_func[2][3])(tuple2));
					for(auto iter3 = range3.first; iter3 != range3.second; ++iter3){
						ind[3] = iter3->second;
						if(ind[3]/a[3] >= round)
							continue;
						auto tuple3 = table[3][ind[3]];
						calc_ripple_ci_test(K, ind, x, w, q, sample_values, n, a, total_size, col_id);
						sum += sample_values[ind[col_id]];
						test_count++;
					}
				}
			}
			}
		}
		//ik == n from table 1
		for(size_t i = 0; i < a[1]; ++i)
		{
		std::vector<size_t> ind(K);
		ind[1] = round * a[1] + i;
		auto tuple1 = table[1][ind[1]];
		{
			auto range0 = h_table[0][1].equal_range((*key_func[1][0])(tuple1));
			for(auto iter0 = range0.first; iter0 != range0.second; iter0++){
				ind[0] = iter0->second;
				if(ind[0]/a[0]>= round)
					continue;
				auto tuple0 = table[0][ind[0]];
				auto range2 = h_table[2][1].equal_range((*key_func[1][2])(tuple1));
				for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
					ind[2] = iter2->second;
					if(ind[2]/a[2]>= round)
						continue;
					auto tuple2 = table[2][ind[2]];
					auto range3 = h_table[3][2].equal_range((*key_func[2][3])(tuple2));
					for(auto iter3 = range3.first; iter3 != range3.second; ++iter3){
						ind[3] = iter3->second;
						if(ind[3]/a[3] >= round)
							continue;
						auto tuple3 = table[3][ind[3]];
						calc_ripple_ci_test(K, ind, x, w, q, sample_values, n, a, total_size, col_id);
						sum += sample_values[ind[col_id]];
						test_count++;
					}
				}
			}
		}
		}
		//ik == n from table 2
		for(size_t i = 0; i < a[2]; ++i)
		{
		std::vector<size_t> ind(K);
		ind[2] = a[2]*round + i;
		auto tuple2 = table[2][ind[2]];
		{
			auto range3 = h_table[3][2].equal_range((*key_func[2][3])(tuple2));
			for(auto iter3 = range3.first; iter3 != range3.second; ++iter3){
				ind[3] = iter3->second;
				if(ind[3]/a[3]>= round)
					continue;
				auto tuple3 = table[3][ind[3]];
				auto range1 = h_table[1][2].equal_range((*key_func[2][1])(tuple2));
				for(auto iter1 = range1.first; iter1 != range1.second; iter1++){
					ind[1] = iter1->second;
					if(ind[1]/a[1] >= round)
						continue;
					auto tuple1 = table[1][ind[1]];
					auto range0 = h_table[0][1].equal_range((*key_func[1][0])(tuple1));
					for(auto iter0 = range0.first; iter0 != range0.second; ++iter0){
						ind[0] = iter0->second;
						if(ind[0]/a[0] >= round)
							continue;
						auto tuple0 = table[0][ind[0]];
						calc_ripple_ci_test(K, ind, x, w, q, sample_values, n, a, total_size, col_id);
						sum += sample_values[ind[col_id]];
						test_count++;
					}
				}
			}
		}
		}
		//! i_n from table 3
		for(size_t i = 0; i < a[3]; ++i)
		{
		std::vector<size_t> ind(K);
		ind[3] = round*a[3] + i;
		auto tuple3 = table[3][ind[3]];
		{
			auto range2 = h_table[2][3].equal_range((*key_func[3][2])(tuple3));
			for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
				ind[2] = iter2->second;
				if(ind[2]/a[2]>= round)
					continue;
				auto tuple2 = table[2][ind[2]];
				auto range1 = h_table[1][2].equal_range((*key_func[2][1])(tuple2));
				for(auto iter1 = range1.first; iter1 != range1.second; iter1++){
					ind[1] = iter1->second;
					if(ind[1]/a[1]>= round)
						continue;
					auto tuple1 = table[1][ind[1]];
					auto range0 = h_table[0][1].equal_range((*key_func[1][0])(tuple1));
					for(auto iter0 = range0.first; iter0 != range0.second; ++iter0){
						ind[0] = iter0->second;
						if(ind[0]/a[0] >= round)
							continue;
						auto tuple0 = table[0][ind[0]];
						calc_ripple_ci_test(K, ind, x, w, q, sample_values, n, a, total_size, col_id);
						sum += sample_values[ind[col_id]];
						test_count++;
					}
				}
			}
		}
		}
		z = 0;
		double nk = 1;
		for(size_t i = 0; i < K; ++i)
			nk *= a[i] * n;
		y = sum/nk*total_size;;
		for(size_t i = 0; i < K; ++i)
			z += w[i] / (n * a[i]-1);
		auto cur_time = maxtimer.get_elapsed();
		ci = zp*std::sqrt(z) / sqrt(n);
		if(ci > 0 && ci < VALUE)
			return;
		if(verbose && cur_time - old_time> STEP){
			tool::report(cur_time, n, test_count, 0, 0, sum/nk*total_size, ci, prob);
			old_time = cur_time;
			report_round++;
		}
		round++;
	}
}
void index_ripple_query10(\
		std::vector<raw_customer> &vec_c,\
		std::vector<raw_orders> &vec_o,\
		std::vector<raw_lineitem> &vec_l,\
		date_t *p2_start,\
		date_t *p2_end,\
		unsigned char *flag,\
		std::vector<std::vector<key_func_type>> &key_func,\
		std::function<double(base_raw*, base_raw*, base_raw*)> result_func,\
		const double STEP,\
		const size_t &MAX,\
		const double &prob,\
		const double &VALUE,\
		bool verbose){
	/* build data structure */
	std::vector<size_t> a(3, 1);
	const size_t K = 3;
	double total_size = vec_c.size() * 1.0 * vec_o.size() * vec_l.size();
	auto h_orders_orderdate = new range_tree<raw_orders, date_t>(vec_o, key_func_orders_orderdate);
	auto h_lineitem_returnflag = new hash_tree<raw_lineitem, unsigned char>(vec_l, key_func_lineitem_returnflag);
	std::vector<std::vector<base_raw *> > table;
	std::mt19937 gen;
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	gen.seed(seed);

	for(size_t i = 0; i < K; ++i){
		std::vector<base_raw *> tmp;
		table.emplace_back(tmp);
	}
	{
		std::uniform_int_distribution<uint64_t> dis(0, vec_c.size() - 1);
		size_t ind = dis(gen);
		table[0].push_back(&vec_c[ind]);
	}
	if(p2_start == nullptr){
		std::uniform_int_distribution<uint64_t> dis(0, vec_o.size() - 1);
		size_t ind = dis(gen);
		table[1].push_back(&vec_o[ind]);
	}
	else{
		size_t d = 0;
		auto iter = h_orders_orderdate->get_sample(*p2_start, *p2_end, d);
		total_size = total_size / vec_o.size() * d;
		table[1].push_back(iter);
	}
	if(flag == nullptr){
		std::uniform_int_distribution<uint64_t> dis(0, vec_l.size() - 1);
		size_t ind = dis(gen);
		table[2].push_back(&vec_l[ind]);
	}
	else{
		size_t d = 0;
		auto iter = h_lineitem_returnflag->get_sample(*flag, d);
		total_size = total_size / vec_l.size() * d;
		table[2].push_back(iter);
	}

	size_t max_round = vec_l.size();
	std::vector<double> sample_values;

	size_t test_count = 0;
	typedef unordered_multimap<uint32_t, size_t> index_type;
	std::vector<std::vector<index_type> > h_table;
	for(size_t i = 0; i < K; ++i){
		std::vector<index_type> tmp;
		h_table.push_back(tmp);
		for(size_t j = 0; j < K; ++j){
			index_type x;
			h_table[i].push_back(x);
		}
	}
	std::vector<double> q[K]; 
	std::vector<double> w(K, 0);
	std::vector<double> x(K, 0);
	for (size_t k = 0; k < K; ++k) q[k].resize(max_round);
	for (size_t k = 0; k < K; ++k)
		for(size_t i = 0; i < q[k].size(); ++i)
			q[k][i] = 0;

	cout<<total_size<<endl;
	const double zp = tool::xql_erf_inv(prob);

	double y = 0.0;
	double ci = 1.0 / 0.0;
	size_t round = 0;
	/* round 1 */
	for(size_t i = 0; i < K; ++i){
		for(size_t j = 0; j < K; ++j){
			if(key_func[i][j] != nullptr){
					h_table[i][j].insert(std::make_pair((*key_func[i][j])(table[i][round]), round)); 
			}
		}
	}
	double sum = 0;
	double v = 0;
	sample_values.push_back(result_func(table[0][round], table[1][round], table[2][round]));	
	std::vector<size_t> ind(K);
	for(size_t i = 0; i < K; ++i)
		ind[i] = round;
	if((*key_func[0][1])(table[0][ind[0]]) == (*key_func[1][0])(table[1][ind[1]]) &&
			(*key_func[1][2])(table[1][ind[1]]) == (*key_func[2][1])(table[2][ind[2]])){
				v = sample_values[round];
				y = v *total_size;
				q[0][ind[0]] += y;
				q[1][ind[1]] += y;
				q[2][ind[2]] += y;
				sum+=v;
				test_count++;
	}
	for(size_t i = 0; i < K; ++i){
		y = sum/a[i] * total_size;
	}
	for(size_t i = 0; i < K; ++i){
		x[i] = q[i][round];
	}
	for(size_t i = 0; i < K; ++i){
		w[i] = 0;
	}
	size_t n = 1;
	double z = 0;
	uint32_t n_rejected = 0;
	size_t report_round = 1;
	timer maxtimer;
	timer steptimer;
	round = 1;
	size_t old_time = 0;
	size_t cur_time;
	while(maxtimer.get_elapsed() < MAX && round < max_round){
		{
			std::uniform_int_distribution<uint64_t> dis(0, vec_c.size() - 1);
			size_t ind = dis(gen);
			table[0].push_back(&vec_c[ind]);
		}
		if(p2_start == nullptr){
			std::uniform_int_distribution<uint64_t> dis(0, vec_o.size() - 1);
			size_t ind = dis(gen);
			table[1].push_back(&vec_o[ind]);
		}
		else{
			size_t d = 0;
			auto iter = h_orders_orderdate->get_sample(*p2_start, *p2_end, d);
			table[1].push_back(iter);
		}
		if(flag == nullptr){
			std::uniform_int_distribution<uint64_t> dis(0, vec_l.size() - 1);
			size_t ind = dis(gen);
			table[2].push_back(&vec_l[ind]);
		}
		else{
			size_t d = 0;
			auto iter = h_lineitem_returnflag->get_sample(*flag, d);
			table[2].push_back(iter);
		}
		n = round + 1;
		for(size_t i = 0; i < K; ++i){
			for(size_t j = 0; j < K; ++j){
				if(key_func[i][j] != nullptr)
						h_table[i][j].insert(std::make_pair((*key_func[i][j])(table[i][round]), round)); 
			}
		}
		sample_values.push_back(result_func(table[0][round], table[1][round], table[2][round]));	
		y = y * std::pow(double(round)/(round+1), K-1);
		std::vector<size_t> ind(K, round);
		if((*key_func[0][1])(table[0][ind[0]]) == (*key_func[1][0])(table[1][ind[1]]) &&
				(*key_func[1][2])(table[1][ind[1]]) == (*key_func[2][1])(table[2][ind[2]])){
					v = sample_values[round];
					auto tmp = v *total_size;
					q[0][ind[0]] += tmp;
					q[1][ind[1]] += tmp;
					q[2][ind[2]] += tmp;
					sum+=v;
					test_count++;
		}
		for(size_t i = 0; i < K; ++i){
			auto tmp = 0;
			for(size_t j = 0; j < a[i]; ++j){
				tmp += q[i][round + j];
			}
			x[i] = y + tmp/n/n;
		}

		for(size_t k = 0; k < K; ++k){
			q[k][round] = 0;
			w[k] = w[k] * std::pow(double(round) / (round+1), 4) + y * y /(round+1)/round;
		}
		//ik == n from table 0
		for(size_t i = 0; i < a[0]; ++i)
		{
		std::vector<size_t> ind(K, round);
		auto tuple0 = table[0][ind[0]];
		{
			auto range1 = h_table[1][0].equal_range((*key_func[0][1])(tuple0));
			for(auto iter1 = range1.first; iter1 != range1.second; iter1++){
				ind[1] = iter1->second;
				if(ind[1] >= round)
					continue;
				auto tuple1 = table[1][ind[1]];
				auto range2 = h_table[2][1].equal_range((*key_func[1][2])(tuple1));
				for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
					ind[2] = iter2->second;
					if(ind[2] >= round)
						continue;
					auto tuple2 = table[2][ind[2]];
					calc_ripple_ci_test(K, ind, x, w, q, sample_values, n, a, total_size);
					sum += sample_values[ind[2]];
					test_count++;
				}
			}
			}
		}
		//ik == n from table 1
		for(size_t i = 0; i < a[1]; ++i)
		{
		std::vector<size_t> ind(K, round);
		auto tuple1 = table[1][ind[1]];
		{
			auto range0 = h_table[0][1].equal_range((*key_func[1][0])(tuple1));
			for(auto iter0 = range0.first; iter0 != range0.second; iter0++){
				ind[0] = iter0->second;
				if(ind[0]>= round)
					continue;
				auto tuple0 = table[0][ind[0]];
				auto range2 = h_table[2][1].equal_range((*key_func[1][2])(tuple1));
				for(auto iter2 = range2.first; iter2 != range2.second; ++iter2){
					ind[2] = iter2->second;
					if(ind[2] >= round)
						continue;
					auto tuple2 = table[2][ind[2]];
					calc_ripple_ci_test(K, ind, x, w, q, sample_values, n, a, total_size);
					sum += sample_values[ind[2]];
					test_count++;
				}
			}
		}
		}
		//ik == n from table 2
		for(size_t i = 0; i < a[2]; ++i)
		{
		std::vector<size_t> ind(K, round);
		auto tuple2 = table[2][ind[2]];
		{
			auto range1 = h_table[1][2].equal_range((*key_func[2][1])(tuple2));
			for(auto iter1 = range1.first; iter1 != range1.second; iter1++){
				ind[1] = iter1->second;
				if(ind[1] >= round)
					continue;
				auto tuple1 = table[1][ind[1]];
				auto range0 = h_table[0][1].equal_range((*key_func[1][0])(tuple1));
				for(auto iter0 = range0.first; iter0 != range0.second; ++iter0){
					ind[0] = iter0->second;
					if(ind[0] >= round)
						continue;
					auto tuple0 = table[0][ind[0]];
					calc_ripple_ci_test(K, ind, x, w, q, sample_values, n, a, total_size);
					sum += sample_values[ind[2]];
					test_count++;
				}
			}
		}
		}
		z = 0;
		double nk=n * a[0] * n *a[1] * n *a[2];
		y = sum*total_size/nk;
		for(size_t i = 0; i < K; ++i)
			z += w[i] / (n*a[i]-1);
		ci = zp*std::sqrt(z) / sqrt(n);
		if(ci > 0 && ci < VALUE)
			return;
		cur_time = maxtimer.get_elapsed();
		if(verbose && cur_time - old_time> STEP){
			tool::report(cur_time, n, test_count, 0, 0, sum*total_size/nk, ci, prob);
			old_time = cur_time;
			report_round++;
		}
		round++;
	}
}
#endif
