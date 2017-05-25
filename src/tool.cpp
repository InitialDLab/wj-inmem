#include "tool.h"
std::vector<std::vector<size_t> > tool::full_permutation(const size_t num){
	std::vector<size_t> order(num);
	std::vector<std::vector<size_t> > result;
	for(size_t i = 0; i < num; ++i)
		order.push_back(i);
	do {
		result.push_back(order);
	}while(std::prev_permutation(order.begin(), order.end()));
	return result;
}
vector<string> tool::split(string &str, const char*c)
{
	char *cstr, *p;
	vector<string> res;
	cstr = new char[str.size() + 1];
	strcpy(cstr, str.c_str());
	p = strtok(cstr, c);
	while(p != NULL)
	{
		res.push_back(p);
		p = strtok(NULL,c);
	}
	return res;
}
void tool::start_report() {
	printf("-----------------------------------------------------------------------------------------------------------\n");
	printf("%16s %16s %16s %16s %16s %16s %16s %16s\n",
			"Time",
			"#. Sample",
			"#. Rejected(join)",
			"#. Rejected(cond)",
			"Frequency",
			"Result",
			"CI",
			"Prob.");
}
void tool::report(const size_t &time, const size_t &n, const size_t &n_rejected_join, const size_t &n_rejected_cond, const double &freq, const double &result, const double &ci, const double &prob){
	printf("%16zu %16zu %16zu %16zu %16.3f %16.3f %16.3f %16.2f\n", time, n, n_rejected_join, n_rejected_cond, freq, result, ci, prob); 
}

void tool::report(const std::list<double> &est, const size_t &time, const size_t &n_rejected_join, const size_t
		&n_rejected_cond, const double &prob){
	size_t n = est.size();
	double y = tool::get_mean(est);
	double s2 =  tool::get_variance(est, y);
	double zp = tool::get_zp(prob);
	double ci =  sqrt(zp * zp *s2/n);
	double freq = double(n_rejected_join + n_rejected_cond)/ n;
	printf("%16zu %16zu %16zu %16zu %16.3f %16.3f %16.3f %16.2f\n", time, n, n_rejected_join, n_rejected_cond, freq, y, ci, prob); 
}
double tool::get_variance(const std::list<double> &est, const double &mean){
	double s = 0;
	for(auto &x:est)
		s += std::pow(x - mean, 2);
	size_t n = est.size();
	s = s / n /(n-1);
	return s;
}
double tool::get_mean(const std::list<double> &est){
	return std::accumulate(est.begin(), est.end(), 0.0) / est.size();
}
double tool::get_zp(const double &p){
	return xql_erf_inv(p) * root_2;
}
double tool::loose_ci(const double &sum_y2, const double &y_hat, const double &round, const double &prob){
		auto sigma = std::sqrt(sum_y2/round + y_hat * y_hat);
		auto ci = sigma / std::sqrt(round) * get_zp(prob); 
		return ci;
}
double tool::calc_var(const double &sum_y2, const double &sum_y, const double &y_hat, const double &round){
		return sum_y2/(round-1) + y_hat * y_hat * round / (round-1) - 2 * sum_y * y_hat/(round-1);
}
double tool::calc_ci(const double &sum_y2, const double &sum_y, const double &y_hat, const double &round, const double &prob){
	if(round <= 1)
		return 0.0;
		auto sigma = std::sqrt(calc_var(sum_y2, sum_y, y_hat, round));
		auto ci = sigma / std::sqrt(round) * get_zp(prob); 
		return ci;
}
double tool::calc_ci(const double &var, const double &round, const double &prob){
	if(round <= 1)
		return 0.0;
		auto sigma = std::sqrt(var);
		auto ci = sigma / std::sqrt(round) * get_zp(prob); 
		return ci;
}
std::vector<std::vector<size_t> > tool::gen_order(const uint32_t &n){
	std::vector<std::vector<size_t>> result;
	std::vector<size_t> tmp(n);
	for(size_t i = 0; i < n; ++i)
		tmp[i] = i;
	do{
		result.push_back(tmp);
	}while(std::next_permutation(tmp.begin(), tmp.end()));
	return result;
}
std::vector<std::vector<size_t> > tool::back_trace(const std::vector<std::vector<size_t> > &order, const std::vector<std::vector<size_t> > &adj){
	std::vector<std::vector<size_t> > result;
	auto  n = adj.size();
	auto m = order.size();
	std::vector<bool> tag(n, false);
	for(size_t i = 0; i < m; ++i){
		for(auto iter = tag.begin(); iter != tag.end(); ++iter)
			*iter = false;
		tag[order[i][0]] = true;
		for(size_t j = 1; j < n; ++j){
			size_t num_k = adj[order[i][j]].size();
			for(size_t k = 0; k < num_k; ++k){
				if(tag[adj[order[i][j]][k]] == true){
					tag[order[i][j]] = true;
					break;
				}
			}
		}
		if(std::all_of(tag.begin(), tag.end(), [](const uint32_t &input){ return input == true;}))
			result.push_back(order[i]);
	}
	return result;
}

void tool::read_date_t(date_t &x, date_t &y){
	std::ifstream fin;
	fin.open("config.txt",std::ifstream::in);
	size_t tmp;
	{
		fin>>tmp;
		x.year = tmp;
		fin>>tmp;
		x.month = tmp;
		fin>>tmp;
		x.day = tmp;
		fin>>tmp;
		y.year = tmp;
		fin>>tmp;
		y.month = tmp;
		fin>>tmp;
		y.day = tmp;
	}
	fin.close();
}
/*
date_t tool::read_config(std::string& name, std::string& tag){
	std::ifstream fin;
	size_t tmp;
	std::string t;
	date_t r;
	fin.open(name.c_str(), std::ifstream::in);
	while(!fin.eof()){
		fin>>t;
		while(tag.compare(t)){
			fin>>tmp;
			fin>>tmp;
			fin>>tmp;
			fin>>t;
		}
		fin>>tmp;
		r.year = tmp;
		fin>>tmp;
		r.month = tmp;
		fin>>tmp;
		r.day = tmp;
	}
	return r;
}
*/


