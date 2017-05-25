#include "tool.h"
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
double tool::calc_ci(const double &sum_y2, const double &sum_y, const double &y_hat, const double &round, const double &prob){
		auto sigma = std::sqrt(sum_y2/(round-1) + y_hat * y_hat * round / (round-1) - 2 * sum_y * y_hat/(round-1));
		auto ci = sigma / std::sqrt(round) * get_zp(prob); 
		return ci;
}
