#ifndef _H_JOIN
#define _H_JOIN

#include <iostream>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <random>
#include <chrono>
#include <cstdlib>

using namespace std;

template<class K, class T>
void unordered_multimap_to_vector(const unordered_multimap<K, T> &m, vector<T> &v) {
	v.clear();
	v.reserve(m.size());
	transform(m.begin(), m.end(), back_inserter(v),
			[](const pair<K, T>& p) -> T {
				return p.second;
			});
}

class join{
public:
	static unsigned seed;
	static std::mt19937 gen;


	static void init_rand();
	static uint64_t getSample(const uint64_t m)
	{
		std::uniform_int_distribution<uint64_t> dis(0, m-1);
		return dis(gen);
	}

	void start_report(int report_round, int fflush_round) {
		m_report_round = report_round;
		if (fflush_round == 0)
			m_fflush_round = m_report_round;
		else
			m_fflush_round = fflush_round;

		if (m_report_round != 0)
		printf("%16s %16s %16s %16s %16s\n",
				"round",
				"nrejected",
				"result",
				"prob",
				"time");
	}

	bool need_report(int n) {
		return m_report_round != 0 && n % m_report_round == 0;
	}

	void report(int n, int nrejected, double y, double ci, double time_elapsed) {
		if (m_report_round != 0 && n % m_report_round == 0) {
			printf("%16d %16d %16e %16.6lf %16lf\n", n, nrejected, y, ci, time_elapsed);
			if (n % m_fflush_round == 0) 
				fflush(stdout);
		}
	}

	void end_report() {
		fflush(stdout);
	}

	int m_report_round;
	int m_fflush_round;
};
#endif
