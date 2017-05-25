#include "rawdata.h"
#include "hashjoin.h"
#include "ripplejoin.h"
#include "samplejoin.h"
#include "loaddata.h"
#include "base.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include "timer.h"
using namespace std;

void do_hashjoin(){
	auto vec_c = load_customer();
	auto vec_o = load_orders();
	auto vec_l = load_lineitem();
	const date_t date_condition = {1995, 3, 15};
	cout<<"customer size: "<<vec_c.size()<<endl;
	cout<<"order size: "<<vec_o.size()<<endl;
	cout<<"lineitem size: "<<vec_l.size()<<endl;
	size_t STEP = 10000; //! 10ms step
	size_t MAX = 3000000; //! max time span
	date_t date_cond_orders = {1995, 2, 1};
	date_t date_cond_lineitem = {1995, 3, 15};
	//std::function<bool(const raw_customer &)> cond_c = [&](const raw_customer &t)->bool{ return strcmp(t.c_mktsegment,"BUILDING") == 0;};
	std::function<bool(const raw_customer &)> cond_c = [&](const raw_customer &t)->bool{ return true;};
	//std::function<bool(const raw_orders &)> cond_o = [&](const raw_orders &t) ->bool {return t.o_orderdate <date_cond_orders;};
	std::function<bool(const raw_orders &)> cond_o = [&](const raw_orders &t) ->bool {return true;};
	//std::function<bool(const raw_lineitem &)> cond_l = [&](const raw_lineitem &t) ->bool {return date_cond_lineitem < t.l_shipdate;};
	std::function<bool(const raw_lineitem &)> cond_l = [&](const raw_lineitem &t) ->bool {return true;};
	timer mytimer;
	auto tmp = hash_join<raw_customer, raw_orders, raw_customer_orders, decltype(cond_c), decltype(cond_o),
		 decltype(key_func_customer_custkey), decltype(key_func_orders_custkey), uint32_t, uint32_t>(vec_c, vec_o,
				 cond_c, cond_o, key_func_customer_custkey, key_func_orders_custkey);

	auto result = hash_join_right<raw_customer_orders, raw_lineitem, raw_customer_orders_lineitem,
		 decltype(cond_l),decltype(key_func_customer_orders_orderkey), decltype(key_func_lineitem_orderkey),uint32_t,
		 uint32_t>(tmp, vec_l, cond_l, key_func_customer_orders_orderkey, key_func_lineitem_orderkey);
	cout<<"Time cost: "<< mytimer.get_elapsed()<<endl;
	auto calc = [&](double &t1, raw_customer_orders_lineitem &t2)->double{
		return t1 + t2.dos.l_extendedprice * (1.0 - t2.dos.l_discount);
	};
	cout<<std::accumulate(result.begin(), result.end(), 0.0, calc)<<endl;
}
void do_samplejoin_no_cond(){
		auto  vec_c = load_customer();
		auto vec_o = load_orders();
		auto vec_l = load_lineitem();
		auto result_func = [&](const raw_lineitem &t)->double{return t.l_extendedprice * (1.0 - t.l_discount);};
		size_t STEP = 10000; //! 10ms step
		size_t MAX = 3000000; //! max time span
		samplejoin<raw_customer, \
			raw_orders, \
			raw_lineitem, \
			std::function<std::string(const raw_customer &)>, \
			std::function<uint32_t(const raw_orders)>, \
			std::function<uint32_t(const raw_lineitem)>,\
			std::string,\
			uint32_t,\
			uint32_t,\
			std::function<double(const raw_lineitem &)> >
			(vec_c, \
				vec_o,\
				vec_l,\
				key_func_orders_custkey,\
				key_func_lineitem_orderkey,\
				key_func_customer_custkey,\
				key_func_orders_orderkey,\
				result_func,
				STEP,
				MAX,
				.95);
}
void do_samplejoin(){
		auto result_func = [&](const raw_lineitem &t)->double{return t.l_extendedprice * (1.0 - t.l_discount);};
		auto  vec_c = load_customer();
		auto vec_o = load_orders();
		auto vec_l = load_lineitem();
		date_t date_cond_orders = {1995, 2, 15};
		date_t date_cond_lineitem = {1995, 3, 15};
		auto cond_orders = [&](const raw_orders &t)->bool{return t.o_orderdate < date_cond_orders;};
		auto cond_lineitem = [&](const raw_lineitem &t)->bool{return date_cond_lineitem < t.l_shipdate;};
		auto cond_customer = [&](const raw_customer &t)->bool{return strcmp(t.c_mktsegment, "BUILDING") == 0;};
		size_t STEP = 10000; //! 10ms step
		size_t MAX = 3000000; //! max time span
		cout<<"first orders:"<<endl;
		samplejoin<raw_customer, \
			raw_orders, \
			raw_lineitem, \
			std::function<std::string(const raw_customer &)>, \
			std::function<uint32_t(const raw_orders)>, \
			std::function<uint32_t(const raw_lineitem)>,\
			std::string,\
			uint32_t,\
			uint32_t,\
			std::function<double(const raw_lineitem &)> >
			(vec_c, \
				vec_o,\
				vec_l,\
				key_func_customer_mktsegment,\
				key_func_orders_custkey,\
				key_func_lineitem_orderkey,\
				key_func_customer_custkey,\
				key_func_orders_orderkey,\
				"BUILDING",
				cond_orders,
				cond_lineitem,
				result_func,
				STEP,
				MAX,
				.95);
		samplejoin_range_more<raw_lineitem, \
			raw_orders, \
			raw_customer, \
			std::function<date_t(const raw_lineitem &)>, \
			std::function<uint32_t(const raw_orders &)>, \
			std::function<uint32_t(const raw_customer &)>,\
			date_t,\
			uint32_t,\
			uint32_t,\
			std::function<double(const raw_lineitem &)> >
			(vec_l, \
				vec_o,\
				vec_c,\
				key_func_lineitem_shipdate,\
				key_func_orders_orderkey,\
				key_func_customer_custkey,\
				key_func_lineitem_orderkey,\
				key_func_orders_custkey,\
				date_cond_lineitem,\
				cond_orders, \
				cond_customer, \
				cmp_func_lineitem_shipdate,\
				result_func, \
				STEP, \
				MAX, \
				.95);
		samplejoin_range_less<raw_orders, \
			raw_customer, \
			raw_lineitem, \
			std::function<date_t(const raw_orders &)>, \
			std::function<uint32_t(const raw_customer &)>, \
			std::function<uint32_t(const raw_lineitem &)>,\
			date_t,\
			uint32_t,\
			uint32_t,\
			std::function<double(const raw_lineitem &)> >\
			(vec_o, \
				vec_c,\
				vec_l,\
				key_func_orders_orderdate,\
				key_func_customer_custkey,\
				key_func_lineitem_orderkey,\
				key_func_orders_custkey,\
				key_func_orders_orderkey,\
				date_cond_orders,
				cond_customer,
				cond_lineitem,
				cmp_func_orders_orderdate,
				result_func,
				STEP,
				MAX,
				.95);
};
void do_ripplejoin_no_cond(){
		auto  vec_c = load_customer();
		auto vec_o = load_orders();
		auto vec_l = load_lineitem();
		size_t STEP = 10000; //! 10ms step
		size_t MAX = 3000000; //! max time span
		auto result_func = [&](const raw_customer &t1, const raw_orders &t2, const raw_lineitem &t)->double{return t.l_extendedprice * (1.0 - t.l_discount);};
		ripple_simple<raw_customer, \
			raw_orders, \
			raw_lineitem, \
			uint32_t,\
			uint32_t,\
			uint32_t>\
			(vec_c, \
				vec_o,\
				vec_l,\
				key_func_customer_custkey,\
				key_func_orders_custkey,\
				key_func_lineitem_orderkey,\
				key_func_customer_custkey,\
				key_func_orders_orderkey,\
				result_func,
				STEP,
				MAX,
				.95);
}
void do_ripplejoin(){
		auto  vec_c = load_customer();
		auto vec_o = load_orders();
		auto vec_l = load_lineitem();
		date_t date_cond_orders = {1995, 2, 15};
		date_t date_cond_lineitem = {1995, 3, 15};
		auto cond_orders = [&](const raw_orders &t)->bool{return t.o_orderdate < date_cond_orders;};
		auto cond_lineitem = [&](const raw_lineitem &t)->bool{return date_cond_lineitem < t.l_shipdate;};
		auto cond_customer = [&](const raw_customer &t)->bool{return strcmp(t.c_mktsegment, "BUILDING") == 0;};
		auto result_func = [&](const raw_customer &t1, const raw_orders &t2, const raw_lineitem &t)->double{return t.l_extendedprice * (1.0 - t.l_discount);};
		size_t STEP = 10000; //! 10ms step
		size_t MAX = 3000000; //! max time span
		std::random_shuffle(vec_c.begin(), vec_c.end());
		std::random_shuffle(vec_o.begin(), vec_o.end());
		std::random_shuffle(vec_l.begin(), vec_l.end());
		ripple_simple<raw_customer, \
			raw_orders, \
			raw_lineitem, \
			uint32_t,\
			uint32_t,\
			uint32_t>\
			(vec_c, \
				vec_o,\
				vec_l,\
				key_func_customer_custkey,\
				key_func_orders_custkey,\
				key_func_lineitem_orderkey,\
				key_func_customer_custkey,\
				key_func_orders_orderkey,\
				cond_customer,
				cond_orders,
				cond_lineitem,
				result_func,
				STEP,
				MAX,
				.95);
}
int main(){
	//do_hashjoin();
	//do_samplejoin();
	do_ripplejoin_no_cond();
	return 0;
}
