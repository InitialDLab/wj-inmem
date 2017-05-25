#ifndef _H_RAWDATA
#define _H_RAWDATA
#include <iostream>
#include <cctype>
#include <cstdint>
#include <cstring>
#include "tool.h"
#include "base.h"
using namespace std;


class raw_region:public base_raw{
	public:
	raw_region(){raw_tag = RAW_REGION;};
	uint32_t r_regionkey;
	char r_name[26];
	char r_comment[153];
	raw_region& operator=(const raw_region &t){
		std::strncpy(r_name, t.r_name, 26);
		std::strncpy(r_comment, t.r_comment, 153);
		r_regionkey = t.r_regionkey;
		pre = t.pre;
		return *this;
	}
};
class raw_nation:public base_raw{
	public:
		raw_nation(){raw_tag = RAW_NATION;};
		uint32_t n_nationkey;
		char n_name[25];
		uint32_t n_regionkey;
		char n_comment[152];
		raw_nation& operator=(const raw_nation& t){
			if(this == &t)
				return *this;
			n_nationkey = t.n_nationkey;
			std::strncpy(n_name, t.n_name, 25);
			std::strncpy(n_comment, t.n_comment, 152);
			pre = t.pre;
			return *this;
		}
};
class raw_supplier:public base_raw{
	public:
		raw_supplier(){raw_tag = RAW_SUPPLIER;};
		uint32_t s_suppkey;
		char s_name[26];
		char s_address[41];
		uint32_t s_nationkey;
		char s_phone[16];
		double s_acctbal;
		char s_comment[102];
		raw_supplier& operator=(const raw_supplier &t){
			if(&t == this)
				return *this;
			s_suppkey = t.s_suppkey;
			std::strncpy(s_name, t.s_name, 26);
			std::strncpy(s_address, t.s_address, 41);
			s_nationkey = t.s_nationkey;
			std::strncpy(s_phone, t.s_phone, 16);
			s_acctbal = t.s_acctbal;
			std::strncpy(s_comment, t.s_comment, 102);
			pre = t.pre;
			return *this;
		}
};
class raw_partsupp:public base_raw{
	public:
		raw_partsupp(){raw_tag = RAW_PARTSUPP;};
		uint32_t ps_partkey;
		uint32_t ps_suppkey;
		int32_t ps_availqty;
		double ps_supplycost;
		char ps_comment[200];
		raw_partsupp& operator=(const raw_partsupp &t){
			if(&t == this)
				return *this;
			ps_partkey = t.ps_partkey;
			ps_suppkey = t.ps_suppkey;
			ps_availqty = t.ps_availqty;
			ps_supplycost = t.ps_supplycost;
			std::strncpy(ps_comment, t.ps_comment, 200);
			pre = t.pre;
			return *this;
		}
};
class raw_customer:public base_raw{
	public:
		raw_customer(){raw_tag = RAW_CUSTOMER;};
		uint32_t c_custkey;
		char c_name[26];
		char c_address[41];
		uint32_t c_nationkey;
		char c_phone[16];
		double c_acctbal;
		char c_mktsegment[11];
		char c_comment[118];
		raw_customer& operator=(const raw_customer &t){
			if(this == &t)
				return *this;
			c_custkey = t.c_custkey;
			std::strncpy(c_name, t.c_name, 26);
			std::strncpy(c_address, t.c_address, 41);
			c_nationkey = t.c_nationkey;
			std::strncpy(c_phone, t.c_phone, 16);
			c_acctbal = t.c_acctbal;
			strncpy(c_mktsegment, t.c_mktsegment, 11);
			strncpy(c_comment, t.c_comment, 118);
			pre = t.pre;
			return *this;
		}
};

class raw_orders:public base_raw{
	public:
		raw_orders(){raw_tag = RAW_ORDERS;};
		uint32_t o_orderkey;
		uint32_t o_custkey;
		unsigned char o_orderstatus;
		double o_totalprice;
		date_t o_orderdate;
		char o_orderpriority[16];
		char o_clerk[16];
		int32_t o_shippriority;
		char o_comment[80];
		raw_orders& operator=(const raw_orders &t){
			if(this == &t)
				return *this;
			o_orderkey = t.o_orderkey;
			o_custkey = t.o_custkey;
			o_orderstatus = t.o_orderstatus;
			o_totalprice = t.o_totalprice;
			o_orderdate = t.o_orderdate;
			std::strncpy(o_orderpriority, t.o_orderpriority, 16);
			std::strncpy(o_clerk, t.o_clerk, 16);
			o_shippriority = t.o_shippriority;
			std::strncpy(o_comment, t.o_comment, 80);
			pre = t.pre;
			return *this;
		}
};
class raw_lineitem:public base_raw{
	public:
		raw_lineitem(){raw_tag = RAW_LINEITEM;}
		uint32_t l_orderkey;
		uint32_t l_partkey;
		uint32_t l_suppkey;
		int32_t l_linenumber;
		double l_quantity;
		double l_extendedprice;
		double l_discount;
		double l_tax;
		unsigned char l_returnflag;
		unsigned char l_linestatus;
		date_t l_shipdate;
		date_t l_commitdate;
		date_t l_receiptdate;
		char l_shipinstruct[26];
		char l_shipmode[11];
		char l_comment[45];
		raw_lineitem& operator=(const raw_lineitem &t){
			if(this == &t)
				return *this;
			l_orderkey = t.l_orderkey;
			l_partkey = t.l_partkey;
			l_suppkey = t. l_suppkey;
			l_linenumber = t.l_linenumber;
			l_quantity = t.l_quantity;
			l_extendedprice = t.l_extendedprice;
			l_discount = t.l_discount;
			l_tax = t.l_tax;
			l_returnflag = t.l_returnflag;
			l_linestatus = t.l_linestatus;
			l_shipdate = t.l_shipdate;
			l_commitdate = t.l_commitdate;
			l_receiptdate = t.l_receiptdate;
			std::strncpy(l_shipinstruct, t.l_shipinstruct, 26);
			std::strncpy(l_shipmode, t.l_shipmode, 11);
			std::strncpy(l_comment, t.l_comment, 45);
			pre = t.pre;
			return *this;
		}
};

class raw_part:public base_raw{
	public:
		raw_part(){raw_tag = RAW_PART;};
		uint32_t p_partkey;
		char p_name[56];
		char p_mfgr[26];
		char p_brand[11];
		char p_type[25];
		int32_t p_size;
		char p_container[11];
		double p_retailprice;
		char p_comment[24];
		raw_part& operator=(const raw_part &t){
			if(this == &t)
				return *this;
			p_partkey = t.p_partkey;
			std::strncpy(p_name, t.p_name, 56);
			std::strncpy(p_mfgr, t.p_mfgr, 26);
			std::strncpy(p_brand, t.p_brand, 11);
			std::strncpy(p_type, t.p_type, 25);
			p_size = t.p_size;
			std::strncpy(p_container, t.p_container, 11);
			p_retailprice = t.p_retailprice;
			std::strncpy(p_comment, t.p_comment, 24);
			pre = t.pre;
			return *this;
		}
};


struct raw_customer_orders{
	raw_customer uno;
	raw_orders dos;
	raw_customer_orders(){};
	raw_customer_orders(const raw_customer&t1, const raw_orders &t2){
		uno = t1;
		dos = t2;
	}
	raw_customer_orders& operator=(const raw_customer_orders &t){
		if(this == &t)
			return *this;
		uno = t.uno;
		dos = t.dos;
		return *this;
	}
};
extern raw_customer_orders operator+(const raw_customer &t1, const raw_orders &t2);

struct raw_customer_orders_lineitem{
	raw_customer_orders uno;
	raw_lineitem dos;
	raw_customer_orders_lineitem(){};
	raw_customer_orders_lineitem(const raw_customer_orders &t1, const raw_lineitem &t2){
		uno = t1;
		dos = t2;
	}
	raw_customer_orders_lineitem& operator=(const raw_customer_orders_lineitem &t){
		if(this == &t)
			return *this;
		uno = t.uno;
		dos = t.dos;
		return *this;
	}
};
extern raw_customer_orders_lineitem operator+(const raw_customer_orders &t1, const raw_lineitem &t2);

extern raw_region parse_raw_region(string &);
extern raw_nation parse_raw_nation(string &);
extern raw_supplier parse_raw_supplier(string &);
extern raw_partsupp parse_raw_partsupp(string &);
extern raw_customer parse_raw_customer(string &);
extern raw_orders parse_raw_orders(string &);
extern raw_lineitem parse_raw_lineitem(string &);
extern raw_part parse_raw_part(string &);

extern bool operator==(const raw_region &, const raw_nation &);
extern bool operator==(const raw_customer &, const raw_orders &);
extern bool operator==(const raw_orders &, const raw_lineitem &);
extern bool operator==(const raw_lineitem &, const raw_orders &);
extern bool operator==(const raw_orders &, const raw_customer &);
//! keyfunc for region
extern uint32_t key_func_region_regionkey(const raw_region &);
extern std::string key_func_region_regionname(const raw_region &);

//! keyfunc for nation
extern uint32_t key_func_nation_nationkey(const raw_nation &);
extern uint32_t key_func_nation_regionkey(const raw_nation &);
extern std::string key_func_nation_nationname(const raw_nation &);

//! key func for customer
extern uint32_t key_func_customer_custkey(const raw_customer &);
extern uint32_t key_func_customer_nationkey(const raw_customer &);
extern std::string key_func_customer_mktsegment(const raw_customer &);
//! key func for supplier
extern uint32_t key_func_supplier_suppkey(const raw_supplier &);
extern uint32_t key_func_supplier_nationkey(const raw_supplier &);
//! key func for orders
extern uint32_t key_func_orders_custkey(const raw_orders &);
extern uint32_t key_func_orders_orderkey(const raw_orders &);
extern unsigned char key_func_orders_orderstatus(const raw_orders &);
extern date_t key_func_orders_orderdate(const raw_orders &);
extern bool operator<(const date_t &, raw_orders*);
extern bool operator<(raw_orders*, const date_t &);
//! key func for lineitem
extern uint32_t key_func_lineitem_orderkey(const raw_lineitem &);
extern uint32_t key_func_lineitem_suppkey(const raw_lineitem &);
extern uint32_t key_func_lineitem_partkey(const raw_lineitem &);
extern date_t key_func_lineitem_shipdate(const raw_lineitem &);
extern unsigned char key_func_lineitem_returnflag(const raw_lineitem &);
extern bool operator<(const date_t &, raw_lineitem *);
extern bool operator<(raw_lineitem *, const date_t &);
//! key func for part
extern uint32_t key_func_part_partkey(const raw_part &);
extern std::string key_func_part_partname(const raw_part &);

//! key func for partsupp
extern uint32_t key_func_partsupp_partkey(const raw_partsupp &);
extern uint32_t key_func_partsupp_suppkey(const raw_partsupp &);

extern bool cmp_func_lineitem_shipdate(const raw_lineitem &, const raw_lineitem &);
extern uint32_t key_func_customer_orders_orderkey(const raw_customer_orders &);
extern bool cmp_func_orders_orderdate(const raw_orders &, const raw_orders &);
inline base_raw * increment(base_raw *p){
	if(p->raw_tag == RAW_NATION){
		auto q = dynamic_cast<raw_nation *>(p);
		return ++q;
	}
	if(p->raw_tag == RAW_REGION){
		auto q = dynamic_cast<raw_region *>(p);
		return ++q;
	}
	if(p->raw_tag == RAW_ORDERS){
		auto q = dynamic_cast<raw_orders *>(p);
		return ++q;
	}
	if(p->raw_tag == RAW_CUSTOMER){
		auto q = dynamic_cast<raw_customer *>(p);
		return ++q;
	}
	if(p->raw_tag == RAW_SUPPLIER){
		auto q = dynamic_cast<raw_supplier *>(p);
		return ++q;
	}
	if(p->raw_tag == RAW_LINEITEM){
		auto q = dynamic_cast<raw_lineitem *>(p);
		return ++q;
	}
	return ++p;
}
#endif
