#include "rawdata.h"

raw_region parse_raw_region(string &line){
	auto ele = tool::split(line, "|");
	raw_region tmp;
	tmp.r_regionkey = stoi(ele[0]);
	strcpy(tmp.r_name, ele[1].c_str());
	strcpy(tmp.r_comment, ele[2].c_str());
	return tmp;

};
raw_nation parse_raw_nation(string &line){
	auto ele = tool::split(line, "|");
	raw_nation tmp;
	tmp.n_nationkey = stoi(ele[0]);
	strcpy(tmp.n_name, ele[1].c_str());
	tmp.n_regionkey = stoi(ele[2]);
	strcpy(tmp.n_comment, ele[3].c_str());
	return tmp;
};

raw_supplier parse_raw_supplier(string &line){
	raw_supplier tmp;
	auto ele = tool::split(line, "|");
	tmp.s_suppkey = stoi(ele[0]);
	strcpy(tmp.s_name, ele[1].c_str());
	strcpy(tmp.s_address, ele[2].c_str());
	tmp.s_nationkey = stoi(ele[3]);
	strcpy(tmp.s_phone, ele[4].c_str());
	tmp.s_acctbal = stod(ele[5]);
	strcpy(tmp.s_comment, ele[6].c_str());
	return tmp;
};
raw_partsupp parse_raw_partsupp(string &line){
	auto ele = tool::split(line, "|");
	raw_partsupp tmp;
	tmp.ps_partkey = stoi(ele[0]);
	tmp.ps_suppkey = stoi(ele[1]);
	tmp.ps_availqty = stoi(ele[2]);
	tmp.ps_supplycost = stoi(ele[3]);
	strcpy(tmp.ps_comment, ele[4].c_str());
	return tmp;
};
raw_part parse_raw_part(string &line){
	auto ele = tool::split(line, "|");
	raw_part tmp;
	tmp.p_partkey = stoi(ele[0]);
	std::strncpy(tmp.p_name, ele[1].c_str(), 56);
	std::strncpy(tmp.p_mfgr, ele[2].c_str(), 26);
	std::strncpy(tmp.p_brand, ele[3].c_str(), 11);
	std::strncpy(tmp.p_type, ele[4].c_str(), 25);
	tmp.p_size = stoi(ele[5]);
	std::strncpy(tmp.p_container, ele[6].c_str(), 11);
	tmp.p_retailprice = stoi(ele[7]);
	std::strncpy(tmp.p_comment, ele[8].c_str(), 24);
	return tmp;
}
raw_customer parse_raw_customer(string &line){
	raw_customer cu;
	auto ele = tool::split(line, "|");
	cu.c_custkey = stoi(ele[0]);
	strncpy(cu.c_name, ele[1].c_str(), 26);
	strncpy(cu.c_address, ele[2].c_str(), 41);
	cu.c_nationkey = stoi(ele[3]);
	strncpy(cu.c_phone, ele[4].c_str(), 16);
	cu.c_acctbal = stod(ele[5]);
	strncpy(cu.c_mktsegment, ele[6].c_str(), 11);
	strncpy(cu.c_comment, ele[7].c_str(), 118);
	return cu;
};
raw_orders parse_raw_orders(string &line){
	auto ele = tool::split(line, "|");
	raw_orders o;
	o.o_orderkey = stoi(ele[0]);
	o.o_custkey = stoi(ele[1]);
	o.o_orderstatus = ele[2][0];
	o.o_totalprice = stod(ele[3]);
	o.o_orderdate = strtodate(ele[4]);
	strncpy(o.o_orderpriority, ele[5].c_str(), 16);
	strncpy(o.o_clerk, ele[6].c_str(), 16);
	o.o_shippriority = stoi(ele[7]);
	strncpy(o.o_comment, ele[8].c_str(), 80);
	return o;
};
raw_lineitem parse_raw_lineitem(string &line){
	raw_lineitem lt;
	auto ele = tool::split(line, "|");
	lt.l_orderkey = stoi(ele[0]);
	lt.l_partkey = stoi(ele[1]);
	lt.l_suppkey = stoi(ele[2]);
	lt.l_linenumber = stoi(ele[3]);
	lt.l_quantity = stod(ele[4]);
	lt.l_extendedprice = stod(ele[5]);
	lt.l_discount = stod(ele[6]);
	lt.l_tax = stod(ele[7]);
	lt.l_returnflag = ele[8][0];
	lt.l_linestatus = ele[9][0];
	lt.l_shipdate = strtodate(ele[10]);
	lt.l_commitdate = strtodate(ele[11]);
	lt.l_receiptdate = strtodate(ele[12]);
	strcpy(lt.l_shipinstruct, ele[13].c_str());
	strcpy(lt.l_shipmode,  ele[14].c_str());
	strcpy(lt.l_comment, ele[15].c_str());
	return lt;
};


raw_customer_orders operator+(const raw_customer &t1, const raw_orders &t2){
	raw_customer_orders r;
	r.uno = t1;
	r.dos = t2;
	return r;
}
raw_customer_orders_lineitem operator+(const raw_customer_orders &t1, const raw_lineitem &t2){
	raw_customer_orders_lineitem r;
	r.uno = t1;
	r.dos = t2;
	return r;
}
bool operator==(raw_region &t1, raw_nation &t2){
	return t1.r_regionkey == t2.n_regionkey;
}
bool operator==(raw_customer &t1, raw_orders &t2){
	return t1.c_custkey == t2.o_custkey;
}
bool operator==(raw_orders &t1, raw_lineitem &t2){
	return t1.o_orderkey == t2.l_orderkey;
}
bool operator==(raw_lineitem &t1, raw_orders &t2){
	return t1.l_orderkey == t2.o_orderkey;
}
bool operator==(raw_orders &t1, raw_customer &t2){
	return t1.o_custkey == t2.c_custkey;
}
//! key func for region
uint32_t key_func_region_regionkey(const raw_region &t){return t.r_regionkey;};
std::string key_func_region_regionname(const raw_region &t){return std::string(t.r_name);};
//! key func for nation
uint32_t key_func_nation_nationkey(const raw_nation &t){return t.n_nationkey;};
uint32_t key_func_nation_regionkey(const raw_nation &t){return t.n_regionkey;};
std::string key_func_nation_nationname(const raw_nation &t){return std::string(t.n_name);};
//! key func for partsupp
uint32_t key_func_partsupp_partkey(const raw_partsupp &t){return t.ps_partkey;};
uint32_t key_func_partsupp_suppkey(const raw_partsupp &t){return t.ps_suppkey;};
//! key func for supplier
uint32_t key_func_supplier_suppkey(const raw_supplier &t){return t.s_suppkey;};
uint32_t key_func_supplier_nationkey(const raw_supplier &t){return t.s_nationkey;};
//! key func for customer
uint32_t key_func_customer_custkey(const raw_customer &t ){return t.c_custkey;};
uint32_t key_func_customer_nationkey(const raw_customer &t){return t.c_nationkey;};
std::string key_func_customer_mktsegment(const raw_customer &t){return std::string(t.c_mktsegment);};
//! key func for orders
uint32_t key_func_orders_custkey(const raw_orders &t){return t.o_custkey;};
uint32_t key_func_orders_orderkey(const raw_orders &t){return t.o_orderkey;};
date_t key_func_orders_orderdate(const raw_orders &t){return t.o_orderdate;};
bool operator<(const date_t &x, raw_orders* y){return x < y->o_orderdate;};
bool operator<(raw_orders* y, const date_t &x){return y->o_orderdate < x;};
unsigned char key_func_orders_orderstatus(const raw_orders &t){return t.o_orderstatus;};
//! key func for lineitem
uint32_t key_func_lineitem_orderkey(const raw_lineitem &t){return t.l_orderkey;};
uint32_t key_func_lineitem_suppkey(const raw_lineitem &t){return t.l_suppkey;};
uint32_t key_func_lineitem_partkey(const raw_lineitem &t){return t.l_partkey;};
extern bool operator<(const date_t &x, raw_lineitem *y){return x <
y->l_shipdate;};
extern bool operator<(raw_lineitem *x, const date_t &y){return x->l_shipdate <
y;};
date_t key_func_lineitem_shipdate(const raw_lineitem &t){return t.l_shipdate;};
unsigned char key_func_lineitem_returnflag(const raw_lineitem &t){return t.l_returnflag;};



uint32_t key_func_customer_orders_orderkey(const raw_customer_orders &t){return t.dos.o_orderkey;};
bool cmp_func_lineitem_shipdate(const raw_lineitem &t1, const raw_lineitem &t2){
	return datetoi(t1.l_shipdate) < datetoi(t2.l_shipdate);
}
bool cmp_func_orders_orderdate(const raw_orders &t1, const raw_orders &t2){
	return datetoi(t1.o_orderdate) < datetoi(t2.o_orderdate);
}
