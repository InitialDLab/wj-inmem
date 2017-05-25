#ifndef _H_TABLE
#define _H_TABLE
#include <iostream>
#include <string>
#include "base.h"
#include "rawdata.h"
#include "loaddata.h"
class base_table{
	public:
		base_table(){
			table_tag = EMPTY;
			seed = std::chrono::system_clock::now().time_since_epoch().count();
			gen.seed(seed);
		}
		virtual ~base_table(){};
		virtual void load_data() = 0;
		virtual void build_hash() = 0;
		virtual base_raw* sample_all() = 0;
		virtual plan_result_type sample_range(head_node_type &input) = 0; 
		virtual plan_result_type search(base_raw *input) = 0;
		virtual plan_result_type search_cond(base_raw *input, cond_func_type cond_func) = 0;
		virtual size_t get_size(base_raw *input) = 0;
		virtual size_t size() = 0;
		table_tag_type table_tag;
	protected:
		void set_table_tag(const table_tag_type &t){
			table_tag = t;
		}
		table_tag_type get_table_tag(){
			return table_tag;
		}
		unsigned seed;
		std::mt19937 gen;
};
class table_region:public base_table{
	public:
		typedef raw_region raw_type;
		table_region(){
			table_tag = REGION;
			h_region_regionname = nullptr;
			h_region_regionkey = nullptr;
		}
		void load_data(){
			data = load_region();
		}
		void build_hash(){
			h_region_regionname = new hash_tree<raw_region, std::string>(data, key_func_region_regionname);
			h_region_regionkey = new hash_tree<raw_region, uint32_t>(data, key_func_region_regionkey);
		}
		
		~table_region(){
			delete h_region_regionname;
			delete h_region_regionkey;
		}
		base_raw *sample_all(){
			std::uniform_int_distribution<uint64_t> dis(0, data.size() - 1);
			return &data[dis(gen)];
		}
		plan_result_type sample_range(head_node_type &input){
			plan_result_type result;
			switch(input.col_id){
				case 1:
				{
					auto iter = h_region_regionkey->get_sample(input.key, result.d);
					result.p_raw = (base_raw*) iter;
					return result;
				}
				case 2:
				{
					auto iter = h_region_regionname->get_sample(input.segment, result.d);
					result.p_raw = (base_raw*) iter;
					return result;
				}
				default:
				return result;
			}
		}
		plan_result_type search(base_raw *input){
			auto p = (raw_nation *)input;
			plan_result_type result;
			auto iter = h_region_regionkey->get_sample(p->n_regionkey, result.d);
			result.p_raw = (base_raw*) iter;
			return result;
		}
		plan_result_type search_cond(base_raw *input, cond_func_type cond_func){
			auto p = (raw_nation *)input;
			plan_result_type result;
			auto iter = h_region_regionkey->get_sample(p->n_regionkey, result.d);
			result.p_raw = (base_raw*) iter;
			return result;
		}
		size_t get_size(base_raw *input){
			if(input == nullptr)
				return 0;
			if(input->raw_tag == RAW_NATION){
				auto p = (raw_nation *)input;
				return h_region_regionkey->get_size(p->n_regionkey);
			}
			return 0;
		}
		size_t size(){
			return h_region_regionname->size();
		}
		std::vector<raw_region> data;
	private:
		hash_tree<raw_region, std::string> *h_region_regionname;
		hash_tree<raw_region, uint32_t> *h_region_regionkey;
};
class table_part:public base_table{
	public:
		typedef raw_part raw_type;
		table_part(){
			table_tag = PART;
			h_part_partname = nullptr;
			h_part_partname = nullptr;
		}
		void load_data(){
			data = load_part();
		}
		void build_hash(){
			h_part_partname = new hash_tree<raw_part, std::string>(data, key_func_part_partname);
			h_part_partkey = new hash_tree<raw_part, uint32_t>(data, key_func_part_partkey);
		}
		
		~table_part(){
			delete h_part_partname;
			delete h_part_partkey;
		}
		base_raw *sample_all(){
			std::uniform_int_distribution<uint64_t> dis(0, data.size() - 1);
			return &data[dis(gen)];
		}
		plan_result_type sample_range(head_node_type &input){
			plan_result_type result;
			switch(input.col_id){
				case 1:
				{
					auto iter = h_part_partkey->get_sample(input.key, result.d);
					result.p_raw = (base_raw*) iter;
					return result;
				}
				case 2:
				{
					auto iter = h_part_partname->get_sample(input.segment, result.d);
					result.p_raw = (base_raw*) iter;
					return result;
				}
				default:
				return result;
			}
		}
		plan_result_type search(base_raw *input){
			auto p = (raw_partsupp *)input;
			plan_result_type result;
			auto iter = h_part_partkey->get_sample(p->ps_partkey, result.d);
			result.p_raw = (base_raw*) iter;
			return result;
		}
		plan_result_type search_cond(base_raw *input, cond_func_type cond_func){
			auto p = (raw_partsupp *)input;
			plan_result_type result;
			auto iter = h_part_partkey->get_sample(p->ps_partkey, result.d);
			result.p_raw = (base_raw*) iter;
			return result;
		}
		size_t get_size(base_raw *input){
			if(input == nullptr)
				return 0;
			return 0;
		}
		size_t size(){
			return h_part_partname->size();
		}
		std::vector<raw_part> data;
	private:
		hash_tree<raw_part, std::string> *h_part_partname;
		hash_tree<raw_part, uint32_t> *h_part_partkey;
};
class table_nation:public base_table{
	public:
		typedef raw_nation raw_type;
		table_nation(){
			table_tag = NATION;
			h_nation_nationname = nullptr;
			h_nation_nationkey = nullptr;
		};
		void load_data(){
			data = load_nation();
		}
		void build_hash(){
			h_nation_nationname = new hash_tree<raw_nation, std::string>(data, key_func_nation_nationname);
			h_nation_nationkey = new hash_tree<raw_nation, uint32_t>(data, key_func_nation_nationkey);
			h_nation_regionkey = new hash_tree<raw_nation, uint32_t>(data, key_func_nation_regionkey);
		}
		base_raw *sample_all(){
			std::uniform_int_distribution<uint64_t> dis(0, data.size() - 1);
			auto ind = dis(gen);
			return &data[ind];
		}
		plan_result_type sample_range(head_node_type &input){
			plan_result_type result;
			switch(input.col_id){
				case 1:
				{
					auto iter = h_nation_nationkey->get_sample(input.key, result.d);
					result.p_raw = (base_raw*) iter;
					return result;
				}
				case 2:
				{
					auto iter = h_nation_nationname->get_sample(input.segment, result.d);
					result.p_raw = (base_raw*) iter;
					return result;
				}
				case 3:
				{
					auto iter = h_nation_regionkey->get_sample(input.key, result.d);
					result.p_raw = (base_raw*) iter;
					return result;
				}
				default:
				{
					return result;
				}
			}
		}
		plan_result_type search(base_raw *input){
			plan_result_type result;
			if(input->raw_tag == RAW_EMPTY)
			{
				return result;
			}
			if(input->raw_tag == RAW_REGION)
			{
				auto p = (raw_region *)input;
				auto iter = h_nation_regionkey->get_sample(p->r_regionkey, result.d);
				result.p_raw = (base_raw*) iter;
				return result;
			}
			if(input->raw_tag == RAW_CUSTOMER)
			{
				auto p = (raw_customer *)input;
				auto iter = h_nation_nationkey->get_sample(p->c_nationkey, result.d);
				result.p_raw = (base_raw*) iter;
				return result;
			}
			if(input->raw_tag == RAW_SUPPLIER)
			{
				auto p = (raw_supplier *)input;
				auto iter = h_nation_nationkey->get_sample(p->s_nationkey, result.d);
				result.p_raw = (base_raw*) iter;
				return result;
			}
			return result;
		}
		plan_result_type search_cond(base_raw *input, cond_func_type cond_func){
			plan_result_type result;
			if(input->raw_tag == RAW_EMPTY)
			{
				return result;
			}
			if(input->raw_tag == RAW_REGION)
			{
				auto p = (raw_region *)input;
				std::vector<base_raw *> tmp;
				auto range = h_nation_regionkey->get_all(p->r_regionkey);
				base_raw *q = range.first;
				for(size_t i = 0; i < range.second; ++i){
					if((*cond_func)(q))
						tmp.push_back(q);
					q = increment(q);
				}
				if(tmp.size() == 0)
					return result;
				std::uniform_int_distribution<uint64_t> dis(0, tmp.size() - 1);
				auto ind = dis(gen);
				result.p_raw = tmp[ind];
				result.d = tmp.size();
				return result;
			}
			if(input->raw_tag == RAW_CUSTOMER)
			{
				auto p = (raw_customer *)input;
				std::vector<base_raw *> tmp;
				auto range = h_nation_nationkey->get_all(p->c_nationkey);
				base_raw *q = range.first;
				for(size_t i = 0; i < range.second; ++i){
					if((*cond_func)(q))
						tmp.push_back(q);
					q = increment(q);
				}
				if(tmp.size() == 0)
					return result;
				std::uniform_int_distribution<uint64_t> dis(0, tmp.size() - 1);
				auto ind = dis(gen);
				result.p_raw = tmp[ind];
				result.d = tmp.size();
				return result;
			}
			if(input->raw_tag == RAW_SUPPLIER)
			{
				auto p = (raw_supplier *)input;
				std::vector<base_raw *> tmp;
				auto range = h_nation_nationkey->get_all(p->s_nationkey);
				base_raw *q = range.first;
				for(size_t i = 0; i < range.second; ++i){
					if((*cond_func)(q))
						tmp.push_back(q);
					q = increment(q);
				}
				if(tmp.size() == 0)
					return result;
				std::uniform_int_distribution<uint64_t> dis(0, tmp.size() - 1);
				auto ind = dis(gen);
				result.p_raw = tmp[ind];
				result.d = tmp.size();
				return result;
			}
			return result;
		}
		size_t get_size(base_raw *input){
			if(input == nullptr)
				return 0;
			if(input->raw_tag == RAW_REGION)
			{
				auto p = (raw_region *)input;
				return h_nation_regionkey->get_size(p->r_regionkey);
			}
			if(input->raw_tag == RAW_CUSTOMER)
			{
				auto p = (raw_customer *)input;
				return h_nation_nationkey->get_size(p->c_nationkey);
			}
			if(input->raw_tag == RAW_SUPPLIER)
			{
				auto p = (raw_supplier *)input;
				return h_nation_nationkey->get_size(p->s_nationkey);
			}
			return 0;
		}
		~table_nation(){
			delete h_nation_nationname;
			delete h_nation_nationkey;
			delete h_nation_regionkey;
		}
		size_t size(){
			return h_nation_nationname->size();
		}
		std::vector<raw_nation> data;
	private:
		hash_tree<raw_nation, std::string> *h_nation_nationname;
		hash_tree<raw_nation, uint32_t> *h_nation_nationkey;
		hash_tree<raw_nation, uint32_t> *h_nation_regionkey;
};

class table_customer:public base_table{
	public:
		typedef raw_customer raw_type;
		table_customer(){
			table_tag = CUSTOMER;
			h_customer_custkey = nullptr;
			h_customer_nationkey = nullptr;
			h_customer_mktsegment = nullptr;
		};
		void load_data(){
			data = load_customer();
		};
		void build_hash(){
			h_customer_custkey = new hash_tree<raw_customer, uint32_t>(data, key_func_customer_custkey);
			h_customer_nationkey = new hash_tree<raw_customer, uint32_t>(data, key_func_customer_nationkey);
			h_customer_mktsegment = new hash_tree<raw_customer, std::string>(data, key_func_customer_mktsegment);
		};
		base_raw *sample_all(){
			std::uniform_int_distribution<uint64_t> dis(0, data.size() - 1);
			return &data[dis(gen)];
		}
		plan_result_type sample_range(head_node_type &input){
			plan_result_type result;
			switch(input.col_id){
				case 1:
				{
					auto iter = h_customer_custkey->get_sample(input.key, result.d);
					result.p_raw = (base_raw*) iter;
					return result;
				}
				case 2:
				{
					auto iter = h_customer_nationkey->get_sample(input.key, result.d);
					result.p_raw = (base_raw*) iter;
					return result;
				}
				case 3:
				{
					auto iter = h_customer_mktsegment->get_sample(input.segment, result.d);
					result.p_raw = (base_raw*) iter;
					return result;
				}
				default:
				{
					return result;
				}
			}
		}
		plan_result_type search_cond(base_raw *input, cond_func_type cond_func){
			plan_result_type result;
			if(input->raw_tag == RAW_NATION)
			{
				auto p = (raw_nation *)input;
				auto range = h_customer_nationkey->get_all(p->n_nationkey);
				std::vector<base_raw *> tmp;
				base_raw *q = range.first;
				for(size_t i = 0; i < range.second; ++i){
					if((*cond_func)(q))
						tmp.push_back(q);
					q = increment(q);
				}
				if(tmp.size() == 0)
					return result;
				std::uniform_int_distribution<uint64_t> dis(0, tmp.size() - 1);
				auto ind = dis(gen);
				result.p_raw = tmp[ind];
				result.d = tmp.size();
				return result;
			}
			if(input->raw_tag == RAW_ORDERS)
			{
				auto p = (raw_orders *)input;
				auto range = h_customer_custkey->get_all(p->o_custkey);
				std::vector<base_raw *> tmp;
				base_raw *q = range.first;
				for(size_t i = 0; i < range.second; ++i){
					if((*cond_func)(q))
						tmp.push_back(q);
					q = increment(q);
				}
				if(tmp.size() == 0)
					return result;
				std::uniform_int_distribution<uint64_t> dis(0, tmp.size() - 1);
				auto ind = dis(gen);
				result.p_raw = tmp[ind];
				result.d = tmp.size();
				return result;
			}
			return result;
		}
	plan_result_type search(base_raw *input){
			plan_result_type result;
			if(input->raw_tag == RAW_NATION)
			{
				auto p = (raw_nation *)input;
				auto iter = h_customer_nationkey->get_sample(p->n_nationkey, result.d);
				result.p_raw = (base_raw*) iter;
				return result;
			}
			if(input->raw_tag == RAW_ORDERS)
			{
				auto p = (raw_orders *)input;
				auto iter = h_customer_custkey->get_sample(p->o_custkey, result.d);
				result.p_raw = (base_raw*) iter;
				return result;
			}
			return result;
			}
		size_t get_size(base_raw *input){
			if(input == nullptr)
				return 0;
			plan_result_type result;
			if(input->raw_tag == RAW_NATION)
			{
				auto p = (raw_nation *)input;
				return h_customer_nationkey->get_size(p->n_nationkey);
			}
			if(input->raw_tag == RAW_ORDERS)
			{
				auto p = (raw_orders *)input;
				return h_customer_custkey->get_size(p->o_custkey);
			}
			return 0;
		}
		size_t size(){
			return h_customer_custkey->size();
		}
		std::vector<raw_customer> data;
	private:
		hash_tree<raw_customer, uint32_t> *h_customer_custkey;
		hash_tree<raw_customer, uint32_t> *h_customer_nationkey;
		hash_tree<raw_customer, std::string> *h_customer_mktsegment;
};
class table_partsupp:public base_table{
	public:
		typedef raw_partsupp raw_type;
		table_partsupp(){
			table_tag = PARTSUPP;
			h_partsupp_suppkey = nullptr;
			h_partsupp_partkey = nullptr;
		};
		void load_data(){
			data = load_partsupp();
		}
		void build_hash(){
			h_partsupp_partkey = new hash_tree<raw_partsupp, uint32_t>(data, key_func_partsupp_partkey);
			h_partsupp_suppkey = new hash_tree<raw_partsupp, uint32_t>(data, key_func_partsupp_suppkey);
		}
		base_raw *sample_all(){
			std::uniform_int_distribution<uint64_t> dis(0, data.size() - 1);
			return &data[dis(gen)];
		}
		~table_partsupp(){
			delete h_partsupp_suppkey;
			delete h_partsupp_partkey;
		}
		plan_result_type sample_range(head_node_type &input){
			plan_result_type result;
			switch(input.col_id){
				case 1:
				{
					auto iter = h_partsupp_partkey->get_sample(input.key, result.d);
					result.p_raw = (base_raw*) iter;
					return result;
				}
				case 2:
				{
					auto iter = h_partsupp_suppkey->get_sample(input.key, result.d);
					result.p_raw = (base_raw*) iter;
					return result;
				}
				default:
				{
					return result;
				}
			}
			return result;
		}
		plan_result_type search(base_raw *input){
			plan_result_type result;
			if(input->raw_tag == RAW_PART)
			{
				auto p = (raw_part *)input;
				auto iter = h_partsupp_partkey->get_sample(p->p_partkey, result.d);
				result.p_raw = (base_raw*) iter;
				return result;
			}
			if(input->raw_tag == RAW_SUPPLIER)
			{
				auto p = (raw_supplier *)input;
				auto iter = h_partsupp_suppkey->get_sample(p->s_suppkey, result.d);
				result.p_raw = (base_raw*) iter;
				return result;
			}
			if(input->raw_tag == RAW_LINEITEM)
			{
				auto p = (raw_lineitem *)input;
				auto iter = h_partsupp_suppkey->get_sample(p->l_suppkey, result.d);
				result.p_raw = (base_raw*) iter;
				return result;
			}
			return result;
		}
		plan_result_type search_all(base_raw *input){
			plan_result_type result;
			if(input == nullptr)
				return result;
			if(input == nullptr)
				return result;
			if(input->raw_tag == RAW_PART)
			{
				auto p = (raw_part *)input;
				auto iter = h_partsupp_partkey->get_sample(p->p_partkey, result.d);
				result.p_raw = (base_raw*) iter;
				return result;
			}
			if(input->raw_tag == RAW_SUPPLIER)
			{
				auto p = (raw_supplier *)input;
				auto iter = h_partsupp_suppkey->get_sample(p->s_suppkey, result.d);
				result.p_raw = (base_raw*) iter;
				return result;
			}
			if(input->raw_tag == RAW_LINEITEM)
			{
				auto p = (raw_lineitem *)input;
				auto iter = h_partsupp_suppkey->get_sample(p->l_suppkey, result.d);
				result.p_raw = (base_raw*) iter;
				return result;
			}
			return result;
		}

		plan_result_type search_cond(base_raw *input, cond_func_type cond_func){
            plan_result_type result;
            return result;
        }

		size_t get_size(base_raw *input){
			if(input == nullptr)
				return 0;
			if(input->raw_tag == RAW_PART)
			{
				auto p = (raw_part *)input;
				return h_partsupp_partkey->get_size(p->p_partkey);
			}
			if(input->raw_tag == RAW_SUPPLIER)
			{
				auto p = (raw_supplier *)input;
				return h_partsupp_suppkey->get_size(p->s_suppkey);
			}
			if(input->raw_tag == RAW_LINEITEM)
			{
				auto p = (raw_lineitem *)input;
				return h_partsupp_suppkey->get_size(p->l_suppkey);
			}
			return 0;
		}
		std::vector<raw_partsupp> data;
		size_t size(){
			return h_partsupp_suppkey->size();
		}
	private:
		hash_tree<raw_partsupp, uint32_t> *h_partsupp_suppkey;
		hash_tree<raw_partsupp, uint32_t> *h_partsupp_partkey;
};

class table_supplier:public base_table{
	public:
		typedef raw_supplier raw_type;
		table_supplier(){
			table_tag = SUPPLIER;
			h_supplier_suppkey = nullptr;
			h_supplier_nationkey = nullptr;
		};
		void load_data(){
			data = load_supplier();
		}
		void build_hash(){
			h_supplier_suppkey = new hash_tree<raw_supplier, uint32_t>(data, key_func_supplier_suppkey);
			h_supplier_nationkey = new hash_tree<raw_supplier, uint32_t>(data, key_func_supplier_nationkey);
		}
		base_raw *sample_all(){
			std::uniform_int_distribution<uint64_t> dis(0, data.size() - 1);
			return &data[dis(gen)];
		}
		~table_supplier(){
			delete h_supplier_suppkey;
			delete h_supplier_nationkey;
		}
		plan_result_type sample_range(head_node_type &input){
			plan_result_type result;
			switch(input.col_id){
				case 1:
				{
					auto iter = h_supplier_suppkey->get_sample(input.key, result.d);
					result.p_raw = (base_raw*) iter;
					return result;
				}
				case 2:
				{
					auto iter = h_supplier_nationkey->get_sample(input.key, result.d);
					result.p_raw = (base_raw*) iter;
					return result;
				}
				default:
				{
					return result;
				}
			}
			return result;
		}
		plan_result_type search(base_raw *input){
			plan_result_type result;
			if(input->raw_tag == RAW_NATION)
			{
				auto p = (raw_nation *)input;
				auto iter = h_supplier_nationkey->get_sample(p->n_nationkey, result.d);
				result.p_raw = (base_raw*) iter;
				return result;
			}
			if(input->raw_tag == RAW_PARTSUPP)
			{
				auto p = (raw_partsupp *)input;
				auto iter = h_supplier_suppkey->get_sample(p->ps_suppkey, result.d);
				result.p_raw = (base_raw*) iter;
				return result;
			}
			if(input->raw_tag == RAW_LINEITEM)
			{
				auto p = (raw_lineitem *)input;
				auto iter = h_supplier_suppkey->get_sample(p->l_suppkey, result.d);
				result.p_raw = (base_raw*) iter;
				return result;
			}
			return result;
		}
		plan_result_type search_cond(base_raw *input, cond_func_type cond_func){
			plan_result_type result;
			if(input->raw_tag == RAW_NATION)
			{
				auto p = (raw_nation *)input;
				auto range = h_supplier_nationkey->get_all(p->n_nationkey);
				std::vector<base_raw *> tmp;
				base_raw *q = range.first;
				for(size_t i = 0; i < range.second; ++i){
					if((*cond_func)(q))
						tmp.push_back(q);
					q = increment(q);
				}
				if(tmp.size() == 0)
					return result;
				std::uniform_int_distribution<uint64_t> dis(0, tmp.size() - 1);
				auto ind = dis(gen);
				result.p_raw = tmp[ind];
				result.d = tmp.size();
				return result;
			}
			if(input->raw_tag == RAW_PARTSUPP)
			{
				auto p = (raw_partsupp *)input;
				auto range = h_supplier_suppkey->get_all(p->ps_suppkey);
				std::vector<base_raw *> tmp;
				base_raw *q = range.first;
				for(size_t i = 0; i < range.second; ++i){
					if((*cond_func)(q))
						tmp.push_back(q);
					q = increment(q);
				}
				if(tmp.size() == 0)
					return result;
				std::uniform_int_distribution<uint64_t> dis(0, tmp.size() - 1);
				auto ind = dis(gen);
				result.p_raw = tmp[ind];
				result.d = tmp.size();
				return result;
			}
			if(input->raw_tag == RAW_LINEITEM)
			{
				auto p = (raw_lineitem *)input;
				auto range = h_supplier_suppkey->get_all(p->l_suppkey);
				std::vector<base_raw *> tmp;
				base_raw *q = range.first;
				for(size_t i = 0; i < range.second; ++i){
					if((*cond_func)(q))
						tmp.push_back(q);
					q = increment(q);
				}
				if(tmp.size() == 0)
					return result;
				std::uniform_int_distribution<uint64_t> dis(0, tmp.size() - 1);
				auto ind = dis(gen);
				result.p_raw = tmp[ind];
				result.d = tmp.size();
				return result;
			}
			return result;
		}
		size_t get_size(base_raw *input){
			if(input == nullptr)
				return 0;
			if(input->raw_tag == RAW_NATION)
			{
				auto p = (raw_nation *)input;
				return h_supplier_nationkey->get_size(p->n_nationkey);
			}
			if(input->raw_tag == RAW_PARTSUPP)
			{
				auto p = (raw_partsupp *)input;
				return h_supplier_suppkey->get_size(p->ps_suppkey);
			}
			if(input->raw_tag == RAW_LINEITEM)
			{
				auto p = (raw_lineitem *)input;
				return h_supplier_suppkey->get_size(p->l_suppkey);
			}
			return 0;
		}
		std::vector<raw_supplier> data;
		size_t size(){
			return h_supplier_suppkey->size();
		}
	private:
		hash_tree<raw_supplier, uint32_t> *h_supplier_suppkey;
		hash_tree<raw_supplier, uint32_t> *h_supplier_nationkey;
};

class table_orders:public base_table{
	public:
		typedef raw_orders raw_type;
		table_orders(){
			table_tag = ORDERS;
			h_orders_orderkey = nullptr;
			h_orders_custkey = nullptr;
			h_orders_orderdate = nullptr;
			h_orders_orderstatus = nullptr;
		};
		void load_data(){
			data = load_orders();
		}
		void build_hash(){
			h_orders_orderkey = new hash_tree<raw_orders, uint32_t>(data, key_func_orders_orderkey);
			h_orders_custkey = new hash_tree<raw_orders, uint32_t>(data, key_func_orders_custkey);
			h_orders_orderstatus = new hash_tree<raw_orders, unsigned char>(data, key_func_orders_orderstatus);
			h_orders_orderdate = new range_tree<raw_orders, date_t>(data, key_func_orders_orderdate);
		}
		~table_orders(){
			delete h_orders_orderkey;
			delete h_orders_custkey;
			delete h_orders_orderdate;
			delete h_orders_orderstatus;
		}
		base_raw *sample_all(){
			std::uniform_int_distribution<uint64_t> dis(0, data.size() - 1);
			return &data[dis(gen)];
		}
		plan_result_type sample_range(head_node_type &input){
			plan_result_type result;
			switch(input.col_id){
				case 1:
				{
					auto iter = h_orders_orderkey->get_sample(input.key, result.d);
					result.p_raw = (base_raw*) iter;
					return result;
				}
				case 2:
				{
					auto iter = h_orders_custkey->get_sample(input.key, result.d);
					result.p_raw = (base_raw*) iter;
					return result;
				}
				case 3:
				{
					auto iter = h_orders_orderdate->get_sample(input.start_date, input.end_date, result.d);
					result.p_raw = (base_raw*) iter;
					return result;
				}
				case 4:
				{
					auto iter = h_orders_orderstatus->get_sample(input.flag, result.d);
					result.p_raw = (base_raw*) iter;
					return result;
				}
				default:
				{
					return result;
				}
			}
		}
		plan_result_type search(base_raw *input){
			plan_result_type result;
			if(input->raw_tag == RAW_CUSTOMER)
			{
				auto p = (raw_customer *)input;
				auto iter = h_orders_custkey->get_sample(p->c_custkey, result.d);
				result.p_raw = (base_raw*) iter;
				return result;
			}
			if(input->raw_tag == RAW_LINEITEM)
			{
				auto p = (raw_lineitem *)input;
				auto iter = h_orders_orderkey->get_sample(p->l_orderkey, result.d);
				result.p_raw = (base_raw*) iter;
				return result;
			}
			return result;
		}
		plan_result_type search_cond(base_raw *input, cond_func_type cond_func){
			plan_result_type result;
			if(input->raw_tag == RAW_CUSTOMER)
			{
				auto p = static_cast<raw_customer *>(input);
				auto range = h_orders_custkey->get_all(p->c_custkey);
				if(range.first == nullptr)
					return result;
				std::vector<base_raw *> tmp;
				base_raw *q = range.first;
				//std::cout<<"all: ";
				for(size_t i = 0; i < range.second; ++i){
					//auto debug_q = (raw_orders*) q;
					//std::cout<<debug_q->o_orderdate<<" ";
					if((*cond_func)(q))
						tmp.push_back(q);
					q = increment(q);
				}
				//std::cout<<std::endl;
				//std::cout<<range.second<<"|"<<tmp.size()<<std::endl;
				//std::cout<< h_orders_custkey->get_size(p->c_custkey)<<std::endl;;
				if(tmp.empty())
					return result;
				std::uniform_int_distribution<uint64_t> dis(0, tmp.size() - 1);
				size_t ind = dis(gen);
				result.p_raw = tmp[ind];
				result.d = tmp.size();
				return result;
			}
			if(input->raw_tag == RAW_LINEITEM)
			{
				auto p = (raw_lineitem *)input;
				auto range = h_orders_orderkey->get_all(p->l_orderkey);
				std::vector<base_raw *> tmp;
				base_raw *q = range.first;
				for(size_t i = 0; i < range.second; ++i){
					if((*cond_func)(q))
						tmp.push_back(q);
					q = increment(q);
				}
				if(tmp.size() == 0)
					return result;
				std::uniform_int_distribution<uint64_t> dis(0, tmp.size() - 1);
				auto ind = dis(gen);
				result.p_raw = tmp[ind];
				result.d = tmp.size();
				return result;
			}
			return result;
		}
		size_t get_size(base_raw *input){
			if(input == nullptr)
				return 0;
			if(input->raw_tag == RAW_CUSTOMER)
			{
				auto p = (raw_customer *)input;
				return h_orders_custkey->get_size(p->c_custkey);
			}
			if(input->raw_tag == RAW_LINEITEM)
			{
				auto p = (raw_lineitem *)input;
				return h_orders_orderkey->get_size(p->l_orderkey);
			}
			return 0;
		}
		std::vector<raw_orders> data;
		size_t size(){
			return h_orders_orderkey->size();
		}
	private:
		hash_tree<raw_orders, uint32_t> *h_orders_orderkey;
		hash_tree<raw_orders, uint32_t> *h_orders_custkey;
		hash_tree<raw_orders, unsigned char> *h_orders_orderstatus;
		range_tree<raw_orders, date_t> *h_orders_orderdate;
};

class table_lineitem:public base_table{
	public:
		typedef raw_lineitem raw_type;
		table_lineitem(){
			table_tag = LINEITEM;
			h_lineitem_orderkey = nullptr;
			h_lineitem_suppkey = nullptr;
			h_lineitem_partkey = nullptr;
			h_lineitem_shipdate = nullptr;
			h_lineitem_returnflag = nullptr;
		}
		void load_data(){
			data = load_lineitem();
		}
		void build_hash(){
			h_lineitem_orderkey = new hash_tree<raw_lineitem, uint32_t>(data, key_func_lineitem_orderkey);
			h_lineitem_suppkey = new hash_tree<raw_lineitem, uint32_t>(data, key_func_lineitem_suppkey);
			h_lineitem_partkey = new hash_tree<raw_lineitem, uint32_t>(data, key_func_lineitem_partkey);
			h_lineitem_shipdate = new range_tree<raw_lineitem, date_t>(data, key_func_lineitem_shipdate);
			h_lineitem_returnflag = new hash_tree<raw_lineitem, unsigned char>(data, key_func_lineitem_returnflag);
		}
		base_raw *sample_all(){
			std::uniform_int_distribution<uint64_t> dis(0, data.size() - 1);
			auto ind = dis(gen);
			return &data[ind];
		}
		~table_lineitem(){
			delete h_lineitem_orderkey;
			delete h_lineitem_suppkey;
			delete h_lineitem_partkey;
			delete h_lineitem_shipdate;
			delete h_lineitem_returnflag;
		}
		plan_result_type sample_range(head_node_type &input){
			plan_result_type result;
			switch(input.col_id){
				case 1:
				{
					auto iter = h_lineitem_orderkey->get_sample(input.key, result.d);
					result.p_raw = (base_raw*) iter;
					return result;
				}
				case 2:
				{
					auto iter = h_lineitem_partkey->get_sample(input.key, result.d);
					result.p_raw = (base_raw*) iter;
					return result;
				}
				case 3:
				{
					auto iter = h_lineitem_suppkey->get_sample(input.key, result.d);
					result.p_raw = (base_raw*) iter;
					return result;
				}
				case 4:
				{
					auto iter = h_lineitem_returnflag->get_sample(input.flag, result.d);
					result.p_raw = (base_raw*) iter;
					return result;
				}
				case 5:
				{
					auto iter = h_lineitem_shipdate->get_sample(input.start_date, input.end_date, result.d);
					result.p_raw = (base_raw*) iter;
					return result;
				}
				default:
				{
					return result;
				}
			}
		}
		plan_result_type search(base_raw *input){
			plan_result_type result;
			if(input->raw_tag == RAW_PARTSUPP)
			{
				auto p = (raw_partsupp *)input;
				auto iter = h_lineitem_partkey->get_sample(p->ps_partkey, result.d);
				result.p_raw = (base_raw*) iter;
				return result;
			}
			if(input->raw_tag ==  RAW_ORDERS)
			{
				auto p = (raw_orders *)input;
				auto iter = h_lineitem_orderkey->get_sample(p->o_orderkey, result.d);
				result.p_raw = (base_raw*) iter;
				return result;
			}
			if(input->raw_tag == RAW_SUPPLIER)
			{
				auto p = (raw_supplier *)input;
				auto iter = h_lineitem_suppkey->get_sample(p->s_suppkey, result.d);
				result.p_raw = (base_raw*) iter;
				return result;
			}
			return result;
		}
		plan_result_type search_cond(base_raw *input, cond_func_type cond_func){
			plan_result_type result;
			if(input->raw_tag == RAW_PARTSUPP)
			{
				auto p = (raw_partsupp *)input;
				auto range = h_lineitem_partkey->get_all(p->ps_partkey);
				std::vector<base_raw *> tmp;
				base_raw *q = range.first;
				for(size_t i = 0; i < range.second; ++i){
					if((*cond_func)(q))
						tmp.push_back(q);
					q = increment(q);
				}
				if(tmp.size() == 0)
					return result;
				std::uniform_int_distribution<uint64_t> dis(0, tmp.size() - 1);
				auto ind = dis(gen);
				result.p_raw = tmp[ind];
				result.d = tmp.size();
				return result;
			}
			if(input->raw_tag ==  RAW_ORDERS)
			{
				auto p = (raw_orders *)input;
				auto range = h_lineitem_orderkey->get_all(p->o_orderkey);
				std::vector<base_raw *> tmp;
				base_raw *q = range.first;
				for(size_t i = 0; i < range.second; ++i){
					if((*cond_func)(q))
						tmp.push_back(q);
					q = increment(q);
				}
				if(tmp.size() == 0)
					return result;
				std::uniform_int_distribution<uint64_t> dis(0, tmp.size() - 1);
				auto ind = dis(gen);
				result.p_raw = tmp[ind];
				result.d = tmp.size();
				return result;
			}
			if(input->raw_tag == RAW_SUPPLIER)
			{
				auto p = (raw_supplier *)input;
				auto range = h_lineitem_suppkey->get_all(p->s_suppkey);
				std::vector<base_raw *> tmp;
				base_raw *q = range.first;
				for(size_t i = 0; i < range.second; ++i){
					if((*cond_func)(q))
						tmp.push_back(q);
					q = increment(q);
				}
				if(tmp.size() == 0)
					return result;
				std::uniform_int_distribution<uint64_t> dis(0, tmp.size() - 1);
				auto ind = dis(gen);
				result.p_raw = tmp[ind];
				result.d = tmp.size();
				return result;
			}
			return result;
		}
		size_t get_size(base_raw *input){
			if(input == nullptr)
				return 0;
			if(input->raw_tag == RAW_PARTSUPP)
			{
				auto p = (raw_partsupp *)input;
				return h_lineitem_partkey->get_size(p->ps_partkey);
			}
			if(input->raw_tag ==  RAW_ORDERS)
			{
				auto p = (raw_orders *)input;
				return h_lineitem_orderkey->get_size(p->o_orderkey);
			}
			if(input->raw_tag == RAW_SUPPLIER)
			{
				auto p = (raw_supplier *)input;
				return h_lineitem_suppkey->get_size(p->s_suppkey);
			}
			return 0;
		}
		plan_result_type search_all(base_raw *input){
			plan_result_type result;
			if(input == nullptr)
				return result;
			result.raw_type = RAW_LINEITEM;
			if(input->raw_tag == RAW_PARTSUPP)
			{
				auto p = (raw_partsupp *)input;
				auto tmp = h_lineitem_partkey->get_all(p->ps_partkey);
				result.p_raw = tmp.first;
				result.len = tmp.second;
				result.d = tmp.second;
				return result;
			}
			if(input->raw_tag ==  RAW_ORDERS)
			{
				auto p = (raw_orders *)input;
				auto tmp = h_lineitem_orderkey->get_all(p->o_orderkey);
				result.p_raw = tmp.first;
				result.len = tmp.second;
				result.d = tmp.second;
				return result;
			}
			if(input->raw_tag == RAW_SUPPLIER)
			{
				auto p = (raw_supplier *)input;
				auto tmp = h_lineitem_suppkey->get_all(p->s_suppkey);
				result.p_raw = tmp.first;
				result.len = tmp.second;
				result.d = tmp.second;
				return result;
			}
			return result;
		}
		size_t size(){
			return h_lineitem_orderkey->size();
		}
		std::vector<raw_lineitem> data;
	private:
		hash_tree<raw_lineitem, uint32_t> *h_lineitem_orderkey;
		hash_tree<raw_lineitem, uint32_t> *h_lineitem_suppkey;
		hash_tree<raw_lineitem, uint32_t> *h_lineitem_partkey;
		hash_tree<raw_lineitem, unsigned char> *h_lineitem_returnflag;
		range_tree<raw_lineitem, date_t> *h_lineitem_shipdate;
};
#endif
