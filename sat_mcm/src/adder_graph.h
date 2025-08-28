/* ADDER_GRAPH.h
 * Version:
 * Datum: 20.11.2014
 */

#ifndef ADDER_GRAPH_H
#define ADDER_GRAPH_H

#include <vector>
#include <list>
#include <map>
#include <stdint.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <typeinfo>
#include <climits>
#include <time.h>
#include <string>


using namespace std;

namespace PAGSuite
{

#define DONT_CARE LLONG_MAX

//definition of an abstract node:
  class node_t
  {
  public:
    int stage;                           //stage of the node, input nodes are always located in stage 0
//    std::vector<node_t *> inputs;         //references to input nodes, for common 2-input adders this vector has length two
    int latency;                         //latency of operation

    node_t()
    {
      stage = -1; //unknown stage
      latency = 0; //default latency is zero
    }

    virtual ~node_t()
    {}

    virtual void dummy(void)
    {}
  };

  class adder_graph_base_node_t : public node_t
  {
  public:                      //stage of the node, input nodes are always located in stage 0
    virtual ~adder_graph_base_node_t()
    {}

    std::vector<std::vector<int64_t> > output_factor;  // inner vector = inputs, outer vector = configurations
    std::map<std::string, std::string> specific_parameters;
    int output_right_shift = 0;

	// *************************转换为字符串的方法
    virtual std::string to_string() const {
        std::stringstream ss;
        ss << "Output Factor: [";
        for (const auto& row : output_factor) {
            ss << "[";
            for (const auto& value : row) {
                ss << value << ",";
            }
            ss.seekp(-1, ss.cur); // 去掉最后一个逗号
            ss << "],";
        }
        ss.seekp(-1, ss.cur); // 去掉最后一个逗号
        ss << "]\nSpecific Parameters: {";
        for (const auto& param : specific_parameters) {
            ss << param.first << ": " << param.second << ", ";
        }
        ss.seekp(-2, ss.cur); // 去掉最后一个逗号和空格
        ss << "}";
        return ss.str();
    }

    virtual void dummy(void)
    {}
  };

  class input_node_t : public adder_graph_base_node_t
  {
  public:
    ~input_node_t()
    {}

    input_node_t()
    {
      stage = 0; //input is stage 0
    }
  };

  class zero_node_t : public adder_graph_base_node_t
  {
  public:
    zero_node_t()
    {
      stage = 0; //stage is 0
    }
  };

//definition of an adder node which can represent an adder or subtractor:
  class adder_subtractor_node_t : public adder_graph_base_node_t
  {
  public:
    ~adder_subtractor_node_t()
    {}

    std::vector<adder_graph_base_node_t *> inputs;         //references to input nodes, for common 2-input adders this vector has length two
    std::vector<int64_t> input_shifts;       //shifts of inputs, left shifts are positive while right shifts negative, must have the same length than "inputs"
    std::vector<bool> input_is_negative; //false if the corresponding input is added, true if it is subtracted, must have the same length than "inputs", at least one element has to be false
  };

//definition of a multiplexer node:
  class mux_node_t : public adder_graph_base_node_t
  {
  public:
    ~mux_node_t()
    {}

    std::vector<adder_graph_base_node_t *> inputs;         //references to input nodes, length of this vector corresponds to the number of configurations
    std::vector<int64_t> input_shifts;       //shifts of inputs, left shifts are positive while right shifts negative, must have the same length as "inputs"
  };

  class conf_adder_subtractor_node_t : public adder_subtractor_node_t
  {
  public:
    ~conf_adder_subtractor_node_t()
    {}

    std::vector<std::vector<bool> > input_is_negative; //inner vector: sign of inputs of the adder (2 for the common adder), outer vector: configuration
  };

  class register_node_t : public adder_graph_base_node_t
  {
  public:
    ~register_node_t()
    {}

    adder_graph_base_node_t *input;
    int64_t input_shift;
  };

  class output_node_t : public register_node_t
  {
  public:
    ~output_node_t()
    {}
  };

  class adder_graph_t
  {
  public:
    adder_graph_t();

    ~adder_graph_t();

    //The adder graph itself is just a linked list of (internally) connected nodes:
    std::list<adder_graph_base_node_t *> nodes_list;
    std::map<std::string, std::string> specific_parameters;

    bool quiet;

    void pipeline_graph(); //pipeline adder_graph
    void drawdot(string filename = "./pag_adder_graph.dot"); //plot the adder graph to dot
    void writesyn(ostream &graphoutstream);

    void writesyn(string filename = "./pag_adder_graph.txt"); // write the graph down in mat syntax
    void print_graph(const std::vector<float>& constant);

    void get_cost(); //modify return type!
    void print_cost();

    void replace_node(adder_graph_base_node_t *original, adder_graph_base_node_t *replacement);

    bool parse_to_graph(string commandLine, bool ignore_outnodes = true);

    void add_input_nodes(int noOfConfigurations, int noOfInputs);

    bool normalize_node(adder_subtractor_node_t *node);  //normalize corrects order of signs of the node such that an adder's first input is always positive (if possible), returns true node could be normalized

    bool normalize_graph();  //normalize corrects order of signs of the graph such that an adder's first input is always positive (if possible), returns true when all nodes could be normalized

    void clear();

    void check_and_correct(string graphstring = ""); //check adder graph

    string get_adder_graph_as_string();

    static std::vector<std::vector<int64_t> > abs_vec(std::vector<std::vector<int64_t> > &factor); //computes the absolute value of a vector

    static bool factor_is_zero(std::vector<std::vector<int64_t> > &factor); //returns true if all elements are zero

    static std::vector<int64_t> normalize(std::vector<int64_t> &row);

    static std::vector<std::vector<int64_t> > normalize(std::vector<std::vector<int64_t> > &factor); //normalizes a vector such that the first element is positive

    static std::vector<std::vector<int64_t> >fundamental(std::vector<std::vector<int64_t> > &factor); //computes a fundamental of a vector

  private:
    adder_graph_base_node_t *parse_node(string nodeStr);

    bool parse_factor(string factorStr, std::vector<std::vector<int64_t> > *factor);

    string convert_old_syntax(string commandLine); //ASCII-Convert old syntax graph

    adder_graph_base_node_t *get_node_from_output_factor_in_stage(std::vector<std::vector<int64_t> > &output_factor, int stage);

    string matrix_to_string(const std::vector<std::vector<int64_t> > &matrix);

    string node_to_string(const adder_graph_base_node_t *node);

    bool input_nodes_have_been_added;

  };

  extern int computeWordSize(std::vector<std::vector<int64_t> > &output_factor, int wIn);

  extern int computeWordSize(adder_graph_base_node_t *node, int wIn);

  template<class T>
  inline bool is_a(const adder_graph_base_node_t &obj)
  {
    return typeid(T) == typeid(obj);
  }

  std::ostream &operator<<(std::ostream &stream, const std::vector<std::vector<int64_t> > &matrix);

}

std::ostream &operator<<(std::ostream &stream, const PAGSuite::adder_graph_t &adder_graph);


#endif // ADDER_GRAPH_H
