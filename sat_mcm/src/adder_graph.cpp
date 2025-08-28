/* ADDER_GRAPH.cpp
 * Version:
 * Datum: 20.11.2014
 */

#include "adder_graph.h"
//#include "paglib_copa.h"
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <assert.h>
#include <sstream>
#include <variant> 
#include <iomanip> 
#include <unordered_set>
// 声明全局变量
extern int num_matrix;

namespace PAGSuite
{

  adder_graph_t::adder_graph_t()
  {
    quiet = true;
  }


  adder_graph_t::~adder_graph_t()
  {

  }

  void adder_graph_t::clear()
  {
    while (!nodes_list.empty())
    {
      delete nodes_list.back();
      nodes_list.pop_back();
    }
  }

/* pipeline the adder graph KONRAD*/
  void adder_graph_t::pipeline_graph()
  {
    adder_graph_t mygraph = *this;
    std::map<pair<vector<vector<int64_t> >, int>, adder_graph_base_node_t *> tmp_reg_list;

    for (std::list<adder_graph_base_node_t *>::iterator it = mygraph.nodes_list.begin(), it_end = mygraph.nodes_list.end();
         it != it_end; ++it)
    {
      if (is_a<adder_subtractor_node_t>(*(*it)))
      {
        adder_subtractor_node_t *t = (adder_subtractor_node_t *) (*it); //
        for (int k = 0, k_end = (int) t->inputs.size(); k < k_end; ++k)
        {
          int diff = (*it)->stage - t->inputs[k]->stage; //Differenz der Stufen
          for (int j = 1; j < diff; ++j)
          {
            register_node_t *tmp = new register_node_t;
            tmp->stage = t->stage - j;
            tmp->output_factor = t->inputs[k]->output_factor;
            if (j == 1)
            {
              if (tmp_reg_list[pair<vector<vector<int64_t> >, int>(t->inputs[k]->output_factor, tmp->stage)] ==
                  NULL) //Register existiert noch nicht
              {
                tmp->input = t->inputs[k];
                t->inputs[k] = tmp;
                tmp_reg_list[pair<vector<vector<int64_t> >, int>(t->inputs[k]->output_factor, tmp->stage)] = tmp;
              }
              else
              {
                t->inputs[k] = tmp_reg_list[pair<vector<vector<int64_t> >, int>(t->inputs[k]->output_factor, tmp->stage)];
              }
            }
            else if (j == diff - 1)
            {

              if (tmp_reg_list[pair<vector<vector<int64_t> >, int>(t->inputs[k]->output_factor, tmp->stage)] == NULL)
              {
                tmp->input = ((register_node_t *) (tmp_reg_list[pair<vector<vector<int64_t> >, int>(
                  t->inputs[k]->output_factor, (tmp->stage) + 1)]))->input;
                tmp_reg_list[pair<vector<vector<int64_t> >, int>(t->inputs[k]->output_factor, tmp->stage)] = tmp;
                ((register_node_t *) (tmp_reg_list[pair<vector<vector<int64_t> >, int>(t->inputs[k]->output_factor,
                                                                                       (tmp->stage) + 1)]))->input = tmp;
              }
              else
              {
                ((register_node_t *) (tmp_reg_list[pair<vector<vector<int64_t> >, int>(t->inputs[k]->output_factor,
                                                                                       (tmp->stage) +
                                                                                       1)]))->input = tmp_reg_list[pair<vector<vector<int64_t> >, int>(
                  t->inputs[k]->output_factor, tmp->stage)];
              }
            }
            else
            {
              if (tmp_reg_list[pair<vector<vector<int64_t> >, int>(t->inputs[k]->output_factor, tmp->stage)] == NULL)
              {
                tmp->input = ((register_node_t *) (tmp_reg_list[pair<vector<vector<int64_t> >, int>(
                  t->inputs[k]->output_factor, (tmp->stage) + 1)]))->input;
                tmp_reg_list[pair<vector<vector<int64_t> >, int>(t->inputs[k]->output_factor, tmp->stage)] = tmp;
                ((register_node_t *) (tmp_reg_list[pair<vector<vector<int64_t> >, int>(t->inputs[k]->output_factor,
                                                                                       (tmp->stage) + 1)]))->input = tmp;
              }
              else
              {
                ((register_node_t *) (tmp_reg_list[pair<vector<vector<int64_t> >, int>(t->inputs[k]->output_factor,
                                                                                       (tmp->stage) +
                                                                                       1)]))->input = tmp_reg_list[pair<vector<vector<int64_t> >, int>(
                  t->inputs[k]->output_factor, tmp->stage)];
              }
            }
          }
        }
      }
      else if (is_a<mux_node_t>(*(*it)))
      {
        mux_node_t *t = (mux_node_t *) (*it); //
        for (int k = 0, k_end = (int) t->inputs.size(); k < k_end; ++k)
        {
          if (t->inputs[k] != NULL)
          {
            int diff = (*it)->stage - t->inputs[k]->stage; //Differenz der Stufen
            for (int j = 1; j < diff; ++j)
            {
              register_node_t *tmp = new register_node_t;
              tmp->stage = t->stage - j;
              tmp->output_factor = t->inputs[k]->output_factor;
              if (j == 1)
              {
                if (tmp_reg_list[pair<vector<vector<int64_t> >, int>(t->inputs[k]->output_factor, tmp->stage)] ==
                    NULL) //Register existiert noch nicht
                {
                  tmp->input = t->inputs[k];
                  t->inputs[k] = tmp;
                  tmp_reg_list[pair<vector<vector<int64_t> >, int>(t->inputs[k]->output_factor, tmp->stage)] = tmp;
                }
                else
                {
                  t->inputs[k] = tmp_reg_list[pair<vector<vector<int64_t> >, int>(t->inputs[k]->output_factor,
                                                                                  tmp->stage)];
                }
              }
              else if (j == diff - 1)
              {

                if (tmp_reg_list[pair<vector<vector<int64_t> >, int>(t->inputs[k]->output_factor, tmp->stage)] == NULL)
                {
                  tmp->input = ((register_node_t *) (tmp_reg_list[pair<vector<vector<int64_t> >, int>(
                    t->inputs[k]->output_factor, (tmp->stage) + 1)]))->input;
                  tmp_reg_list[pair<vector<vector<int64_t> >, int>(t->inputs[k]->output_factor, tmp->stage)] = tmp;
                  ((register_node_t *) (tmp_reg_list[pair<vector<vector<int64_t> >, int>(t->inputs[k]->output_factor,
                                                                                         (tmp->stage) +
                                                                                         1)]))->input = tmp;
                }
                else
                {
                  ((register_node_t *) (tmp_reg_list[pair<vector<vector<int64_t> >, int>(t->inputs[k]->output_factor,
                                                                                         (tmp->stage) +
                                                                                         1)]))->input = tmp_reg_list[pair<vector<vector<int64_t> >, int>(
                    t->inputs[k]->output_factor, tmp->stage)];
                }
              }
              else
              {
                if (tmp_reg_list[pair<vector<vector<int64_t> >, int>(t->inputs[k]->output_factor, tmp->stage)] == NULL)
                {
                  tmp->input = ((register_node_t *) (tmp_reg_list[pair<vector<vector<int64_t> >, int>(
                    t->inputs[k]->output_factor, (tmp->stage) + 1)]))->input;
                  tmp_reg_list[pair<vector<vector<int64_t> >, int>(t->inputs[k]->output_factor, tmp->stage)] = tmp;
                  ((register_node_t *) (tmp_reg_list[pair<vector<vector<int64_t> >, int>(t->inputs[k]->output_factor,
                                                                                         (tmp->stage) +
                                                                                         1)]))->input = tmp;
                }
                else
                {
                  ((register_node_t *) (tmp_reg_list[pair<vector<vector<int64_t> >, int>(t->inputs[k]->output_factor,
                                                                                         (tmp->stage) +
                                                                                         1)]))->input = tmp_reg_list[pair<vector<vector<int64_t> >, int>(
                    t->inputs[k]->output_factor, tmp->stage)];
                }
              }
            }
          }
        }
      }
      else if (is_a<conf_adder_subtractor_node_t>(*(*it)))
      {
        conf_adder_subtractor_node_t *t = (conf_adder_subtractor_node_t *) (*it); //conf_adder_subtractor_node_t*
        for (int k = 0, k_end = (int) t->inputs.size(); k < k_end; ++k)
        {
          int diff = (*it)->stage - t->inputs[k]->stage; //Differenz der Stufen
          for (int j = 1; j < diff; ++j)
          {
            register_node_t *tmp = new register_node_t;
            tmp->stage = t->stage - j;
            tmp->output_factor = t->inputs[k]->output_factor;
            if (j == 1)
            {
              if (tmp_reg_list[pair<vector<vector<int64_t> >, int>(t->inputs[k]->output_factor, tmp->stage)] ==
                  NULL) //Register existiert noch nicht
              {
                tmp->input = t->inputs[k];
                t->inputs[k] = tmp;
                tmp_reg_list[pair<vector<vector<int64_t> >, int>(t->inputs[k]->output_factor, tmp->stage)] = tmp;
              }
              else
              {
                t->inputs[k] = tmp_reg_list[pair<vector<vector<int64_t> >, int>(t->inputs[k]->output_factor, tmp->stage)];
              }
            }
            else if (j == diff - 1)
            {

              if (tmp_reg_list[pair<vector<vector<int64_t> >, int>(t->inputs[k]->output_factor, tmp->stage)] == NULL)
              {
                tmp->input = ((register_node_t *) (tmp_reg_list[pair<vector<vector<int64_t> >, int>(
                  t->inputs[k]->output_factor, (tmp->stage) + 1)]))->input;
                tmp_reg_list[pair<vector<vector<int64_t> >, int>(t->inputs[k]->output_factor, tmp->stage)] = tmp;
                ((register_node_t *) (tmp_reg_list[pair<vector<vector<int64_t> >, int>(t->inputs[k]->output_factor,
                                                                                       (tmp->stage) + 1)]))->input = tmp;
              }
              else
              {
                ((register_node_t *) (tmp_reg_list[pair<vector<vector<int64_t> >, int>(t->inputs[k]->output_factor,
                                                                                       (tmp->stage) +
                                                                                       1)]))->input = tmp_reg_list[pair<vector<vector<int64_t> >, int>(
                  t->inputs[k]->output_factor, tmp->stage)];
              }
            }
            else
            {
              if (tmp_reg_list[pair<vector<vector<int64_t> >, int>(t->inputs[k]->output_factor, tmp->stage)] == NULL)
              {
                tmp->input = ((register_node_t *) (tmp_reg_list[pair<vector<vector<int64_t> >, int>(
                  t->inputs[k]->output_factor, (tmp->stage) + 1)]))->input;
                tmp_reg_list[pair<vector<vector<int64_t> >, int>(t->inputs[k]->output_factor, tmp->stage)] = tmp;
                ((register_node_t *) (tmp_reg_list[pair<vector<vector<int64_t> >, int>(t->inputs[k]->output_factor,
                                                                                       (tmp->stage) + 1)]))->input = tmp;
              }
              else
              {
                ((register_node_t *) (tmp_reg_list[pair<vector<vector<int64_t> >, int>(t->inputs[k]->output_factor,
                                                                                       (tmp->stage) +
                                                                                       1)]))->input = tmp_reg_list[pair<vector<vector<int64_t> >, int>(
                  t->inputs[k]->output_factor, tmp->stage)];
              }
            }
          }
        }

      }
    }

    for (std::map<pair<vector<vector<int64_t> >, int>, adder_graph_base_node_t *>::iterator iter = tmp_reg_list.begin(), iter_end = tmp_reg_list.end();
         iter != iter_end; ++iter)
    {
      mygraph.nodes_list.push_back((*iter).second);
    }
    std::cout << "Finished Pipelining" << std::endl;
  }

  std::vector<std::vector<int64_t> > adder_graph_t::abs_vec(std::vector<std::vector<int64_t> > &factor)
  {
    if (factor.size() == 0) return factor;

    std::vector<std::vector<int64_t> > factor_abs(factor.size(), std::vector<int64_t>(factor[0].size()));

    for (int r = 0; r < factor.size(); r++)
    {
      for (int c = 0; c < factor[r].size(); c++)
      {
        factor_abs[r][c] = abs(factor[r][c]);
      }
    }

    return factor_abs;
  }

  bool adder_graph_t::factor_is_zero(std::vector<std::vector<int64_t> > &factor)
  {
    for (int r = 0; r < factor.size(); r++)
    {
      for (int c = 0; c < factor[r].size(); c++)
      {
        if(factor[r][c] != 0) return false;
      }
    }

    return true;
  }


  std::vector<int64_t> adder_graph_t::normalize(std::vector<int64_t> &row)
  {
    bool found_leading_nonzero = false;

    std::vector<int64_t> row_norm(row.size());

    for (int c = 0; c < row.size(); c++)
    {
      if (!found_leading_nonzero)
      {
        row_norm[c] = row[c];
        if (row[c] != 0)
        {
          found_leading_nonzero = true;
          if (row[c] > 0)
          {
            return row;
          }
          else
          {
            row_norm[c] = -row[c];
          }
        }
      }
      else
      {
        //now, row needs normalization
        row_norm[c] = -row[c];
      }
    }
    return row_norm;
  }

  std::vector<std::vector<int64_t> > adder_graph_t::fundamental(std::vector<std::vector<int64_t> > &matrix)
  {
    if (matrix.size() == 0) return matrix;

    std::vector<std::vector<int64_t> > fundamental(matrix.size(), std::vector<int64_t>(matrix[0].size()));
    fundamental = matrix; //copy matrix

    while(true)
    {
      //check if there is any odd number, return if true
      for(int r = 0; r < fundamental.size(); r++)
      {
        for (int c = 0; c < fundamental[r].size(); c++)
        {
          if((fundamental[r][c] % 2) != 0)
            return fundamental; //we found one element that is odd, so we are done!
        }
      }

      //this part can only be reached if all elements are even, let's shift them one to the right and check again
      for(int r = 0; r < fundamental.size(); r++)
      {
        for (int c = 0; c < fundamental[r].size(); c++)
        {
          fundamental[r][c] >>= 1;
        }
      }
    }

  }

  std::vector<std::vector<int64_t> > adder_graph_t::normalize(std::vector<std::vector<int64_t> > &factor)
  {
    if (factor.size() == 0) return factor;

    std::vector<std::vector<int64_t> > factor_norm(factor.size(), std::vector<int64_t>(factor[0].size()));

    for (int r = 0; r < factor.size(); r++)
    {
      factor_norm[r] = normalize(factor[r]);
    }

    return factor_norm;
  }

  string adder_graph_t::matrix_to_string(const std::vector<std::vector<int64_t> > &matrix)
  {
    std::stringstream stream;

    for (int r = 0; r < matrix.size(); r++)
    {
      for (int c = 0; c < matrix[r].size(); c++)
      {
        if (matrix[r][c] < 0) stream << "m";

        stream << abs(matrix[r][c]);
        if (c < matrix[r].size() - 1) stream << "_";
      }
      if (r < matrix.size() - 1) stream << "__";
    }
    return stream.str();
  }

  string adder_graph_t::node_to_string(const adder_graph_base_node_t *node)
  {
    std::stringstream stream;

    if (is_a<input_node_t>(*node))
    {
      stream << "I";
    }
    else if ((is_a<adder_subtractor_node_t>(*node)) || (is_a<conf_adder_subtractor_node_t>(*node)))
    {
      stream << "A";
    }
    else if (is_a<register_node_t>(*node))
    {
      stream << "R";
    }
    else if (is_a<output_node_t>(*node))
    {
      stream << "O";
    }
    else if (is_a<mux_node_t>(*node))
    {
      stream << "M";
    }
    else if (is_a<zero_node_t>(*node))
    {
      stream << "Z"; //zero node!
    }
    else if (is_a<adder_graph_base_node_t>(*node))
    {
      stream << "U"; //unknown node!
    }
    else
    {
      throw runtime_error("node_to_string: unknown node of type " + string(typeid(*node).name()));
    }
    stream << matrix_to_string(node->output_factor) << "_s" << node->stage;
    return stream.str();
  }

  adder_graph_base_node_t *adder_graph_t::get_node_from_output_factor_in_stage(std::vector<std::vector<int64_t> > &output_factor, int stage)
  {
    std::vector<std::vector<int64_t> > output_factor_norm = normalize(output_factor);
    if (!quiet) cout << "searching for node with factor " << output_factor << " (normalized:" << output_factor_norm << ") in stage " << stage << endl;

    adder_graph_base_node_t *node_found = nullptr;
    for (adder_graph_base_node_t *node: nodes_list)
    {
      if(is_a<output_node_t>(*node))
      {
        if (!quiet) cout << "  skipping node " << node << " with factor " << node->output_factor << " as it is an output node" << endl;
        continue;
      }
      std::vector<std::vector<int64_t> > node_output_factor_norm = normalize(node->output_factor);
      if (!quiet) cout << "  checking node " << node << " with factor " << node->output_factor << " (normalized " << node_output_factor_norm << ") and stage=" << node->stage << endl;
      if ((node_output_factor_norm == output_factor_norm) && (node->stage == stage))
      {
        if (!quiet) cout << "  found node " << node << endl;
        if(node->output_factor == output_factor)
        {
          //signs are identical, so we can be sure that this is the right node
          return node;
        }
        else
        {
          //signs are different, we keep searching if there is a node with identical sign, if not we take this one
          node_found = node;
        }
      }
    }
    if (!quiet && node_found == nullptr) cout << "  node not found" << endl;
    return node_found;
  }

/* plot the adder graph to dot KONRAD*/
  void adder_graph_t::drawdot(string filename)
  {
    adder_graph_t mygraph = *this;
    //string filename = "./k_adder_graph.dot";

    FILE *graphfilepointer = NULL;

    graphfilepointer = fopen(filename.c_str(), "w");

    time_t mytime;
    time(&mytime);
    fprintf(graphfilepointer, "// File generated:  %s\n", ctime(&mytime));
    fprintf(graphfilepointer, "digraph DAG {\n");
    std::map<string, int> realized_edges;

    // plot nodes
    for (std::list<adder_graph_base_node_t *>::iterator it = mygraph.nodes_list.begin(), it_end = mygraph.nodes_list.end(); it != it_end; ++it)
    {
      stringstream node_string;
      node_string << node_to_string((*it)) << "[label=\"";
      node_string << (*it)->output_factor;
      double height = 0.3 * (*it)->output_factor.size();
      node_string << "\",fontsize=12,shape=";
      if (is_a<input_node_t>(*(*it)))
      {
        node_string << "ellipse];\n";
      }
      else if (is_a<adder_subtractor_node_t>(*(*it)))
      {
//      node_string << "box,height = " << height << ",width=.4];\n";
        node_string << "box];\n";
      }
      else if (is_a<mux_node_t>(*(*it)))
      {
        node_string << "polygon,sides=4,distortion=.7,fixedsize=\"true\",height = " << height - 0.3 * height << ",width=1.2];\n";
      }
      else if (is_a<conf_adder_subtractor_node_t>(*(*it)))
      {
        node_string << "ellipse,height =" << height << ",width=.4];\n";
      }
      else if (is_a<register_node_t>(*(*it)))
      {
        node_string << "box];\n";
      }
      else if (is_a<output_node_t>(*(*it)) || is_a<zero_node_t>(*(*it)) || is_a<adder_graph_base_node_t>(*(*it)))
      {
        node_string << "none];\n";
      }
      else
      {
        adder_graph_base_node_t *p = *it;
        throw runtime_error("unknown node of type " + string(typeid(*p).name()));
      }
      fprintf(graphfilepointer, "%s", node_string.str().c_str());
    }

    for (std::list<adder_graph_base_node_t *>::iterator it = mygraph.nodes_list.begin(), it_end = mygraph.nodes_list.end(); it != it_end; ++it)
    {
      if (is_a<adder_subtractor_node_t>(*(*it)))
      {
        for (int i = 0, i_end = (int) ((adder_subtractor_node_t *) (*it))->inputs.size(); i < i_end; ++i)
        {
          stringstream edge_string;
          edge_string << node_to_string(((adder_subtractor_node_t *) (*it))->inputs[i]) << " -> " << node_to_string((*it));
          edge_string << " [label=\"";
          if (((adder_subtractor_node_t *) (*it))->input_shifts[i] != 0)
            edge_string << ((adder_subtractor_node_t *) (*it))->input_shifts[i];
          if (((adder_subtractor_node_t *) (*it))->input_is_negative[i])
            edge_string << "(-)";
          edge_string << "\",fontsize=12]\n";
          if (realized_edges.find(edge_string.str().c_str()) == realized_edges.end())
          {
            realized_edges[edge_string.str().c_str()] = 1;
            fprintf(graphfilepointer, "%s", edge_string.str().c_str());
          }
        }
      }
      else if (is_a<mux_node_t>(*(*it)))
      {
        for (int i = 0, i_end = (int) ((mux_node_t *) (*it))->inputs.size(); i < i_end; ++i)
        {
          stringstream edge_string;
          if (((mux_node_t *) (*it))->inputs[i] != 0)
            edge_string << node_to_string(((mux_node_t *) (*it))->inputs[i]) << " -> " << node_to_string(*it);

          if (((mux_node_t *) (*it))->input_shifts[i] != 0 && ((mux_node_t *) (*it))->input_shifts[i] != DONT_CARE)
          {
            edge_string << " [label=\"";
            edge_string << ((mux_node_t *) (*it))->input_shifts[i];
            edge_string << "\",fontsize=12]\n";
          }
          else
            edge_string << "\n";
          if (realized_edges.find(edge_string.str().c_str()) == realized_edges.end())
          {
            realized_edges[edge_string.str().c_str()] = 1;
            fprintf(graphfilepointer, "%s", edge_string.str().c_str());
          }
        }
      }
      else if (is_a<conf_adder_subtractor_node_t>(*(*it)))
      {
        for (int i = 0, i_end = (int) ((conf_adder_subtractor_node_t *) (*it))->inputs.size(); i < i_end; ++i)
        {
          stringstream edge_string;
          edge_string << node_to_string(((conf_adder_subtractor_node_t *) (*it))->inputs[i]) << " -> " << node_to_string((*it));
          edge_string << " [label=\"";
          if (((conf_adder_subtractor_node_t *) (*it))->input_shifts[i] != 0)
            edge_string << ((conf_adder_subtractor_node_t *) (*it))->input_shifts[i] << "\\n";
          for (int j = 0, j_end = (int) ((conf_adder_subtractor_node_t *) (*it))->input_is_negative.size(); j < j_end; ++j)
          {
            if (((conf_adder_subtractor_node_t *) (*it))->input_is_negative[j][i])
              edge_string << "-" << "\\n";
            else
              edge_string << "+" << "\\n";
          }
          edge_string << "\",fontsize=12]\n";
          if (realized_edges.find(edge_string.str().c_str()) == realized_edges.end())
          {
            realized_edges[edge_string.str().c_str()] = 1;
            fprintf(graphfilepointer, "%s", edge_string.str().c_str());
          }
        }
      }
      else if (is_a<register_node_t>(*(*it)) || is_a<output_node_t>(*(*it)))
      {
        stringstream edge_string;
        edge_string << node_to_string(((register_node_t *) (*it))->input) << " -> " << node_to_string(*it);
        if (((register_node_t *) (*it))->input_shift > 0)
        {
          edge_string << "[label=\"" << ((register_node_t *) (*it))->input_shift << "\",fontsize=12]";
        }
        edge_string << "\n";
        if (realized_edges.find(edge_string.str().c_str()) == realized_edges.end())
        {
          realized_edges[edge_string.str().c_str()] = 1;
          fprintf(graphfilepointer, "%s", edge_string.str().c_str());
        }
      }
      else if (is_a<input_node_t>(*(*it)) || is_a<zero_node_t>(*(*it)))
      {
        //do nothing
      }
      else if (is_a<adder_graph_base_node_t>(*(*it)))
      {
        //do nothing
        //if (!quiet)
        std::cout << "WARNING: unspecified node (adder_graph_base_node_t)!" << std::endl;
      }
      else
      {
        throw runtime_error("node type not supported");
      }

    }

    fprintf(graphfilepointer, "}\n");
    fflush(graphfilepointer);
    fclose(graphfilepointer);
    if (!quiet) std::cout << "DOT-File generated" << std::endl;
  }

/* write the graph down in mat syntax KONRAD*/

  void adder_graph_t::writesyn(string filename)
  { //, string filename
    ofstream graphfilestream;
    graphfilestream.open(filename.c_str(), ofstream::out | ofstream::trunc);
    writesyn(graphfilestream);
    graphfilestream.close();
  }

  void adder_graph_t::writesyn(ostream &graphoutstream)
  {
    adder_graph_t mygraph = *this;

    int current_id = -1;
    time_t mytime;
    time(&mytime);

    graphoutstream << "{";
    // plot nodes
    for (std::list<adder_graph_base_node_t *>::iterator it = mygraph.nodes_list.begin(), it_end = mygraph.nodes_list.end(); it != it_end; ++it)
    {

      //node_t* tmp = (node_t*)(*it);
      if (!is_a<input_node_t>(*(*it)))
      {
        current_id++;
        if (current_id > 0) graphoutstream << ",";

        stringstream node_string;
        stringstream pre_node_string;
        stringstream type_string;
        type_string << "{";
        node_string << "[";
        for (int i = 0, i_end = (int) (*it)->output_factor.size(); i < i_end; ++i)
        {
          for (int j = 0, j_end = (int) (*it)->output_factor[i].size(); j < j_end; ++j)
          {
            if (j == (int) (*it)->output_factor[i].size() - 1)
            {
              if ((*it)->output_factor[i][j] != DONT_CARE)
                node_string << (*it)->output_factor[i][j];
              else
                node_string << "NaN";
            }
            else
            {
              if ((*it)->output_factor[i][j] != DONT_CARE)
                node_string << (*it)->output_factor[i][j] << ",";
              else
                node_string << "NaN,";
            }

          }
          if (i < (int) (*it)->output_factor.size() - 1)
            node_string << ";";
        }
        node_string << "]," << (*it)->stage << ",";

        if (is_a<adder_subtractor_node_t>(*(*it)) || is_a<conf_adder_subtractor_node_t>(*(*it)))
        {
          type_string << "'A',";
          for (int k = 0, k_end = (int) ((adder_subtractor_node_t *) (*it))->inputs.size(); k < k_end; ++k)
          {
            if (k == 0)
              pre_node_string << "[";
            else
              pre_node_string << ",[";
            //set outputs of pre_node
            for (int i = 0, i_end = (int) ((adder_subtractor_node_t *) (*it))->inputs[k]->output_factor.size();
                 i < i_end; ++i)
            {
              for (int j = 0, j_end = (int) ((adder_subtractor_node_t *) (*it))->inputs[k]->output_factor[i].size();
                   j < j_end; ++j)
              {
                if (j == (int) ((adder_subtractor_node_t *) (*it))->inputs[k]->output_factor[i].size() - 1)
                {
                  if (((adder_subtractor_node_t *) (*it))->inputs[k]->output_factor[i][j] != DONT_CARE)
                  {
                    bool negate = false;
                    if (is_a<adder_subtractor_node_t>(*(*it)))
                    {
                      if (((adder_subtractor_node_t *) (*it))->input_is_negative[k])
                        negate = true;
                    }
                    else
                    {
                      if (((conf_adder_subtractor_node_t *) (*it))->input_is_negative[i][k])
                        negate = true;
                    }
                    if (negate)
                      pre_node_string << (-1 * ((adder_subtractor_node_t *) (*it))->inputs[k]->output_factor[i][j]);
                    else pre_node_string << ((adder_subtractor_node_t *) (*it))->inputs[k]->output_factor[i][j];
                  }
                  else
                    pre_node_string << "NaN";
                }
                else
                {
                  if (((adder_subtractor_node_t *) (*it))->inputs[k]->output_factor[i][j] != DONT_CARE)
                  {
                    bool negate = false;
                    if (is_a<adder_subtractor_node_t>(*(*it)))
                    {
                      if (((adder_subtractor_node_t *) (*it))->input_is_negative[k])
                        negate = true;
                    }
                    else
                    {
                      if (((conf_adder_subtractor_node_t *) (*it))->input_is_negative[i][k])
                        negate = true;
                    }
                    if (negate)
                      pre_node_string << (-1 * ((adder_subtractor_node_t *) (*it))->inputs[k]->output_factor[i][j])
                                      << ",";
                    else pre_node_string << ((adder_subtractor_node_t *) (*it))->inputs[k]->output_factor[i][j] << ",";
                  }
                  else
                    pre_node_string << "NaN,";
                }

              }
              if (i < (int) (*it)->output_factor.size() - 1)
                pre_node_string << ";";
            }
            //get input shift
            pre_node_string << "]," << ((adder_subtractor_node_t *) (*it))->inputs[k]->stage << ","
                            << ((adder_subtractor_node_t *) (*it))->input_shifts[k];
          }
          //}
        }
        if (is_a<mux_node_t>(*(*it)))
        {
          type_string << "'M',";
          map<adder_graph_base_node_t *, bool> already_seen;
          already_seen.clear();

          bool k_first = true;
          for (int k = 0, k_end = (int) ((mux_node_t *) (*it))->inputs.size(); k < k_end; ++k)
          {
            if (already_seen.find(((mux_node_t *) (*it))->inputs[k]) == already_seen.end() &&
                ((mux_node_t *) (*it))->inputs[k] != NULL)
            {
              already_seen[((mux_node_t *) (*it))->inputs[k]] = true;
              if (k_first)
              {
                pre_node_string << "[";
                k_first = false;
              }
              else
                pre_node_string << ",[";
              //get outputs of pre_node
              for (int i = 0, i_end = (int) ((mux_node_t *) (*it))->inputs[k]->output_factor.size(); i < i_end; ++i)
              {
                for (int j = 0, j_end = (int) ((mux_node_t *) (*it))->inputs[k]->output_factor[i].size(); j < j_end; ++j)
                {
                  if (j == (int) ((mux_node_t *) (*it))->inputs[k]->output_factor[i].size() - 1)
                  {
                    if (((mux_node_t *) (*it))->inputs[k]->output_factor[i][j] != DONT_CARE)
                    {
                      pre_node_string << ((mux_node_t *) (*it))->inputs[k]->output_factor[i][j];
                    }
                    else
                    {
                      pre_node_string << "NaN";
                    }
                  }
                  else
                  {
                    if (((mux_node_t *) (*it))->inputs[k]->output_factor[i][j] != DONT_CARE)
                    {
                      pre_node_string << ((mux_node_t *) (*it))->inputs[k]->output_factor[i][j] << ",";
                    }
                    else
                    {
                      pre_node_string << "NaN,";
                    }
                  }
                }
                if (i < (int) (*it)->output_factor.size() - 1)
                  pre_node_string << ";";
              }
              pre_node_string << "]," << ((mux_node_t *) (*it))->inputs[k]->stage << ",[";
              //get input shifts
              for (int i = 0, i_end = (int) ((mux_node_t *) (*it))->input_shifts.size(); i < i_end; ++i)
              {

                if (((mux_node_t *) (*it))->input_shifts[k] != DONT_CARE && ((mux_node_t *) (*it))->inputs[k] != NULL)
                {
                  bool same = true;
                  for (unsigned int j = 0, j_end = ((mux_node_t *) (*it))->output_factor[i].size(); j < j_end; ++j)
                  {
                    if (((mux_node_t *) (*it))->output_factor[i][j] !=
                        ((mux_node_t *) (*it))->inputs[k]->output_factor[i][j] << ((mux_node_t *) (*it))->input_shifts[i])
                      same = false;
                  }

                  if (!same)
                    pre_node_string << "NaN";
                  else
                    pre_node_string << ((mux_node_t *) (*it))->input_shifts[i];
                }
                else
                  pre_node_string << "NaN";
                if (i < (int) ((mux_node_t *) (*it))->input_shifts.size() - 1)
                  pre_node_string << ";";
              }
              pre_node_string << "]";
            }
          }
        }

        if (is_a<register_node_t>(*(*it)) || is_a<output_node_t>(*(*it)))
        {
          if (is_a<register_node_t>(*(*it)))
            type_string << "'R',";
          else
            type_string << "'O',";

          pre_node_string << "[";
          //set outputs of pre_node
          for (int i = 0, i_end = (int) ((register_node_t *) (*it))->input->output_factor.size(); i < i_end; ++i)
          {
            for (int j = 0, j_end = (int) ((register_node_t *) (*it))->input->output_factor[i].size(); j < j_end; ++j)
            {
              if (j == (int) ((register_node_t *) (*it))->input->output_factor[i].size() - 1)
              {
                if (((register_node_t *) (*it))->input->output_factor[i][j] != DONT_CARE)
                {
                  pre_node_string << ((register_node_t *) (*it))->input->output_factor[i][j];
                }
                else
                  pre_node_string << "NaN";
              }
              else
              {
                if (((register_node_t *) (*it))->input->output_factor[i][j] != DONT_CARE)
                  pre_node_string << ((register_node_t *) (*it))->input->output_factor[i][j] << ",";
                else
                  pre_node_string << "NaN,";
              }
            }
            if (i < (int) (*it)->output_factor.size() - 1)
            {
              pre_node_string << ";";
            }
          }
          pre_node_string << "]," << ((register_node_t *) (*it))->input->stage;
          if (((register_node_t *) (*it))->input_shift != 0)
            pre_node_string << "," << ((register_node_t *) (*it))->input_shift;
        }

        pre_node_string << "}";


        if (!(*it)->specific_parameters.empty())
        {
          pre_node_string << ":{";
          for (std::map<std::string, std::string>::iterator its = (*it)->specific_parameters.begin();
               its != (*it)->specific_parameters.end(); ++its)
          {
            if (!(*its).first.empty() && !(*its).second.empty())
            {
              if (its != (*it)->specific_parameters.begin())
                pre_node_string << ",";

              pre_node_string << (*its).first << "=" << (*its).second;
            }
          }
          pre_node_string << "}";
        }

        graphoutstream << type_string.str() << node_string.str() << pre_node_string.str();
        //fprintf(graphfilepointer,"%s%s%s",type_string.str().c_str(),node_string.str().c_str(),pre_node_string.str().c_str());

      }
    }
    graphoutstream << "}";
    if (!specific_parameters.empty())
    {
      graphoutstream << ":{";
      for (std::map<std::string, std::string>::iterator its = specific_parameters.begin();
           its != specific_parameters.end(); ++its)
      {
        if (!(*its).first.empty() && !(*its).second.empty())
        {
          if (its != specific_parameters.begin())
            graphoutstream << ",";

          graphoutstream << (*its).first << "=" << (*its).second;
        }
      }
      graphoutstream << "}";
    }

    //fprintf(graphfilepointer,"}\n");

    graphoutstream.flush();

    //fflush(graphfilepointer);
    //fclose(graphfilepointer);
    if (!quiet) std::cout << "Finished Syntax Output" << std::endl;
  }

  /**
   * replace_node replaces original by replacement by correcting all dependent inputs
   * Note: The function does not manipulate the node_list !
   *
   * */
  void adder_graph_t::replace_node(adder_graph_base_node_t *original, adder_graph_base_node_t *replacement)
  {
    if (!quiet) cout << "replacing " << original << " by new replacement " << replacement << endl;

    for (adder_graph_base_node_t *n: nodes_list) //go through all nodes and search for occurrences of the original node
    {
      //each node type has to be treated a bit different
      if (is_a<adder_subtractor_node_t>(*n) || is_a<conf_adder_subtractor_node_t>(*n))
      {
        for (int i = 0; i < ((adder_subtractor_node_t *) n)->inputs.size(); i++)
        {
          if (((adder_subtractor_node_t *) n)->inputs[i] == original)
          {
            ((adder_subtractor_node_t *) n)->inputs[i] = replacement; //replace replacement
            if (!quiet) cout << "replacing input " << i << " of replacement " << n << ": replacement " << original << " is replaced by " << replacement << endl;

            if(original->output_factor != replacement->output_factor)
            {
              cout << "signs of output factors are different. Correcting sign from " << ((adder_subtractor_node_t *) n)->input_is_negative[i] << " to " << !(((adder_subtractor_node_t *) n)->input_is_negative[i]) << endl;
              ((adder_subtractor_node_t *) n)->input_is_negative[i] = !(((adder_subtractor_node_t *) n)->input_is_negative[i]);

            }


          }
        }
      }
      else if (is_a<mux_node_t>(*n))
      {
        for (int i = 0; i < ((mux_node_t *) n)->inputs.size(); i++)
        {
          if (((mux_node_t *) n)->inputs[i] == original)
          {
            ((mux_node_t *) n)->inputs[i] = replacement; //replace replacement
            if (!quiet) cout << "replacing input " << i << " of replacement " << n << ": replacement " << original << " is replaced by " << replacement << endl;
          }
        }
      }
      else if (is_a<mux_node_t>(*n))
      {
        for (int i = 0; i < ((mux_node_t *) n)->inputs.size(); i++)
        {
          if (((mux_node_t *) n)->inputs[i] == original)
          {
            ((mux_node_t *) n)->inputs[i] = replacement; //replace replacement
            if (!quiet) cout << "replacing input " << i << " of replacement " << n << ": replacement " << original << " is replaced by " << replacement << endl;
          }
        }
      }
      else if ((is_a<register_node_t>(*n)) || (is_a<output_node_t>(*n)))
      {
        if (((register_node_t *) n)->input == original)
        {
          if (!quiet) cout << "replacing input of replacement " << n << ": replacement " << original << " is replaced by " << replacement << endl;
          ((register_node_t *) n)->input = replacement;
        }
      }
      else if ((is_a<adder_graph_base_node_t>(*n)) || (is_a<input_node_t>(*n)))
      {
        //ignore temporary nodes, those don't have an assigned input
      }
      else
      {
        throw runtime_error("unsupported replacement type " + string(typeid(*n).name()));
      }
    }
    delete original;
  }

  bool adder_graph_t::parse_to_graph(string commandLine, bool ignore_outnodes)
  {
    if (!quiet) cout << "parsing string: " << commandLine << endl;

    typedef enum
    {
      START, PARSING_NODE_BEGIN, PARSING_NODE_END, PARSING_NODE, PARSING_NODE_DELIM
    } state_t;
    state_t state = START;

//    noOfInputs = -1;
//    noOfConfigurations = -1;
    input_nodes_have_been_added = false;

    int pos = 0;
    bool complete = false;

    adder_graph_base_node_t *node;
    adder_graph_base_node_t *node_tmp;
    try
    {
      do
      {
        int nodeEndPos;
        string nodeStr;

        switch (state)
        {
          case START:
            if (commandLine.at(pos) == '{')
            {
              state = PARSING_NODE_BEGIN;
            }
            else
              throw runtime_error("{ expected");

            break;
          case PARSING_NODE_BEGIN:
            if (commandLine.at(pos) == '{')
            {
              state = PARSING_NODE;
            }
            else
              throw runtime_error("'{' expected (got '" + commandLine.substr(pos, 1) + "')");

            break;
          case PARSING_NODE:
            nodeEndPos = commandLine.find('}', pos);
            nodeStr = commandLine.substr(pos, nodeEndPos - pos);
            node = parse_node(nodeStr);


            if(factor_is_zero(node->output_factor))
            {
              node_tmp = nullptr; //ignore output nodes as there will be a node with same output factor
            }
            else
            if (!is_a<output_node_t>(*node))
            {
              node_tmp = get_node_from_output_factor_in_stage(node->output_factor, node->stage);
            }
            else
            {
              node_tmp = nullptr; //ignore output nodes as there will be a node with same output factor
            }

            if (node_tmp != nullptr)
            {
              if (is_a<adder_graph_base_node_t>(*node_tmp))
              {
                //temporary node found, remove that node and replace by new node:
                if(node_tmp->output_factor == node->output_factor)
                {
                  if (!quiet) cout << "found identical output factor, replace node " << node_tmp << " by new node " << node << endl;
                }
                else
                {
                  if (!quiet) cout << "found negated output factor" << endl;
                }
                if (!quiet) cout << "replace node " << node_tmp << " by new node " << node << endl;
                replace_node(node_tmp, node);
                if (!quiet) cout << "removing node " << node_tmp << " to list" << endl;
                nodes_list.remove(node_tmp);
                if (!quiet) cout << "inserting node " << node << " to list" << endl;
                nodes_list.push_back(node);
              }
              else if (is_a<output_node_t>(*node_tmp))
              {
                //ignore output nodes as those will have the same value as other nodes and are not temporary nodes
              }
              else
              {
                if(node_tmp->output_factor != node->output_factor)
                {
                  if(!quiet) cout << "node with different sign was already there, inserting new node " << node << " to list" << endl;
                  nodes_list.push_back(node); // add node to node list
                }
                else
                {
                  throw runtime_error("unsupported node type of temporary node " + string(typeid(*node_tmp).name()));
                }
              }
            }
            else
            {
              if (!quiet) cout << "inserting node " << node << " to list" << endl;
              nodes_list.push_back(node); // add node to node list
            }

            pos += nodeEndPos - pos - 1;
            state = PARSING_NODE_END;
            break;
          case PARSING_NODE_END:
            if (commandLine.at(pos) == '}')
            {
              state = PARSING_NODE_DELIM;
            }
            else
              throw runtime_error("'}' expected (got '" + commandLine.substr(pos, 1) + "')");

            break;
          case PARSING_NODE_DELIM:
            if (commandLine.at(pos) == ',')
            {
              state = PARSING_NODE_BEGIN;
            }
            else if (commandLine.at(pos) == '}')
            {
              complete = true;
            }
            else
              throw runtime_error("',' or '}' expected (got '" + commandLine.substr(pos, 1) + "')");
            break;
          default:
            throw runtime_error("invalid state in parse_to_graph()");
        }
        pos++;
      } while (!complete);
    }
    catch (const std::exception &e)
    {
      // cout << commandLine << endl;
      // for (int i = 0; i < pos; i++)
      //   cout << "-";
      // cout << "^" << endl;
      // cout << "parsing error: " << e.what() << endl;
      // exit(-1);//!!!
    }

    if (!quiet) cout << "parsing successful!" << endl;
    return true;
  }


  adder_graph_base_node_t *adder_graph_t::parse_node(string nodeStr)
  {
    if (!quiet) cout << "=== parsing node {" << nodeStr << "} ===" << endl;

    typedef enum
    {
      NODE_ID_DELIMITER1,      // first "'"
      NODE_ID,                 // "A", "R", "O" or "M"
      NODE_ID_DELIMITER2,      // second "'"
      NODE_ELEMENT_DELIMITER,  // the "," separating the contents
      NODE_OUT_FACTOR,         // the output factor, format is e.g. "[1,2;3,4]"
      NODE_OUT_STAGE,          // output stage (int value)
      NODE_OUT_RIGHT_SHIFT,    //加法后位移
      NODE_ARG_VALUE,          // an argument value of a node, format is e.g. "[1,2;3,4]"
      NODE_ARG_STAGE,          // input stage (int value)
      NODE_ARG_SHIFT           // input shift (int value)
    } state_t;
    state_t state = NODE_ID_DELIMITER1;
    state_t stateNext;

    int pos = 0;
    bool complete = false;
    char nodeId;
    int elemEndPos, elemBeginPos;
    string elemStr;
    int stage;
    int shift;
    int argNo = 0;
    adder_graph_base_node_t *node;
    bool requiresReconfAddSub = false; //flag indicating that this node is an adder is reconfigured at runtime

    try
    {
      std::vector<std::vector<int64_t> > factor_norm;
      std::vector<std::vector<int64_t> > factor;

      adder_graph_base_node_t *input_node = nullptr;

      do
      {

        switch (state)
        {
          case NODE_ID_DELIMITER1:
            if (nodeStr.at(pos) == '\'')
            {
              state = NODE_ID;
            }
            else
              throw runtime_error("'\'' expected (got '" + nodeStr.substr(pos, 1) + "')");

            break;
          case NODE_ID:
            nodeId = nodeStr.at(pos);
            if (!quiet) cout << "State NODE_ID: found node of type " << nodeId << endl;
            switch (nodeId)
            {
              case 'A':
                node = new adder_subtractor_node_t();
                break;
              case 'R':
                node = new register_node_t();
                break;
              case 'O':
                node = new output_node_t();
                break;
              case 'M':
                node = new mux_node_t();
                break;
              default:
                throw runtime_error("node with identifier '" + nodeStr.substr(pos, 1) + "' unknown or not supported");
            }
            if (!quiet) cout << "creating node " << node << endl;
            state = NODE_ID_DELIMITER2;

            break;
          case NODE_ID_DELIMITER2:
            if (nodeStr.at(pos) == '\'')
            {
              state = NODE_ELEMENT_DELIMITER;
              stateNext = NODE_OUT_FACTOR;
            }
            else
              throw runtime_error("'\'' expected (got '" + nodeStr.substr(pos, 1) + "')");

            break;
          case NODE_ELEMENT_DELIMITER:
            if (!quiet) cout << "State NODE_ELEMENT_DELIMITER" << endl;
            if (nodeStr.at(pos) == ',')
            {
              state = stateNext;
            }
            else
              throw runtime_error("',' expected (got '" + nodeStr.substr(pos, 1) + "')");

            break;

          case NODE_OUT_FACTOR:
            if (!quiet) cout << "State NODE_OUT_FACTOR" << endl;
            elemBeginPos = nodeStr.find('[', pos);

            if (elemBeginPos == string::npos)
              throw runtime_error("'[' expected (got '" + nodeStr.substr(pos, 1) + "')");

            elemEndPos = nodeStr.find(']', pos + 1);
            if (elemEndPos == string::npos)
              throw runtime_error("']' expected");

            elemStr = nodeStr.substr(elemBeginPos, elemEndPos - elemBeginPos + 1);
            parse_factor(elemStr, &(node->output_factor));
            pos += elemEndPos - elemBeginPos;
            state = NODE_ELEMENT_DELIMITER;
            stateNext = NODE_OUT_STAGE;
            break;

          case NODE_OUT_STAGE:
            if (!quiet) cout << "State NODE_OUT_STAGE" << endl;
            elemEndPos = nodeStr.find(',', pos);
            elemStr = nodeStr.substr(pos, elemEndPos - pos);

            node->stage = stoi(elemStr);

            if (!quiet) cout << "State NODE_OUT_STAGE: stage of node is " << node->stage << endl;

            pos += elemEndPos - pos - 1;
            state = NODE_ELEMENT_DELIMITER;
            stateNext = NODE_OUT_RIGHT_SHIFT;
            argNo = 1;
            break;

          //加法器右移
          case NODE_OUT_RIGHT_SHIFT:
            if (!quiet) cout << "State NODE_OUT_RIGHT_SHIFT" << endl;
            elemEndPos = nodeStr.find(',', pos);
            elemStr = nodeStr.substr(pos, elemEndPos - pos);

            node->output_right_shift = stoi(elemStr);

            if (!quiet) cout << "State NODE_OUT_RIGHT_SHIFT: stage of node is " << node->output_right_shift << endl;

            pos += elemEndPos - pos - 1;
            state = NODE_ELEMENT_DELIMITER;
            stateNext = NODE_ARG_VALUE;
            argNo = 1;
            break;

          case NODE_ARG_VALUE:
            if (!quiet) cout << "State NODE_ARG_VALUE" << endl;
            elemBeginPos = nodeStr.find('[', pos);

            if (elemBeginPos == string::npos)
              throw runtime_error("'[' expected (got '" + nodeStr.substr(pos, 1) + "')");

            elemEndPos = nodeStr.find(']', pos + 1);
            if (elemEndPos == string::npos)
              throw runtime_error("']' expected");

            elemStr = nodeStr.substr(elemBeginPos, elemEndPos - elemBeginPos + 1);

            if (!quiet) cout << "State NODE_ARG_VALUE: arg value is " << elemStr << endl;
            parse_factor(elemStr, &factor);

            factor_norm = normalize(factor);

            if (!quiet) cout << "normalized factor is " << factor_norm << endl;

            //If node is an adder node, check if signs of elements are different, if true replace adder_subtractor_node by conf_adder_subtractor_node
            if (nodeId == 'A')
            {
              if(!requiresReconfAddSub) //only check when it was not checked to be true in another argument before
              {
                enum sign_type {UNDECIDED, NEGATIVE, POSITIVE} sign=UNDECIDED;
                for (int c = 0; c < factor.size(); c++) //loop configurations
                {
                  for (int i = 0; i < factor[c].size(); i++) //loop inputs
                  {
                    if(factor[c][i] != 0)
                    {
                      if(factor[c][i] == factor_norm[c][i])
                      {
                        if(sign == UNDECIDED)
                        {
                          sign = POSITIVE; //first non-zero element found showing that sign is positive
                        }
                        else if(sign == NEGATIVE) //not all elements are positive, this is only possible with an add/sub
                        {
                          requiresReconfAddSub = true;
                          break;
                        }
                      }
                      else if(factor[c][i] == -factor_norm[c][i])
                      {
                        if(sign == UNDECIDED)
                        {
                          sign = NEGATIVE; //first non-zero element found showing that sign is positive
                        }
                        else if(sign == POSITIVE) //not all elements are negative, this is only possible with an add/sub
                        {
                          requiresReconfAddSub = true;
                          break;
                        }
                      }
                      else throw runtime_error("Abs value does not be equal between normalized and non-normalized vector (should never appear)");
                    }
                  }
                }

                //the following has to be done only once per node:
                if(requiresReconfAddSub)
                {
//                  factor_norm = abs_vec(factor); //normalization requires the absolute of all values as signs can change in each configuration

                  if (!quiet) cout << "Node requires a reconfigurable adder/subtractor, replacing node" << endl;

                  //replace adder_subtractor_node by conf_adder_subtractor_node:
                  adder_graph_base_node_t* node_orig = node;
                  node = new conf_adder_subtractor_node_t();
                  node->output_factor = node_orig->output_factor;
                  node->stage = node_orig->stage;

                  //copy the inputs assigned so far:
                  if(((adder_subtractor_node_t*) node_orig)->inputs.size() > 0 && !quiet) cout << "copy the inputs assigned so far" << endl;
                  ((conf_adder_subtractor_node_t*) node)->inputs.resize(((adder_subtractor_node_t*) node_orig)->inputs.size());

                  for(int i=0; i < ((adder_subtractor_node_t*) node_orig)->inputs.size(); i++)
                  {
                    //copy input:
                    ((adder_subtractor_node_t*) node)->inputs[i] = ((adder_subtractor_node_t*) node_orig)->inputs[i];

                    //copy signs:
                    ((conf_adder_subtractor_node_t*) node)->input_is_negative.resize(((adder_subtractor_node_t*) node_orig)->input_is_negative.size());

                    for(int j=0; j < ((adder_subtractor_node_t*) node_orig)->input_is_negative.size(); j++) //loop over the different configurations
                    {
                      ((conf_adder_subtractor_node_t*) node)->input_is_negative[j].resize(((adder_subtractor_node_t*) node_orig)->input_is_negative.size());
                      for(int k=0; k < ((conf_adder_subtractor_node_t*) node)->input_is_negative[j].size(); k++) //loop over the inputs
                      {
                        ((conf_adder_subtractor_node_t*) node)->input_is_negative[j][k] = ((adder_subtractor_node_t*) node_orig)->input_is_negative[i]; // so far, there was no reconfiguration, hence, it is save to copy all signs
                      }
                    }

                  }

                  //copy the shifts assigned so far:
                  ((conf_adder_subtractor_node_t*) node)->input_shifts.resize(((adder_subtractor_node_t*) node_orig)->input_shifts.size());
                  for(int i=0; i < ((adder_subtractor_node_t*) node_orig)->input_shifts.size(); i++)
                  {
                    ((conf_adder_subtractor_node_t*) node)->input_shifts[i] = ((adder_subtractor_node_t*) node_orig)->input_shifts[i];
                  }
                  replace_node(node_orig, node);
                }
              }
            }
            pos += elemEndPos - elemBeginPos;
            state = NODE_ELEMENT_DELIMITER;
            stateNext = NODE_ARG_STAGE;
            break;

          case NODE_ARG_STAGE:
            if (!quiet) cout << "State NODE_ARG_STAGE" << endl;
            elemEndPos = nodeStr.find(',', pos);
            if (elemEndPos > 0)
            {
              elemStr = nodeStr.substr(pos, elemEndPos - pos); //a comma was found
            }
            else
            {
              elemStr = nodeStr.substr(pos); //no comma was found, it is the last argument (e.g., for a register)...
              complete = true; //... then we are done
            }

            stage = stoi(elemStr);

            if (!quiet) cout << "State NODE_ARG_STAGE: stage of arg is " << stage << endl;

            input_node = get_node_from_output_factor_in_stage(factor, stage);

            if (input_node == nullptr)
            {
              if(factor_is_zero(factor_norm))
              {
                //create zero node:
                input_node = new zero_node_t();
                input_node->output_factor = factor_norm;
                input_node->stage = stage;
                if (!quiet) cout << "adding zero node " << input_node << " with factor " << input_node->output_factor << " at stage " << input_node->stage << endl;
              }
              else
              {
                //create temporary node:
                input_node = new adder_graph_base_node_t();
                input_node->output_factor = factor;
                input_node->stage = stage;
                if (!quiet) cout << "adding temporary node " << input_node << " with factor " << input_node->output_factor << " at stage " << input_node->stage << endl;
              }
              nodes_list.push_back(input_node);
            }
            if (nodeId == 'A')
            {
              int noOfConfigurations = node->output_factor.size();
              if (!quiet) cout << "adding input to node " << node << endl;
              ((adder_subtractor_node_t *) node)->inputs.push_back(input_node);
              if (is_a<adder_subtractor_node_t>(*node))
              {
                if (!quiet) cout << "setting input sign to " << (input_node->output_factor == factor ? "+" : "-") << endl;
                ((adder_subtractor_node_t *) node)->input_is_negative.push_back(input_node->output_factor != factor);

              }
              else if (is_a<conf_adder_subtractor_node_t>(*node))
              {
                if(!quiet) cout << "setting input_is_negative of conf_adder_subtractor_node_t" << endl;
                ((conf_adder_subtractor_node_t *) node)->input_is_negative.resize(noOfConfigurations);
                if(!quiet) cout << "input_is_negative = [";
                for (int c = 0; c < noOfConfigurations; c++) //loop configurations
                {
                  if(((conf_adder_subtractor_node_t *) node)->input_is_negative[c].size() < argNo)
                    ((conf_adder_subtractor_node_t *) node)->input_is_negative[c].resize(argNo);

                  ((conf_adder_subtractor_node_t *) node)->input_is_negative[c][argNo-1] = (factor[c] != factor_norm[c]); //check the whole vector (for all inputs)
                  if(!quiet)
                  {
                    cout << ((conf_adder_subtractor_node_t *) node)->input_is_negative[c][argNo-1];
                    if(c < ((conf_adder_subtractor_node_t *) node)->inputs.size()-1) cout << ",";
                  }
                }
                if(!quiet) cout << "]" << endl;
              }
              else throw runtime_error("Adder node not of type adder_subtractor_node_t or conf_adder_subtractor_node_t!");
            }
            else if (nodeId == 'M')
            {
              //the inputs are added later when the input shifts are known
            }
            else if ((nodeId == 'R') || (nodeId == 'O'))
            {
              if (!quiet) cout << "adding input to node " << node << endl;
              ((register_node_t *) node)->input = input_node;
            }

            pos += elemEndPos - pos - 1;
            state = NODE_ELEMENT_DELIMITER;
            stateNext = NODE_ARG_SHIFT;

            /*
            if (nodeId == 'R')
              complete = true; //for registers, we are done
            */
            break;

          case NODE_ARG_SHIFT:
            if (!quiet) cout << "State NODE_ARG_SHIFT" << endl;
            elemEndPos = nodeStr.find(',', pos); //first, check if more arguments are present. If this is the case, parse until the next ','
            if (elemEndPos < 0)
            {
              //no ',' was found, so take the shift until the end of the string, and ...
              elemStr = nodeStr.substr(pos);
              complete = true; //... then we are done!
            }
            else
            {
              elemStr = nodeStr.substr(pos, elemEndPos - pos);
            }

            if (nodeId == 'M')
            {
              std::vector<std::vector<int64_t> > f;
              parse_factor(elemStr, &f);
              assert(f.size() > 0);

              int noOfConfigurations = factor.size();
              ((mux_node_t *) node)->inputs.resize(noOfConfigurations);
              ((mux_node_t *) node)->input_shifts.resize(noOfConfigurations);

              for (int c = 0; c < f.size(); c++)
              {
                if (f[c][0] != DONT_CARE) //only add input if it is not a don't care
                {
                  if (!quiet) cout << "setting input of node " << node << " to " << input_node << " with shift " << f[c][0] << " for configuration " << c << endl;
                  ((mux_node_t *) node)->inputs[c] = input_node;
                  ((mux_node_t *) node)->input_shifts[c] = f[c][0];
                }
              }
            }
            else if (nodeId == 'A')
            {
              shift = stoi(elemStr);
              ((adder_subtractor_node_t *) node)->input_shifts.push_back(shift);
            }
            else if ((nodeId == 'R') || (nodeId == 'O'))
            {
              shift = stoi(elemStr);
              ((register_node_t *) node)->input_shift = shift;
              complete = true; //for registers, we are done
            }
            else
            {
              throw runtime_error("node type not handled");
            }

            pos += elemEndPos - pos - 1;
            state = NODE_ELEMENT_DELIMITER;
            stateNext = NODE_ARG_VALUE;
            argNo++;


            if (pos + 1 == nodeStr.length() - 1)
              complete = true; //per default, we are done when reaching the end of the string
            break;

          default:
            cout << "invalid state in parse_node()" << endl;
            throw runtime_error("invalid state");
        }
        pos++;
      } while (!complete);
    }
    catch (const std::exception &e)
    {
      cout << "node parsing error: " << e.what() << endl;
      cout << nodeStr << endl;
      for (int i = 0; i < pos + 1; i++)
        cout << "-";
      cout << "^" << endl;
      exit(-1);
    }
    return node;
  }

  bool adder_graph_t::parse_factor(string factorStr, std::vector<std::vector<int64_t> > *factor)
  {
    int noOfConfigurations=0;
    int noOfInputs=0;

    if (!quiet) cout << "parsing factor " << factorStr << endl;

    factorStr = factorStr.substr(1, factorStr.length() - 2); //remove closing brackets

    size_t p=0;
    do
    {
      p = factorStr.find_first_of(';',p);
      noOfConfigurations++;
      if(p != string::npos) p++;
    } while(p != string::npos);
    //we counted one semicolon too much and we have one noOfConfigurations = noOfSemicolons + 1 -> nothing further to do!

    int noOfCommas=0;
    p=0;
    do
    {
      p = factorStr.find_first_of(',',p);
      noOfCommas++;
      if(p != string::npos) p++;
    } while(p != string::npos);
    noOfCommas--; //we counted one comma too much -> decrement

    noOfInputs = noOfCommas/noOfConfigurations+1;
    if (!quiet) cout << "  factor has " << noOfConfigurations << " configuration(s)" << " and " << noOfInputs << " input(s)" << endl;
    assert(noOfCommas % noOfConfigurations == 0);

    if(!input_nodes_have_been_added) //the first factor we find determines the noOfInputs -> add input nodes, such that they can be found to avoid temp nodes
    {
      add_input_nodes(noOfConfigurations,noOfInputs);
      input_nodes_have_been_added=true;
    }


    typedef enum
    {
      FACTOR_DELIMITER, FACTOR_VALUE
    } state_t;
    state_t state = FACTOR_VALUE;

    int row = 0; //row of factor matrix
    int col = 0; //column of factor matrix

    int pos = 0;
    bool complete = false;
    int commaPos;
    int semicolonPos;
    int stringEndPos;
    int minPos;
    int64_t value;

    string elemStr;

    if ((noOfInputs == -1) && (noOfConfigurations == -1))
    {
      factor->resize(1); //allocate the first row
    }
    else
    {
      factor->resize(noOfConfigurations, std::vector<int64_t>(noOfInputs)); //allocate complete matrix
    }

    try
    {
      do
      {
        switch (state)
        {
          case FACTOR_VALUE:
            commaPos = factorStr.find(',', pos);
            semicolonPos = factorStr.find(';', pos);
            stringEndPos = factorStr.length();

            commaPos < 0 ? commaPos = INT_MAX : commaPos;
            semicolonPos < 0 ? semicolonPos = INT_MAX : semicolonPos;

            minPos = min(commaPos, semicolonPos);
            minPos = min(minPos, stringEndPos);

            elemStr = factorStr.substr(pos, minPos - pos);

            flush(cout);
            if (elemStr == "NaN")
            {
              value = DONT_CARE;
            }
            else
            {
              value = stoll(elemStr);
            }

            //resize matrix if necessary:
            if ((*factor).size() <= row)
              (*factor).resize(row + 1);

            if ((*factor)[row].size() <= col)
              (*factor)[row].resize(col + 1);

            (*factor)[row][col] = value;

            pos = minPos - 1;

            state = FACTOR_DELIMITER;
            break;
          case FACTOR_DELIMITER:

            if (pos == factorStr.length())
            {
              complete = true;
              break;
            }
            else if (factorStr[pos] == ' ') //ignore spaces
              pos++;
            else if (factorStr[pos] == ',')
            {
              col++;
              state = FACTOR_VALUE;
            }
            else if (factorStr[pos] == ';')
            {
              row++;
              col=0;
              state = FACTOR_VALUE;
            }
            break;

          default:
            cout << "invalid state in parse_factor()" << endl;
            throw runtime_error("invalid state");
        }
        pos++;
      } while (!complete);
    }
    catch (const std::exception &e)
    {
      cout << "factor parsing error: " << e.what() << endl;
      cout << factorStr << endl;
      for (int i = 0; i < pos + 1; i++)
        cout << "-";
      cout << "^" << endl;
      exit(-1);
    }

/*
    if ((noOfInputs == -1) && (noOfConfigurations == -1))
    {
      noOfConfigurations = (*factor).size();
      noOfInputs = (*factor)[0].size();
      if (!quiet) cout << "We have " << noOfConfigurations << " configuration(s)" << endl;
      if (!quiet) cout << "We have " << noOfInputs << " input(s)" << endl;

      for (int i = 0; i < noOfInputs; i++)
      {
        input_node_t *input_node = new input_node_t();
        input_node->output_factor.resize(noOfConfigurations);
        for (int r = 0; r < input_node->output_factor.size(); r++)
        {
          input_node->output_factor[r].resize(noOfInputs);
          for (int c = 0; c < input_node->output_factor[r].size(); c++)
          {
            input_node->output_factor[r][c] = (i == c);
          }
        }
        nodes_list.push_back(input_node);
      }
    }
    else
    {
      assert(noOfConfigurations == (*factor).size());
      assert(noOfInputs == (*factor)[0].size());
    }
*/
    assert(noOfConfigurations == (*factor).size());
    assert(noOfInputs == (*factor)[0].size());
    if (!quiet) cout << "factor is a " << (*factor).size() << " x " << (*factor)[0].size() << " matrix: " << (*factor) << endl;
    return complete;

  }

  void adder_graph_t::add_input_nodes(int noOfConfigurations, int noOfInputs)
  {
    for (int i = 0; i < noOfInputs; i++)
    {
      input_node_t *input_node = new input_node_t();
      input_node->output_factor.resize(noOfConfigurations);
      for (int r = 0; r < input_node->output_factor.size(); r++)
      {
        input_node->output_factor[r].resize(noOfInputs);
        for (int c = 0; c < input_node->output_factor[r].size(); c++)
        {
          input_node->output_factor[r][c] = (i == c);
        }
      }
      if(!quiet) cout << "adding input node " << input_node->output_factor << endl;
      nodes_list.push_back(input_node);
    }
  }


  void adder_graph_t::get_cost()
  {
  }


  void adder_graph_t::print_cost()
  {
  }

void save_matrix_length_to_txt(const std::vector<std::vector<std::vector<float>>>& matrix_adj, const std::string& filename) {
    std::ofstream outfile(filename, std::ios::app);
    if (!outfile.is_open()) {
        std::cerr << "无法打开文件: " << filename << std::endl;
        return;
    }

    // 保存元素到文本文件
    outfile << matrix_adj[0].size();  // 用空格分隔每个元素
    outfile << "\n";  // 每一行后换行
    }

void save_matrix_to_txt(const std::vector<std::vector<std::vector<float>>>& matrix_adj, const std::string& filename) {
    std::ofstream outfile(filename, std::ios::app);
    if (!outfile.is_open()) {
        std::cerr << "无法打开文件: " << filename << std::endl;
        return;
    }

    // 遍历三维矩阵，逐层保存元素到文本文件
    for (const auto& matrix : matrix_adj) {
        for (const auto& row : matrix) {
            for (const auto& element : row) {
                outfile << element << " ";  // 用空格分隔每个元素
            }
            outfile << "\n";  // 每一行后换行
        }
        outfile << "\n";  // 每个子矩阵后换行
    }

    outfile.close();
    num_matrix += 1;
    //std::cout << "矩阵已保存到文本文件: " << filename << std::endl;
}

//*********************判断是否存在相同的矩阵**********************/
bool areVectorsEqual(const std::vector<float>& a, const std::vector<float>& b) {
    return a == b;
}

void compareMatrices(const std::vector<std::vector<float>>& matrix, const std::vector<std::vector<float>>& all_row,
                      const std::vector<std::vector<std::vector<float>>>& matrix_adj,const std::vector<float>& factors,const int& count) {
    int found = 0;
    int num_row = matrix.size();
    for (const auto& row1 : matrix) {
        for (const auto& row2 : all_row) {
            if (areVectorsEqual(row1, row2)) {
                found += 1;
                break;
            }
        }
    }
    if (found == num_row) {
        //std::cout << "相同矩阵" << std::endl;
    } else {
        //添加矩阵
        save_matrix_to_txt(matrix_adj, "matrix_output.txt");
        //   std::cout << "Add Matrix" << std::endl;
        // //打印矩阵内容（每个元素是空的向量）
        // std::cout << "	" ;
        // for (const auto& f : factors) {
        //   std::cout << f << "	";
        //   }
        // std::cout << endl;
        // for (int i = 0; i < count; ++i) {
        // std::cout << factors[i] << "	";
        //     for (int j = 0; j < count; ++j) {
        //     std::cout <<'[';
        //             for (const auto& s : matrix_adj[i][j]) {
        //       std::cout << s << " ";
        //       }
        //     std::cout << ']' <<'	';
        //         }
            
        //     std::cout << std::endl;
        // }
          }
}

std::vector<std::vector<float>> flattenMatricesTo1D(const std::vector<std::vector<std::vector<float>>>& matrix_adj) {
    std::vector<std::vector<float>> flattened_matrices;

    for (const auto& matrix : matrix_adj) {
        std::vector<float> flattened_row;
        for (const auto& row : matrix) {
            flattened_row.insert(flattened_row.end(), row.begin(), row.end());
        }
        flattened_matrices.push_back(flattened_row);
    }

    return flattened_matrices;
}

std::vector<std::vector<float>> read_matrix_from_txt(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "无法打开文件: " << filename << std::endl;
        return {};
    }
    std::vector<std::vector<float>> all_row;
    std::vector<float> current_row;
    std::string line;
    
    while (std::getline(infile, line)) {
        line = line.substr(0, line.find_last_not_of(" \n\r\t") + 1);  // 去除行尾空白
        
        if (line.empty()) {
            // 处理空行
                all_row.push_back(current_row);
                current_row.clear();
        } else {
            // 处理非空行
            std::istringstream iss(line);
            int value;
            while (iss >> value) {
                current_row.push_back(value);
            }
        }
    }

    // 处理文件结束后剩余的矩阵和行
    if (!current_row.empty()) {
        all_row.push_back(current_row);
    }

    infile.close();
    return all_row;
}
/****************************************************************************/
  void adder_graph_t::print_graph(const std::vector<float>& constant)
  {
    adder_graph_t mygraph = *this;


    if (!mygraph.specific_parameters.empty())
    {
      std::cout << "Attached parameters to this graph:" << std::endl;
      for (std::map<std::string, std::string>::iterator its = mygraph.specific_parameters.begin();
           its != mygraph.specific_parameters.end(); ++its)
      {
        if (!(*its).first.empty() && !(*its).second.empty())
        {
          std::cout << "\t" << (*its).first << " = " << (*its).second << ";" << std::endl;
        }
      }
    }
	//得到节点的个数
	int count = 0;
	for (std::list<adder_graph_base_node_t *>::iterator it = mygraph.nodes_list.begin(), it_end = mygraph.nodes_list.end(); it != it_end; ++it)
    {
			count++;
	}
	//std::cout << "The numbers of nodes are: "<<count << endl;

	// 定义邻接矩阵矩阵类型，矩阵中的每个元素是一个空的向量
  std::vector<std::vector<std::vector<float>>> matrix_adj(count, std::vector<std::vector<float>>(count,std::vector<float>(1,0)));
	std::vector<float> factors;
	std::vector<adder_graph_base_node_t *> nodes;
  std::unordered_set<int> seen;
  std::vector<int> uniqueVec;
  int found_same = 0;
    for (std::list<adder_graph_base_node_t *>::iterator it = mygraph.nodes_list.begin(), it_end = mygraph.nodes_list.end(); it != it_end; ++it)
    {
      adder_graph_base_node_t *p = *it;
      //std::cout << "\nNode type of node " << (*it) << " in stage " << (*it)->stage << " is: " << typeid(*p).name() << std::endl;
      //std::cout << "The factor  is" << (*it)->output_factor << endl;
      //std::cout << "The outputs are: ";
		int output;
	  
      //std::cout << "Size:"<< (*it)->output_factor.size();
      for (unsigned int k = 0, k_end = (*it)->output_factor.size(); k < k_end; ++k)
      {
        if (k > 0) std::cout << "; ";
        for (unsigned int f = 0, f_end = (*(*it)).output_factor[k].size(); f < f_end; ++f)
        {
          if (f > 0) std::cout << ", ";

          //std::cout << (*(*it)).output_factor[k][f] << " ";
		  factors.push_back((*(*it)).output_factor[k][f]);
		  output = (*(*it)).output_factor[k][f];
        }
		
      }
      //std::cout << std::endl;
      //std::cout << "加法器右移" << (*it)->output_right_shift << endl;

      if (is_a<input_node_t>(*(*it)))
      {
        //std::cout << "No input" << std::endl;
      }
      else if (is_a<adder_subtractor_node_t>(*(*it)))
      {
        adder_subtractor_node_t *t = (adder_subtractor_node_t *) (*it);

        for (int k = 0, k_end = (int) t->inputs.size(); k < k_end; ++k)
          if (t->inputs[k] == nullptr)
          {
            std::cout << "Input: unconnected" << std::endl;
          }
          else
          {
            int input=  t->inputs[k]->output_factor[0][0];
            // std::cout << "Input: " << t->inputs[k] << ", shift: " << t->input_shifts[k] << ", in stage "
            // << t->inputs[k]->stage << ", is negative: " << t->input_is_negative[k] << std::endl;
            // std::cout << "Input: " << t->inputs[k]->output_factor << ", shift: " << t->input_shifts[k] << ", in stage "
            // << t->inputs[k]->stage << ", is negative: " << t->input_is_negative[k] << std::endl;
            //找输入节点在邻接矩阵的索引
            auto it_in = std::find(factors.begin(), factors.end(), input);
            int index_in = std::distance(factors.begin(), it_in);
            //找输出节点在邻接矩阵的索引
            //auto it_out = std::find(factors.begin(), factors.end(), output);
            //int index_out = std::distance(factors.begin(), it_out);
            int index_out = factors.size() - 1;
            
            //输入是否取负
            if(t->input_is_negative[k] == 0){
              if (matrix_adj[index_in][index_out].size() == 1 && matrix_adj[index_in][index_out].back() == 0) {
                // 如果向量只有一个元素且该元素为0，则删除该元素
                matrix_adj[index_in][index_out].pop_back();
              }
              matrix_adj[index_in][index_out].push_back(pow(2,(t->input_shifts[k])));
            }
            else{
              if (matrix_adj[index_in][index_out].size() == 1 && matrix_adj[index_in][index_out][0] == 0) {
                // 如果向量只有一个元素且该元素为0，则删除该元素
                matrix_adj[index_in][index_out].clear();
              }
              matrix_adj[index_in][index_out].push_back(-pow(2,(t->input_shifts[k])));
            }
          }

          //factors去重
          // for (const int& num : factors) {
          //     // 如果 num 不在 seen 中，添加到 uniqueVec，并标记为 seen
          //     //if (seen.insert(num).second) {
          //         uniqueVec.push_back(num);
          //     //}
          // }
          //如果有加法后右移
          //找输入节点在邻接矩阵的索引
          //找输出节点在邻接矩阵的索引
          // auto it_out = std::find(factors.begin(), factors.end(), output);
          // int index_out = std::distance(factors.begin(), it_out);
          int index_out = factors.size() - 1;
          //如果有加法后右移，往对角线上添加小数
          if(t->output_right_shift != 0){
            if (matrix_adj[index_out][index_out].size() == 1 && matrix_adj[index_out][index_out][0] == 0) {
              // 如果向量只有一个元素且该元素为0，则删除该元素
                matrix_adj[index_out][index_out].clear();
                matrix_adj[index_out][index_out].push_back(pow(2,-(t->output_right_shift)));
            }
          }
      }
      else if (is_a<mux_node_t>(*(*it)))
      {
        mux_node_t *t = (mux_node_t *) (*it);
        for (int k = 0, k_end = (int) t->inputs.size(); k < k_end; ++k)
          if (t->inputs[k] != nullptr)
          {
            std::cout << "Input: " << t->inputs[k] << ", shift: ";
            if (t->input_shifts[k] != DONT_CARE)
              std::cout << t->input_shifts[k];
            else
              std::cout << "NaN";
            std::cout << ", in stage " << t->inputs[k]->stage << std::endl;
          }
      }
      else if (is_a<conf_adder_subtractor_node_t>(*(*it)))
      {
        conf_adder_subtractor_node_t *t = (conf_adder_subtractor_node_t *) (*it);
        for (int k = 0, k_end = (int) t->inputs.size(); k < k_end; ++k)
        {
          if (t->inputs[k] == nullptr)
          {
            std::cout << "Input: unconnected" << std::endl;
          }
          else
          {
            std::cout << "Input: " << t->inputs[k] << ", shift: " << t->input_shifts[k] << ", in stage "
                      << t->inputs[k]->stage << ", is negative: ";
            for (int j = 0, j_end = (int) t->input_is_negative.size(); j < j_end; ++j)
            {
              std::cout << t->input_is_negative[j][k];
              if(j < j_end-1) std::cout << ", ";
            }
            std::cout << std::endl;
          }
        }
      }
      else if ((is_a<register_node_t>(*(*it))) || (is_a<output_node_t>(*(*it))))
      {
        register_node_t *t = (register_node_t *) (*it);
        if (t->input == nullptr)
        {
          std::cout << "Input: unconnected" << std::endl;
        }
        else
        {
          std::cout << "Input: " << t->input << ", shift: " << t->input_shift << ", in stage " << t->input->stage << std::endl;
        }
      }
      else if (is_a<zero_node_t>(*(*it)))
      {
        //nothing to do
      }
      else
      {
        stringstream ss;
        ss << "Node type of node " << *it << " unknown";
        throw runtime_error(ss.str());
      }

      if (!(*it)->specific_parameters.empty())
      {
        std::cout << "Attached parameters:" << std::endl;
        for (std::map<std::string, std::string>::iterator its = (*it)->specific_parameters.begin();
             its != (*it)->specific_parameters.end(); ++its)
        {
          if (!(*its).first.empty() && !(*its).second.empty())
          {
            std::cout << "\t" << (*its).first << " = " << (*its).second << ";" << std::endl;
          }
        }
      }
    }
    //std::cout << "Printed all" << std::endl;

	//节点
	// std::cout << "节点：";
	// for (size_t i = 0; i < factors.size(); ++i) {
	// std::cout << factors[i] << " ";
  //   }
  //   std::cout << std::endl;

	// // 打印矩阵内容（每个元素是空的向量）
	// std::cout << "	" ;
	// for (const auto& f : factors) {
	// 	std::cout << f << "	";
	// 	}
	// std::cout << endl;

  //std::cout <<constant/factors[-1]<<endl;

  //右下角是否需要左移
    // if (matrix_adj[count-1][count-1].size() == 1 && matrix_adj[count-1][count-1][0] != 0) {
    //       // 如果右下角出书存在加法后右移
    //       float temp = matrix_adj[count-1][count-1][0];
    //       matrix_adj[count-1][count-1].clear();
		// 		  matrix_adj[count-1][count-1].push_back(temp*(constant/(factors.back())));
		// 		}
		// if (matrix_adj[count-1][count-1].size() == 1 && matrix_adj[count-1][count-1][0] == 0) {
    //       // 添加新元素
    //       matrix_adj[count-1][count-1].clear();
		// 		  matrix_adj[count-1][count-1].push_back(constant/(factors.back()));
		// 		}
  
  /************************MCM部分****************************/
  //将加法后右移转移到input shift上
  int n = matrix_adj.size(); // 获取矩阵的大小

  for (int i = 0; i < n-1; ++i) {
      // 获取对角线元素
      float diag_element = matrix_adj[i][i][0];
      // 检查对角线元素是否在 (0, 1) 之间
        if (diag_element > 0 && diag_element < 1) {
            // 遍历第 i 行的右侧元素
            for (int j = i + 1; j < n; ++j) {
              std::vector<float> element = matrix_adj[i][j];
              // 逐个元素相乘
              std::vector<float> vec;
              for (const auto& elem : element) {
                vec.push_back(elem*diag_element);
              }
              matrix_adj[i][j] = vec;
            }
        }
      }
  
  //判断哪些节点是输出节点，记录节点坐标index_out和系数坐标index_c
  //对c1 c2的取值无限定，无需假设c1 != 2^n*c2 n=0,1,2,....
  std::vector<std::pair<int, int>> result;
  for (size_t j = 0; j < constant.size(); ++j) {
  //for (size_t j = constant.size() - 1; j > 0; --j) {
    int found = 0;
    //for (size_t i = 0; i < factors.size(); ++i) {
      for (size_t i = factors.size() -1 ; i > 0; --i) {
      int f = factors[i];
      int c = constant[j];
      // 检查 f 是否等于 c
      if (f == c) {
          result.emplace_back(i, j);
          found = 1;
      }
      // 检查 f 是否能整除 c
      else if (c % f == 0) {
          int k = c / f;
          // 检查 k 是否是 2 的幂
          if (k > 0 && (k & (k - 1)) == 0) {
              result.emplace_back(i, j);
              found = 1;
          }
      }
      if(found == 1){//c找到对应的输出节点
        break;
      }
    }
  }
  //1.直接输出 2.输出左移
  //在输出节点所在列的对角线上修改数据为1或者C/factors
  float temp = 0;
  std::vector<int> index_of_diag_element;
  for (int i = 0; i < n-1; ++i) {
      // 获取对角线元素
      float diag_element = matrix_adj[i][i][0];
      // 检查对角线元素是否在 (0, 1) 之间
        if (diag_element != 0) {
            // 记录对角线非零元素的坐标
            index_of_diag_element.push_back(i);
        }
      }
  std::vector<std::vector<std::vector<float>>> matrix_adj_origin = matrix_adj;
  for (const auto& pair : result) {
    //std::cout << "Factor index: " << pair.first << ", Constant index: " << pair.second << std::endl;
    if (std::find(index_of_diag_element.begin(), index_of_diag_element.end(), pair.first) != index_of_diag_element.end()) {
        // 存在c1 != 2^n*c2 n=0,1,2,....
        temp = matrix_adj_origin[pair.first][pair.first][0];
        // matrix_adj[pair.first][pair.first].clear();
        matrix_adj[pair.first][pair.first].push_back(temp*(constant[pair.second]/(factors[pair.first])));
      }
    else{
      if (matrix_adj[pair.first][pair.first].size() == 1 && matrix_adj[pair.first][pair.first][0] == 0) {
        // 直接输出
        matrix_adj[pair.first][pair.first].clear();
        }
        matrix_adj[pair.first][pair.first].push_back(constant[pair.second]/(factors[pair.first]));
      }
    //如果输出节点的对角线上存在非零元素
    // if (std::find(index_of_diag_element.begin(), index_of_diag_element.end(), pair.first) != index_of_diag_element.end()) {
    //   temp = matrix_adj[pair.first][pair.first][0];
    // }
  }

  //如果加法后右移没有作为输出节点，则对角线元素归0
  for (int i = 1; i < n-1; ++i) {
    int flag = 0;
    //如果该节点是输出节点，则跳过
    for(const auto& pair : result){
      if(pair.first == i)
        flag = 1;
    }
    if(flag == 1)
      continue;

    // 获取对角线元素
    float diag_element = matrix_adj[i][i][0];
    // 检查对角线元素是否在 (0, 1) 之间
    if (diag_element > 0 && diag_element < 1) {
          matrix_adj[i][i][0] = 0;
        }
    }

    // for (int i = 0; i < count; ++i) {
		// std::cout << factors[i] << "	";
    //     for (int j = 0; j < count; ++j) {
		// 		std::cout <<'[';
    //             for (const auto& s : matrix_adj[i][j]) {
		// 			std::cout << s << " ";
		// 			}
		// 		std::cout << ']' <<'	';
    //         }
        
    //     std::cout << std::endl;
    // }

    //保存加法器个数
    save_matrix_length_to_txt(matrix_adj, "matrix_length.txt");
    // //读取矩阵
    // auto all_row = read_matrix_from_txt("matrix_output.txt");
    // // 执行转换
    // std::vector<std::vector<float>> flattened_matrices = flattenMatricesTo1D(matrix_adj);
    // // 比较矩阵
    // compareMatrices(flattened_matrices, all_row,matrix_adj,factors,count);
    save_matrix_to_txt(matrix_adj, "matrix_output.txt");
  
  }


  void adder_graph_t::check_and_correct(string graphstring)
  {
    if (!quiet) cout << "entering adder_graph_t::check_and_correct" << endl;

    if (graphstring == "") graphstring = "STRINGSYN NOT AVAILABLE";

    list<string> errorList;
    map<int, vector<string> > configurationCountMap;
    map<int, vector<string> > inputCountMap;
    configurationCountMap.clear();
    inputCountMap.clear();

    //check for errors
    if (!quiet) cout << "checking for errors" << endl;
    for (list<adder_graph_base_node_t *>::iterator it = nodes_list.begin(), it_end = nodes_list.end(); it != it_end; ++it)
    {
      stringstream node_name;
      if (is_a<input_node_t>(*(*it))) node_name << "INP@";
      else if (is_a<adder_subtractor_node_t>(*(*it))) node_name << "ADD@";
      else if (is_a<mux_node_t>(*(*it))) node_name << "MUX@";
      else if (is_a<conf_adder_subtractor_node_t>(*(*it))) node_name << "CONF_ADD@";
      else if (is_a<register_node_t>(*(*it))) node_name << "REG@";
      node_name << "[";
      for (unsigned int cfg = 0, cfg_end = (*it)->output_factor.size(); cfg < cfg_end; ++cfg)
      {
        if (cfg > 0) node_name << ";";
        for (unsigned int v = 0, v_end = (*it)->output_factor[cfg].size(); v < v_end; ++v)
        {
          if (v > 0) node_name << ",";
          if ((*it)->output_factor[cfg][v] == DONT_CARE) node_name << "-";
          else node_name << (*it)->output_factor[cfg][v];
        }
      }
      node_name << "],stage " << (*it)->stage;

      if (!quiet) cout << "checking presence of output factor" << endl;
      if ((*it)->output_factor.size() == 0)
      {
        errorList.push_back("OUTPUTFACTOR MISSING: " + node_name.str());
        continue;
      }
      if ((*it)->output_factor[0].size() == 0)
      {
        errorList.push_back("OUTPUTFACTOR INNER MISSING: " + node_name.str());
        continue;
      }

      map<int, vector<string> >::iterator pr1 = configurationCountMap.find((*it)->output_factor.size());
      if (pr1 != configurationCountMap.end())
      {
        (*pr1).second.push_back(node_name.str());
      }
      else
      {
        vector<string> newvec;
        newvec.push_back(node_name.str());
        configurationCountMap.insert(make_pair((*it)->output_factor.size(), newvec));
      }

      map<int, vector<string> >::iterator pr2 = inputCountMap.find((*it)->output_factor[0].size());
      if (pr2 != inputCountMap.end())
      {
        (*pr2).second.push_back(node_name.str());
      }
      else
      {
        vector<string> newvec;
        newvec.push_back(node_name.str());
        inputCountMap.insert(make_pair((*it)->output_factor[0].size(), newvec));
      }

      if (!quiet) cout << "checking nodes" << endl;
      if (is_a<input_node_t>(*(*it)))
      {

      }
      else if (is_a<adder_subtractor_node_t>(*(*it)))
      {
        if (!quiet) cout << "checking adder node" << node_name.str() << endl;
        adder_subtractor_node_t *t = (adder_subtractor_node_t *) (*it);
        //check for NULL Inputs
        int inputs_notnull = 0;
        bool not_found_error = false;
        int max_neg_shift = 0;
        for (unsigned int i = 0, i_end = t->inputs.size(); i < i_end; ++i)
        {
          if (t->inputs[i] != NULL)
          {
            inputs_notnull++;
            if (std::find(nodes_list.begin(), nodes_list.end(), t->inputs[i]) == nodes_list.end())
            {
              errorList.push_back("NODE_LIST DOES NOT CONTAIN INPUT: " + node_name.str());
              not_found_error = true;
            }
            if (t->input_shifts[i] < 0 && t->input_shifts[i] < max_neg_shift)
              max_neg_shift = t->input_shifts[i];
          }
        }

        max_neg_shift = abs(max_neg_shift);
        if (not_found_error) continue;
        if (inputs_notnull == 0)
        {
          errorList.push_back("ALL INPUTS NULL: " + node_name.str());
          continue;
        }
        else if (inputs_notnull == 1)
        {
          errorList.push_back("INPUT NULL: " + node_name.str());
          continue;
        }

        if (t->input_is_negative.size() != t->inputs.size())
        {
          errorList.push_back("ISNEGATIVE VECTOR SIZE DOESNT MATCH INPUT SIZE: " + node_name.str());
          continue;
        }

        //check output values
        for (unsigned int cfg = 0, cfg_end = t->output_factor.size(); cfg < cfg_end; ++cfg)
        {
          for (unsigned int v = 0, v_end = t->output_factor[cfg].size(); v < v_end; ++v)
          {
            int64_t computed_output = 0;
            for (unsigned int i = 0, i_end = t->inputs.size(); i < i_end; ++i)
            {
              if (t->inputs[i] != NULL)
              {
                if (t->input_is_negative[i])
                {
                  computed_output -= (t->inputs[i]->output_factor[cfg][v] << (t->input_shifts[i] + max_neg_shift));
                }
                else
                {
                  computed_output += (t->inputs[i]->output_factor[cfg][v] << (t->input_shifts[i] + max_neg_shift));
                }
              }
            }
            if (t->output_factor[cfg][v] != DONT_CARE && t->output_factor[cfg][v] != 0)
            {
              if ((computed_output >> max_neg_shift) != t->output_factor[cfg][v])
              {
                stringstream tmp;
                tmp << "\t(computed:" << (computed_output >> abs(max_neg_shift)) << " != given:"
                    << t->output_factor[cfg][v]
                    << ")";
                errorList.push_back("OUTPUTFACTOR MISSMATCH 1: " + node_name.str() + tmp.str());
                continue;
              }
            }
          }
        }
      }
      else if (is_a<mux_node_t>(*(*it)))
      {
        if (!quiet) cout << "checking mux node" << node_name.str() << endl;
        mux_node_t *t = (mux_node_t *) (*it);
        int muxInputCount = 0;
        bool not_found_error = false;
        for (unsigned int i = 0, i_end = t->inputs.size(); i < i_end; ++i)
        {
          if (t->inputs[i] != NULL)
          {
            muxInputCount++;
            if (std::find(nodes_list.begin(), nodes_list.end(), t->inputs[i]) == nodes_list.end())
            {
              errorList.push_back("NODE_LIST DOES NOT CONTAIN INPUT: " + node_name.str());
              not_found_error = true;
            }
          }
        }
        if (not_found_error) continue;
        if (muxInputCount == 0)
        {
          errorList.push_back("ALL INPUTS NULL: " + node_name.str());
          continue;
        }
        else
        {
          for (unsigned int cfg = 0, cfg_end = t->output_factor.size(); cfg < cfg_end; ++cfg)
          {
            for (unsigned int v = 0, v_end = t->output_factor[cfg].size(); v < v_end; ++v)
            {
              if (cfg >= t->inputs.size())
              {
                errorList.push_back("INPUT SIZE MISSMATCH: " + node_name.str() + " (size is " + to_string(t->inputs.size()) + ", should be " + to_string(t->output_factor.size()) + ")");
              }
              else
              {
                if (t->inputs[cfg] == NULL)
                {
                  if (t->output_factor[cfg][v] != 0 && t->output_factor[cfg][v] != DONT_CARE)
                  {
                    cout << "t=" << t << endl;
                    errorList.push_back("OUTPUTFACTOR MISSMATCH INPUT NULL: " + node_name.str() + " for configuration " + to_string(cfg) + " and input " + to_string(v));
                  }
                }
                else
                {
                  if ((t->inputs[cfg]->output_factor[cfg][v] << t->input_shifts[cfg]) != t->output_factor[cfg][v] && t->output_factor[cfg][v] != DONT_CARE)
                  {

                    stringstream tmp;
                    tmp << "\t(" << t->inputs[cfg]->output_factor[cfg][v] << " != " << t->input_shifts[cfg] << ")";
                    errorList.push_back("OUTPUTFACTOR MISSMATCH 2: " + node_name.str() + " cfg " + to_string(cfg) + " input " + to_string(v) + tmp.str());
                  }
                }
              }
            }
          }
        }
      }
      else if (is_a<conf_adder_subtractor_node_t>(*(*it)))
      {
        if (!quiet) cout << "checking conf. add/sub node" << node_name.str() << endl;
        conf_adder_subtractor_node_t *t = (conf_adder_subtractor_node_t *) (*it);

        bool not_found_error = false;
        int max_neg_shift = 0;
        for (unsigned int i = 0, i_end = t->inputs.size(); i < i_end; ++i)
        {
          if (t->inputs[i] != NULL)
          {
            if (std::find(nodes_list.begin(), nodes_list.end(), t->inputs[i]) == nodes_list.end())
            {
              errorList.push_back("NODE_LIST DOES NOT CONTAIN INPUT: " + node_name.str());
              not_found_error = true;
            }
            if (t->input_shifts[i] < 0 && t->input_shifts[i] < max_neg_shift)
              max_neg_shift = t->input_shifts[i];
          }
        }
        max_neg_shift = abs(max_neg_shift);
        if (not_found_error) continue;

        if (t->input_is_negative.size() != t->output_factor.size())
        {
          errorList.push_back("OUTER ISNEGATIVE VECTOR SIZE DOESNT MATCH INPUT SIZE: " + node_name.str());
          continue;
        }

        for (unsigned int cfg = 0; cfg < t->input_is_negative.size(); cfg++)
        {
          if (t->input_is_negative[cfg].size() != t->inputs.size())
          {
            errorList.push_back("INNER ISNEGATIVE VECTOR SIZE DOESNT MATCH INPUT SIZE: " + node_name.str() + ", ISNEGATIVE VECTOR SIZE=" + to_string(t->input_is_negative[cfg].size()) + ", INPUT SIZE=" + to_string(t->inputs.size()));
            continue;
          }

        }

        for (unsigned int cfg = 0, cfg_end = t->output_factor.size(); cfg < cfg_end; ++cfg)
        {
          for (unsigned int v = 0, v_end = t->output_factor[cfg].size(); v < v_end; ++v)
          {
            //int64_t neg_shift=0;
            int64_t computed_output = 0;
            for (unsigned int i = 0, i_end = t->inputs.size(); i < i_end; ++i)
            {
              if (t->inputs[i] != NULL)
              {
                if (t->input_is_negative[cfg][i])
                {
                  computed_output -= (t->inputs[i]->output_factor[cfg][v] << (t->input_shifts[i] + max_neg_shift));
                }
                else
                {
                  computed_output += (t->inputs[i]->output_factor[cfg][v] << (t->input_shifts[i] + max_neg_shift));
                }
              }
            }
            if ((computed_output >> abs(max_neg_shift)) != t->output_factor[cfg][v] &&
                t->output_factor[cfg][v] != DONT_CARE)
            {
              stringstream tmp;
              tmp << "\t(computed:" << (computed_output >> abs(max_neg_shift)) << " != given:" << t->output_factor[cfg][v]
                  << ")";
              errorList.push_back("OUTPUTFACTOR MISSMATCH 3: " + node_name.str() + tmp.str());
            }
          }
        }
      }
      else if (is_a<register_node_t>(*(*it)) || is_a<output_node_t>(*(*it)))
      {
        if (!quiet) cout << "checking register/output node" << node_name.str() << endl;
        register_node_t *t = (register_node_t *) (*it);
        for (unsigned int cfg = 0, cfg_end = t->output_factor.size(); cfg < cfg_end; ++cfg)
        {
          for (unsigned int v = 0, v_end = t->output_factor[cfg].size(); v < v_end; ++v)
          {
            if (t->input == NULL)
            {
              errorList.push_back("INPUT OF REG IS NULL: " + node_name.str());
            }
            else
            {
              if (std::find(nodes_list.begin(), nodes_list.end(), t->input) == nodes_list.end())
              {
                errorList.push_back("NODE_LIST DOES NOT CONTAIN INPUT: " + node_name.str());
                v = v_end;
                cfg = cfg_end;
                break;
              }

            }
            if (t->output_factor[cfg][v] != DONT_CARE && t->output_factor[cfg][v] != 0)
            {
//              if ((t->input->output_factor[cfg][v] << t->input_shift) != t->output_factor[cfg][v])
              if (abs(t->input->output_factor[cfg][v] << t->input_shift) != abs(t->output_factor[cfg][v]))
              {
                errorList.push_back("OUTPUTFACTOR MISSMATCH 4: " + node_name.str() + " " + to_string(t->input->output_factor[cfg][v]) + "<<" + to_string(t->input_shift) + " != " + to_string(t->output_factor[cfg][v]));
              }
            }
          }
        }
      }
    }

    if (!quiet) cout << "checking configurations" << endl;
    if (configurationCountMap.size() > 1)
    {
      stringstream error;
      error << "CONFIGURATION COUNT MISSMATCH:" << endl;
      for (map<int, vector<string> >::iterator pr1 = configurationCountMap.begin(), pr1_end = configurationCountMap.end();
           pr1 != pr1_end;
           ++pr1)
      {
        error << "\tNODES WITH " << (*pr1).first << " CFG(S):" << endl;
        for (unsigned int i = 0, i_end = (*pr1).second.size(); i < i_end; ++i)
        {
          error << "\t\t" << (*pr1).second[i] << endl;
        }
      }
      errorList.push_back(error.str());
    }
    if (!quiet) cout << "checking input counts" << endl;
    if (inputCountMap.size() > 1)
    {
      stringstream error;
      error << "INPUT COUNT MISSMATCH:" << endl;
      for (map<int, vector<string> >::iterator pr1 = inputCountMap.begin(), pr1_end = inputCountMap.end();
           pr1 != pr1_end;
           ++pr1)
      {
        error << "\tNODES WITH " << (*pr1).first << " INPUT(S):" << endl;
        for (unsigned int i = 0, i_end = (*pr1).second.size(); i < i_end; ++i)
        {
          error << "\t\t" << (*pr1).second[i] << endl;
        }
      }
      errorList.push_back(error.str());
    }

    if (errorList.size() > 0)
    {
      for (list<string>::iterator iter = errorList.begin(), iter_end = errorList.end(); iter != iter_end; ++iter)
      {
        cout << *iter << endl;
      }
      throw runtime_error("Parse check failed");
    }
  }

  bool adder_graph_t::normalize_node(adder_subtractor_node_t *node)
  {
    bool allNodesNormalized=true;
    if(is_a<adder_subtractor_node_t>(*node))
    {
      adder_subtractor_node_t &anode = *((adder_subtractor_node_t *) node);
      if(anode.input_is_negative.size() >= 2)
      {
        if(anode.input_is_negative[0])
        {
          int indexOfNextPositiveInput = 1;
          while(anode.input_is_negative[indexOfNextPositiveInput])
          {
            if(indexOfNextPositiveInput < anode.input_is_negative.size()-1)
              indexOfNextPositiveInput++;
            else
            {
              if(!quiet) cout << "Found adder " << anode.output_factor << " in stage " << anode.stage << " with input 0 negative, swapping with other inputs is impossible as they are also negative (this may be more costly implementing those)" << endl;
              allNodesNormalized=false;
              indexOfNextPositiveInput = -1; // -1 means no positive input found (costly to realize but may happen)
            }
          }
          if(indexOfNextPositiveInput > 0)
          {
            if(!quiet) cout << "Found adder " << anode.output_factor << " in stage " << anode.stage << " with input 0 negative, swap it with input " << indexOfNextPositiveInput << " which is positive" << endl;

            //swap nodes:
            bool tmp_input_is_negative = anode.input_is_negative[indexOfNextPositiveInput];
            adder_graph_base_node_t *tmp_input = anode.inputs[indexOfNextPositiveInput];
            int64_t tmp_shift = anode.input_shifts[indexOfNextPositiveInput];

            anode.input_is_negative[indexOfNextPositiveInput] = anode.input_is_negative[0];
            anode.inputs[indexOfNextPositiveInput] = anode.inputs[0];
            anode.input_shifts[indexOfNextPositiveInput] = anode.input_shifts[0];

            anode.input_is_negative[0] = tmp_input_is_negative;
            anode.inputs[0] = tmp_input;
            anode.input_shifts[0] = tmp_shift;
          }
        }
      }
      else throw runtime_error("Found adder with less than 2 inputs, this does not make sense");
    }
    else if(is_a<conf_adder_subtractor_node_t>(*node))
    {
      conf_adder_subtractor_node_t &cnode = *((conf_adder_subtractor_node_t *) node);

      //search in all configurations for a negative input:
      bool foundNegativeFirstInput=false;
      for(int c = 0; c < cnode.input_is_negative.size(); c++)
      {
        if(cnode.input_is_negative[c].size() >= 2)
        {
          if(cnode.input_is_negative[c][0])
          {
            foundNegativeFirstInput=true;
          }
        }
      }

      //if negative input is found, try to normalize by swapping the inputs
      if(foundNegativeFirstInput)
      {
        if(!quiet) cout << "Found configurable adder " << cnode.output_factor << " in stage " << cnode.stage << " with input 0 negative" << endl;

        //now search for the next input that is positive in all configurations (if existing)
        int indexOfNextPositiveInput = 1;
        bool foundNegativeNextInput;
        while(indexOfNextPositiveInput < cnode.input_is_negative[0].size())
        {
          foundNegativeNextInput=false;
          for(int c = 0; c < cnode.input_is_negative.size(); c++)
          {
            if(cnode.input_is_negative[c][indexOfNextPositiveInput])
            {
              foundNegativeNextInput=true;
            }
          }
          if(!foundNegativeNextInput)
          {
            break; //done, we found an input for which all the configurations are positive
          }

          indexOfNextPositiveInput++;
        }

        if(indexOfNextPositiveInput >= cnode.input_is_negative[0].size())
        {
          allNodesNormalized=false;
          cout << "All other inputs are negative, swapping with other input is impossible (this may be more costly implementing those)" << endl;
        }
        else
        {
          if(!quiet) cout << "Found configurable adder with input 0 negative, swap it with input " << indexOfNextPositiveInput << " which is positive" << endl;

          //swap nodes:
          for(int ci = 0; ci < cnode.input_is_negative.size(); ci++)
          {
            bool tmp_input_is_negative = cnode.input_is_negative[ci][indexOfNextPositiveInput];
            cnode.input_is_negative[ci][indexOfNextPositiveInput] = cnode.input_is_negative[ci][0];
            cnode.input_is_negative[ci][0] = tmp_input_is_negative;
          }
          adder_graph_base_node_t *tmp_input = cnode.inputs[indexOfNextPositiveInput];
          cnode.inputs[indexOfNextPositiveInput] = cnode.inputs[0];
          cnode.inputs[0] = tmp_input;

          int64_t tmp_shift = cnode.input_shifts[indexOfNextPositiveInput];
          cnode.input_shifts[indexOfNextPositiveInput] = cnode.input_shifts[0];
          cnode.input_shifts[0] = tmp_shift;
        }
      }
    }
    return allNodesNormalized;
  }

  bool adder_graph_t::normalize_graph()
  {
    bool allNodesNormalized=true;
    if(!quiet) cout << "normalizing graph ..." << endl;
    for(auto node: nodes_list)
    {
      if(is_a<adder_subtractor_node_t>(*node) || is_a<conf_adder_subtractor_node_t>(*node))
      {
        allNodesNormalized &= normalize_node((adder_subtractor_node_t*) node);
      }
    }
    return allNodesNormalized;
  }

  string adder_graph_t::convert_old_syntax(string commandLine)
  {
    bool replace_active = false;
    for (unsigned int i = 0, i_end = commandLine.size(); i < i_end; ++i)
    {
      if (commandLine[i] == '[') replace_active = true;
      if (commandLine[i] == ']') replace_active = false;
      if (replace_active && commandLine[i] == ';') commandLine[i] = ',';
    }

    string zeroString;
    {
      int p1 = commandLine.find("["), p2 = commandLine.find("]");
      string firstOutput = commandLine.substr(p1 + 1, p2 - p1 - 1);

      int cmm_size = 1;
      for (unsigned int f = 0, f_end = firstOutput.size(); f < f_end; ++f)
      {
        if (firstOutput[f] == ',') cmm_size++;
      }
      zeroString = "[";
      for (int i = 0; i < cmm_size; ++i)
      {
        if (i > 0) zeroString += ",";
        zeroString += "0";
      }
      zeroString += "]";
    }

    stringstream new_k_adder_graph;
    new_k_adder_graph << "{";

    commandLine = commandLine.substr(commandLine.find("{") + 1);
    commandLine = commandLine.substr(0, commandLine.find_last_of("}"));
    bool first = true;
    while (commandLine.find("{") != string::npos)
    {
      if (!first)
      { new_k_adder_graph << ","; }
      first = false;
      new_k_adder_graph << "{";
      int i1 = commandLine.find("{") + 1, i2 = commandLine.find("}");
      string node_string = commandLine.substr(i1, i2 - i1);

      if (node_string.find(zeroString.c_str()) != string::npos)
      {
        int i3 = node_string.find(zeroString.c_str());
        new_k_adder_graph << "'R',";
        string tmp_node_string = node_string.substr(0, i3 - 1);
        if (tmp_node_string.substr(tmp_node_string.find_last_of(",") + 1) == "0")
        {
          new_k_adder_graph << tmp_node_string.substr(0, tmp_node_string.find_last_of(","));
        }
        else
        {
          new_k_adder_graph << tmp_node_string;
        }
      }
      else
      {
        new_k_adder_graph << "'A',";
        new_k_adder_graph << node_string;
      }
      new_k_adder_graph << "}";
      commandLine = commandLine.substr(i2 + 1);
    }
    new_k_adder_graph << "}";
    return new_k_adder_graph.str();
  }

/* use ::writesyn() instead !!
  string adder_graph_t::get_adder_graph_as_string()
  {
    stringstream graphStr;

    graphStr << "{";

    bool firstIteration = true;
    for (adder_graph_base_node_t *node: nodes_list)
    {
      if (is_a<input_node_t>(*node))
      {
        //nothing to do
      }
      else
      {
        if (!firstIteration)
        {
          graphStr << ",";
        }
        firstIteration = false;

        graphStr << "{";
        if (is_a<adder_subtractor_node_t>(*node))
        {
          graphStr << "'A',";
          graphStr << "[" << node->output_factor << "],";
          graphStr << node->stage << ",";
          int i = 0;
          bool firstInInnerIteration = true;
          for (adder_graph_base_node_t *input_node: ((adder_subtractor_node_t *) node)->inputs)
          {
            if (!firstInInnerIteration)
            {
              graphStr << ",";
            }
            graphStr << "[" << (((adder_subtractor_node_t *) node)->input_is_negative[i] ? "-" : "") << input_node->output_factor << "],";
            graphStr << input_node->stage << ",";
            graphStr << ((adder_subtractor_node_t *) node)->input_shifts[i];
            i++;
            firstInInnerIteration = false;
          }
        }
        else if (is_a<output_node_t>(*node) || is_a<register_node_t>(*node))
        {
          if (is_a<output_node_t>(*node))
            graphStr << "'O',";
          else
            graphStr << "'R',";

          graphStr << "[" << node->output_factor << "],";
          graphStr << node->stage << ",";
          graphStr << "[" << ((register_node_t *) node)->input->output_factor << "],";
          graphStr << ((register_node_t *) node)->input->stage << ",";
          graphStr << ((register_node_t *) node)->input_shift;
        }
        else if (is_a<mux_node_t>(*node))
        {
          graphStr << "'M',";

          graphStr << "[" << node->output_factor << "],";
          graphStr << node->stage << ",";
          int noOfConfigurations = node->output_factor[0].size();
//          vector<int64_t> shifts(noOfConfigurations);
          bool firstInInnerIteration = true;
          int i = 0;
          for (adder_graph_base_node_t *input_node: ((mux_node_t *) node)->inputs)
          {
            if (!firstInInnerIteration)
            {
              graphStr << ",";
            }
            graphStr << "[" << input_node->output_factor << "],";
            graphStr << input_node->stage << ",";
            graphStr << ((mux_node_t *) node)->input_shifts[i];
            i++;
            firstInInnerIteration = false;
          }
//          throw runtime_error("output for MUX node type not supported so far, sorry");
        }
        else if (is_a<conf_adder_subtractor_node_t>(*node))
        {
          throw runtime_error("output for conf. add/subtract node type not supported so far, sorry");
        }
        graphStr << "}";
      }
    }
    graphStr << "}";
    return graphStr.str();
  }
*/

  std::ostream &operator<<(std::ostream &stream, const std::vector<std::vector<int64_t> > &matrix)
  {
    for (int r = 0; r < matrix.size(); r++)
    {
      for (int c = 0; c < matrix[r].size(); c++)
      {
        if (matrix[r][c] == DONT_CARE)
          stream << "-";
        else
          stream << matrix[r][c];
        if (c < matrix[r].size() - 1) stream << ",";
      }
      if (r < matrix.size() - 1) stream << ";";
    }
    return stream;
  }

}

std::ostream &operator<<(std::ostream &stream, PAGSuite::adder_graph_t &adder_graph)
{
//  stream << adder_graph.get_adder_graph_as_string();
  adder_graph.writesyn(stream);
  return stream;
}