#pragma once

#include "gurobi_c++.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <tuple>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <new>
#include <vector>
#include <map>
#include <tuple>
#include <iterator>
#include <utility> 
#include <cstddef>
#include <type_traits>
#include <chrono>
#include <cfloat>
#include <time.h>
#include <cstring>
#include <math.h>
#include <cmath>
#include <ctgmath>
#include <stdlib.h>

#include <lemon/list_graph.h>
#include <lemon/network_simplex.h>
//#include <lemon/capacity_scaling.h>
#include <lemon/concepts/digraph.h>
#include <lemon/time_measure.h>
#include <lemon/adaptors.h>
#include <lemon/dijkstra.h>
#include <lemon/bfs.h>
#include <lemon/maps.h>
#include <lemon/cycle_canceling.h>
#include <lemon/capacity_scaling.h>
#include <lemon/cost_scaling.h>
#include <lemon/tolerance.h>
#include <lemon/preflow.h>

/* ************************************************************* */

#define comparator
#ifndef comparator
	//#define Primera_Vez
#endif

#define maximization_problem


/* ************************************************************* */


const double EPS1 = 1E-1;
const double EPS2 = 1E-2;
const double EPS3 = 1E-3;
const double EPS4 = 1E-4;
const double EPS5 = 1E-5;
const double EPS6 = 1E-6;
const double EPS7 = 1E-7;
const double EPS8 = 1E-8;

using namespace lemon;
using namespace lemon::concepts;
using namespace std;

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;


static int feas (0);
static int opti (0);
static int feas_2 (0);
static int opti_2 (0);


typedef duration<double> seconds_type;
//typedef chrono::duration<double,ratio<1,1000>> milliseconds_type;
typedef duration<int,std::milli> milliseconds_type;
typedef duration<double,ratio<3600,1>> hours_type;

seconds_type subproblem_time(0);

#ifndef vectors_typedef
#define vectors_typedef

typedef   ListDigraph::Node   						 Node;
typedef   ListDigraph::Arc     						 Arc;
typedef   ListDigraph::ArcMap<double>				 ArcDou;
typedef   ListDigraph::NodeMap<double>				 NodeDou;
typedef   ListDigraph::ArcMap<int>			    	 ArcInt;
typedef   ListDigraph::ArcMap<bool>			    	 ArcBool;
typedef   ListDigraph::NodeMap<int>			    	 NodeInt;
typedef   ListDigraph::NodeMap<bool>			     NodeBool;
typedef	  ListDigraph::OutArcIt 					 OutArcIt;
typedef   ListDigraph::NodeIt 						 NodeIt;
typedef   ReverseDigraph<ListDigraph>::NodeIt 		 NodeIt_rev;
typedef   ReverseDigraph<ListDigraph>::NodeMap<bool> NodeBool_rev;
typedef   ListDigraph::ArcIt 					     ArcIt;


//typedef  typename ReverseGraph<ListDigraph>::NodeIt 						 NodeIt_r;
//typedef  typename ReverseGraph<ListDigraph>::ArcIt 					         ArcIt_r;

using NS = NetworkSimplex<ListDigraph, double, double>;
using NS_F = NetworkSimplex<FilterArcs<ListDigraph>, double, double>;

using MaxFlow = Preflow<SubDigraph<ListDigraph>, SubDigraph<ListDigraph>::ArcMap<int>>;

typedef vector<vector<vector<vector<double> > > >  				 vec4Dou;
typedef vector<vector<vector<double> > >              			 vec3Dou;
typedef vector<vector<double> >                           	     vec2Dou;
typedef vector<double>                                           vec1Dou;
typedef vector<vector<vector<vector<int> > > > 					 vec4Int;
typedef vector<vector<vector<int> > >              				 vec3Int;
typedef vector<vector<int> >                        		     vec2Int;
typedef vector<int>                                      	     vec1Int;

typedef vec1Int::iterator                                          		 iteInt;
typedef vec1Int::const_iterator                                    		 citeInt;
typedef vec1Int::reverse_iterator                                  		 riteInt;
typedef vec1Int::const_reverse_iterator                            		 criteInt;
typedef vec1Dou::iterator                                          		 iteDou;
typedef vec1Dou::const_iterator                                    		 citeDou;
typedef vec1Dou::reverse_iterator                                  		 riteDou;
typedef vec1Dou::const_reverse_iterator                            		 criteDou;



typedef tuple<int,int,int,int>                  	             tup4Int;
typedef vector<tup4Int>                	                 vec_tup4Int;
typedef vector< pair<int,int> >                	            vec_pair;

typedef vector<GRBVar>									         vec1GRBvar;
typedef vector<vec1GRBvar >   			                    	 vec2GRBVar;
typedef vector<vec2GRBVar >  				                     vec3GRBVar;
typedef vector<vec3GRBVar>  				                     vec4GRBVar;


typedef vec_pair::const_iterator                        vec_pair_cit;
typedef vec_pair::iterator                              vec_pair_it;
typedef vec_tup4Int::const_iterator						vec_tup4Int_cit;
typedef vec_tup4Int::iterator							vec_tup4Int_it;


typedef int							Veh;



 #endif

#ifndef DVAP
#define DVAP


struct Instance_data
{
    int I;          
	int J;          
    int T;          
    int V;          
	vec2Int  tau;   
	vec3Dou  cos;   
	vec3Dou  cosf;  
	vec3Dou  pro;   
	vec3Int  m;	  	
	vec3Int  dem;	
	vec3Int  A;	    
	vec3Dou  In;
	vec2Int  Kap;
	
	string name;
	


	

	
};



#endif

static ofstream resultados_OV( "resultados_OV.txt", ios::app);
static ofstream op_sol_LP;
static ofstream op_sol_IP;
static ofstream op_solV_LP;
static ofstream op_solV_IP;


unsigned semilla;

static double sub_runtime_both(0.0);
static double sub_runtime_1st(0.0);
static double sub_runtime_2nd(0.0);

//static double sub_runtime_div(0.0);
//static double sub_runtime_cos(0.0);
//static double sub_runtime_run(0.0);
//static double sub_runtime_cut(0.0);

//static double sub_runtime_1(0.0);
//static double sub_runtime_2(0.0);
//static double sub_runtime_3(0.0);
//static double sub_runtime_4(0.0);
	
int DVAP_read_instance(
	int argc, 
	char *argv[],  
	Instance_data &instance);
  
	
bool Generate_Instance(
	int argc, 
	char *argv[], 
	Instance_data &instance);

//void set_network( Instance_data &instance);

inline int flat_ijtv(const Instance_data &ins, const int i, const int j, const int t, const int v ){
		
	return (i + (j*ins.I) + (t*ins.I*ins.J) + (v*ins.I*ins.J*ins.T));
	
	
}
inline int flat_itv(const Instance_data &ins, const int i, const int t, const int v ){
		
	return (i + (t*ins.I) + (v*ins.I*ins.T));
	
	
}
inline void ijtv_unflat(const Instance_data &ins, const int ind, int &i, int &j, int &t, int &v  ){
	
	
	i = ind%ins.I;
	j = ((ind-i)/ins.I)%ins.J;
	t = ((ind-(j*ins.I)-i)/(ins.I*ins.J))%ins.T;
	v = ((ind-(t*ins.I*ins.J)-(j*ins.I) - i)/(ins.J*ins.I*ins.T))%ins.V;
	
}
inline void itv_unflat(const Instance_data &ins, const int ind, int &i, int &t, int &v  ){
	
	
	i = ind%ins.I;
	t = ((ind-i)/ins.I)%ins.T;
	v = ((ind-(t*ins.I)-i)/(ins.I*ins.T))%ins.V;
	
}
inline int flat_ixt( const Instance_data &ins, const int t, const int i ){
		
	return (  (ins.I * t)+i < (ins.I*(ins.T+1)) ? (ins.I * t)+i:(ins.I * ins.T)+i  );
	
}
inline void ixt_flat(const Instance_data &ins, const int ind, int &t,  int &i ){
		
	if( ind >= (ins.I*(ins.T+1))){ cerr << " -- Index not within bounds to be unflattened --  " << flush; exit(EXIT_FAILURE);}
	i=ind%ins.I;
	t=ind/ins.I;
}
inline int flat_to_jit(const Instance_data &ins, const int t, const int i, const int j ){
	
	return j + (i*ins.J) + (t*ins.I*ins.J);
	
}
inline void jit_to_unflat(const Instance_data &ins, const int ind, int &t, int &i, int &j){
		
	j = ind%ins.J;
	i = ((ind-j)/ins.J)%ins.I;
	t = ((ind-(i*ins.J)-j)/(ins.J*ins.I))%ins.T;	
}
void print_status( GRBModel &model ){
	
	switch (model.get(GRB_IntAttr_Status)) {
	case 1:
		cerr << "\nLOADED: Model is loaded, but no solution information is available." << endl;
	break;
	case 2:		
		cerr << "\nOPTIMAL: Model was solved to optimality (subject to tolerances), and an optimal solution is available." << endl;	
	break;
	case 3:
		cerr << "\nINFEASIBLE: Model was proven to be infeasible." << endl;	
	break;
	case 4:
		cerr << "\nINF_OR_UNBD: Model was proven to be either infeasible or unbounded." << endl;	
	break;
	case 5:
		cerr << "\nUNBOUNDED: Model was proven to be unbounded." << endl;	
	break;
	case 6:
		cerr << "\nCUTOFF: Optimal objective for model was proven to be worse than the value specified in the Cutoff parameter. No solution information is available." << endl;	
	break;
	case 7:
		cerr << "\nITERATION_LIMIT: Optimization terminated because the total number of simplex iterations performed exceeded the value specified in the IterationLimit parameter, or because the total number of barrier iterations exceeded the value specified in the BarIterLimit parameter." << endl;	
	break;
	case 8:
		cerr << "\nNODE_LIMIT: Optimization terminated because the total number of branch-and-cut nodes explored exceeded the value specified in the NodeLimit parameter." << endl;	
	break;
	case 9:
		cerr << "\nTIME_LIMIT: Optimization terminated because the time expended exceeded the value specified in the TimeLimit parameter." << endl;	
	break;
	case 10:
		cerr << "\nSOLUTION_LIMIT: Optimization terminated because the number of solutions found reached the value specified in the SolutionLimit parameter." << endl;	
	break;
	case 11:
		cerr << "\nINTERRUPTED: Optimization was terminated by the user." << endl;	
	break;
	case 12:
		cerr << "\nNUMERIC: Optimization was terminated due to unrecoverable numerical difficulties." << endl;	
	break;
	case 13:
		cerr << "\nSUBOPTIMAL: Unable to satisfy optimality tolerances; a sub-optimal solution is available." << endl;	
	break;
	case 14:
		cerr << "\nINPROGRESS: An asynchronous optimization call was made, but the associated optimization run is not yet complete." << endl;	
	break;
	case 15:
		cerr << "\nUSER_OBJ_LIMIT: User specified an objective limit (a bound on either the best objective or the best bound), and that limit has been reached." << endl;	
	break;
	default:
		break;
	}
}


//inline int flat_ixt(Instance_data &instance, const int t, const int i ){
		
	//return (  (instance.i * t)+i < (instance.i*(instance.t+1)) ? (instance.i * t)+i:(instance.i * instance.t)+i  );
	
//}

template <typename DivMap >
class Visitor_AcSup : public BfsVisitor< const ReverseDigraph<ListDigraph> > {
public:
	Visitor_AcSup(const ReverseDigraph<ListDigraph>& _graph, DivMap& _dirMap, double& _sum_sup)
    : graph(_graph), dirMap(_dirMap), sum_sup(_sum_sup){
		cout << "\n --Calling constructor of Visitor_AcSup -- " << endl;
		sum_sup = 0.0;
	}
    
    void start (const Node &node){
		//cout << "\nstart Node: " << graph.id(node) << flush;
		sum_sup -= dirMap[node];
	}
    void reach (const Node &node){
		//cout << "\nReach Node: " << graph.id(node) << flush;
	}
    void process (const Node &node){
		//cout << "\nProcess Node: " << graph.id(node) << setw(5) << dirMap[node] << flush;
		sum_sup += dirMap[node];
	}
    void discover(const Arc& arc) {
		//cout << "\tDiscover Arc: " << graph.id(arc) << flush;
	}
	void examine(const Arc& arc) {
		//cout << "\tExamine Arc: " << graph.id(arc) << flush;
	}  
    
private:
  const ReverseDigraph<ListDigraph>& graph;  
  DivMap& dirMap;
  double& sum_sup;
};



