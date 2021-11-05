#include "DVAP_par.hpp"

//#define second_benders_BFVAP
/*Activate the code in which backlogging demand constraints are placed within the MP and taken out of the subporblem.*/

//#define sub_solver
/* Activate the code that uses GUROBI for solving the first type of subproblem both in:
 * - BFTBenders
 * - FleetBenders
 * When Deactivated, LEMON-library-graph is used for solving the first type of subproblem both in:
 * - BFTBenders
 * - FleetBenders
 *   In this case, a general graph is built outside the callback function and especializations are used
 *   by means of subgraphs with filted arcs for each vehicle type v. */

//#define solver
/* Activate the code for using GUROBI in solving the first type of subproblem both in:
 * - BFTBenders
 * - FleetBenders
 * When Deactivated, LEMON-library-graph is used for solving the first type of subproblem both in:
 * - BFTBenders
 * - FleetBenders
 *   In this case, one new graph is built for each subproblem of vehicle type v. 
 * NOTE: In this two cases, its Activation depends on the Activation of #define sub_solver      */


//#define primals
/* Activate the code that solves the primal subproblem (Maximization) for both types of subproblem:
 * - Sub I:  decides the repostioning and sizing of empty vehicles for each vehicle type v.
 * - Sub II: decides the movement of backlogged demand for each arc (i,j) along the planning horizon*/

//static ofstream cuts("cuts.txt",ios::app);


class BFTBendersSub: public GRBCallback
{
  public:
   
    const Instance_data &ins;
    const vec4GRBVar    &X;
    const vec1GRBvar    &u;
    #ifndef second_benders_BFVAP
    const vec2GRBVar    &uu;
    #endif
    ListDigraph         &grafo;
    NodeDou             &divergence;
	ArcDou              &costo     ;
	ArcBool             &filter;
	NodeDou             &po;    
	#ifndef second_benders_BFVAP
	vec1Dou             &diver;
	vec1Dou             &dual;
	#endif
    BFTBendersSub( const Instance_data    &ins,  
                   const vec4GRBVar       &X, 
                   const vec1GRBvar       &u 
                   #ifndef second_benders_BFVAP
                   ,const vec2GRBVar       &uu
                   #endif
                   ,ListDigraph 			  &grafo, 
				   NodeDou  			  &divergence, 
				   ArcDou  				  &costo, 
				   ArcBool  			  &filter, 
				   NodeDou  			  &potential
				   #ifndef second_benders_BFVAP
				   ,vec1Dou                &diver,
				   vec1Dou                &dual
				   #endif
				   )
    :ins(ins), X(X), u(u), 
    #ifndef second_benders_BFVAP 
		uu(uu), 
    #endif
    grafo(grafo), divergence(divergence), costo(costo), filter(filter), po(potential)
    #ifndef second_benders_BFVAP
		,diver(diver),dual(dual)
    #endif
    {
		  //cout << "\n\tI:" << ins.I << " J:" << ins.J << " T:" << ins.T << " V:" << ins.V << flush;
    }
    protected:
    void callback() {
      try {
		  
		  if (where == GRB_CB_MIPSOL /*&& where != GRB_CB_PRESOLVE*/) {
			  //cout << "\n -- Entering Callback -- " << flush;
			  
			  high_resolution_clock::time_point tbi = high_resolution_clock::now();
			  string var_name;  
			  int I = ins.I , J = ins.I , T = ins.T , V = ins.V,i,j,k,t;Veh v;
			  string sub;
			  
			  //cout << "\n\n\tI:" << ins.I << " J:" << ins.J << " T:" << ins.T << " V:" << ins.V << flush;
			  
			  double coe;
			  double ZLe(-DBL_MAX);
			  
			  for (v=0; v<V; v++){
				  
				  	  //cout << "\n\t --------- v = " << v << " ----------- " << flush;//while(getchar() != '\n');
					  
					  //vec3Dou X_sol(I,vec2Dou(J, vec1Dou(T,0.0)));//get rid of this at the end
					  //for (t = 0; t < T; t++)
						//for(i = 0; i < I ; i++)
							//for(j = 0; j < J; j++)
								//if( getSolution(X.at(i).at(j).at(t).at(v)) > EPS4 ){
									//X_sol.at(i).at(j).at(t) = getSolution(X.at(i).at(j).at(t).at(v));
									//cout << "\n\t" << X.at(i).at(j).at(t).at(v).get(GRB_StringAttr_VarName)  << " = " << getSolution(X.at(i).at(j).at(t).at(v)) << flush;
									//cout << setw(10) << flat_ixt(ins,t,i) << setw(10) << flat_ixt(ins,(i==j?t+1:t + ins.tau.at(i).at(j)),j) << flush;
								//}
								
					    ZLe= -DBL_MAX;
					    
					    //high_resolution_clock::time_point tb1 = high_resolution_clock::now();
					    
        				double supersup_div(0.0);//cout << "\nOK1" << flush;
        				for (t = 0; t < T; t++){
        					for(i = 0; i < I ; i++){
        							coe=ins.m.at(i).at(t).at(v);															
        							for (j = 0; j < J; j++) if(t<T) coe -= (getSolution(X.at(i).at(j).at(t).at(v)) > EPS4?getSolution(X.at(i).at(j).at(t).at(v)):0.0) ;																											
        							for (j = 0; j < J; j++) if((i != j) && ((t - ins.tau.at(j).at(i)) >= 0) ) coe += (getSolution(X.at(j).at(i).at(t - ins.tau.at(j).at(i)).at(v)) > EPS4?getSolution(X.at(j).at(i).at(t - ins.tau.at(j).at(i)).at(v)):0.0) ;
        						
        							//if(coe > EPS2 || coe < -EPS2 )
        							//cout << "\ncoe[" << flat_ixt(ins,t,i) /*<< "]-[" << i << "][" << t*/ << "]=" << coe << flush;
        							
        							if(std::round(coe) < -EPS4) supersup_div += fabs(std::round(coe));
        							//divergence[grafo.nodeFromId(flat_ixt(ins,t,i))] = std::round(coe);
        							divergence.set(grafo.nodeFromId(flat_ixt(ins,t,i)), std::round(coe));
        						
        							//if(fabs(divergence[grafo.nodeFromId(flat_ixt(ins,t,i))]-coe) > EPS3 ){
        								//cerr << "\nError in setting divergence in subproblem: " << v 
        									 //<< " as LEMON( " << divergence[grafo.nodeFromId(flat_ixt(ins,t,i))] 
        									 //<< " ) and COE( " << coe << " ) does not match " << flush;
        								//exit(EXIT_FAILURE);
        							//}
        							
        							//if(divergence[grafo.nodeFromId(flat_ixt(ins,t,i))] > EPS2 || divergence[grafo.nodeFromId(flat_ixt(ins,t,i))] < -EPS2 )
        							//cout << "\tdiv[" << flat_ixt(ins,t,i) /*<< "]-[" << i << "][" << t*/ << "]=" << divergence[grafo.nodeFromId(flat_ixt(ins,t,i))] << flush;
        									
        					}
        				}
        				divergence[grafo.nodeFromId((I*(T+1)))] = supersup_div;
        				//cout << "\nsupersup: " << supersup_div << endl;
        				
        				//for (ArcIt u(grafo); u != INVALID; ++u){
        					
        					//cout << "\nArcId: " << setw(5) << grafo.id(u)  
        										//<< setw(5) << grafo.id(grafo.source(u))  
        										//<< setw(5) << grafo.id(grafo.target(u))  
        										//<< setw(5) << costo[u]  
        										//<< flush;
        				//}cout << endl;
        															
        				
        				
        				//cout << "\nOK2" << flush;
        				//high_resolution_clock::time_point tb2 = high_resolution_clock::now();
        				
        				int id1,id2,ti,tj;
        				for (ArcIt u(grafo); u != INVALID; ++u){
        					
        					filter.set(u,true);
        					id1= grafo.id(grafo.source(u)); 
        					id2= grafo.id(grafo.target(u));
        					
        					//cout << "\nid1: " << setw(5) << id1 << "      id2: " << setw(5) << id2 << flush;
        					
        					if( (id1 < (I*(T+1))) && (id2 < (I*(T+1))) ){
        						
        						ixt_flat(ins,id1,ti,i );ixt_flat(ins,id2,tj,j );
        						
        						if( ins.A.at(i).at(j).at(v) == 0 ) filter.set(u,false);
        						costo.set(u, ins.cos.at(i).at(j).at(v) );
        						
        						
        					}else if( (id1 == (I*(T+1))) && (id2 < (I*(T+1))) ){
        						ixt_flat(ins,id2,tj,j );
        						costo.set(u, ins.cosf.at(j).at(tj).at(v) );
        					}else{
        						
        						cerr << "\nTHIS OPTION MUST NOT HAPPEN" << endl;exit(EXIT_FAILURE);
        						
        					}
        				
        					
        				}
        			
        				//cout << "\nOK3" << flush;
        				//high_resolution_clock::time_point tb3 = high_resolution_clock::now();
        				
        				FilterArcs<ListDigraph> subgrafo(grafo, filter); //high_resolution_clock::time_point tb3_1 = high_resolution_clock::now();
        				
        				
        				NS_F sol_NS(subgrafo); //high_resolution_clock::time_point tb3_2 = high_resolution_clock::now();
        				sol_NS.costMap(costo).supplyMap(divergence).supplyType(NS_F::LEQ); //high_resolution_clock::time_point tb3_3 = high_resolution_clock::now();
        				NS_F::ProblemType stat_NS = sol_NS.run();
        				
        				//high_resolution_clock::time_point tb4 = high_resolution_clock::now();
        				
						/*
						cerr << "\nOPTIMAL COST=" << sol_NS.totalCost() << endl; 		
						for (ArcIt n(subgrafo); n != INVALID; ++n){
							if( flo[n] != 0 ){
									Node x = subgrafo.source(n);
									Node y = subgrafo.target(n);
									cout << "\nflow[" << subgrafo.id(n) << "]=" << flo[n] << flush;
									cout << "\t[" << subgrafo.id(x) << "][" << subgrafo.id(y) << "]=" << flo[n] << " Costo: " << costo[n] << flush;
			
							}
						}
						for (NodeIt n(subgrafo); n != INVALID; ++n){
							if( po[n] != 0  ){
								cout << "\npot[" << subgrafo.id(n) << "]="<< setprecision(30) << po[n] << flush;
							}
						}*/			
					
						//cout << "\nOK4" << flush;
						if(stat_NS == NS_F::OPTIMAL ){
							
							ZLe = sol_NS.totalCost();
							//ArcDou flo(grafo);
							//sol_NS.flowMap(flo);
							
							sol_NS.potentialMap(po);
							
							if( ZLe > getSolution(u.at(v)) ){
								
								opti++;
								GRBLinExpr lhs = 0;
								for(int i=0;i<I;i++){
									for(int t=0;t<T;t++){
										//cout << "\nOut" << flush;
										for (j = 0; j < J; j++){
											if( po[subgrafo.nodeFromId(flat_ixt(ins,t,i))] > EPS3 )	
												lhs -= po[subgrafo.nodeFromId(flat_ixt(ins,t,i))]*X.at(i).at(j).at(t).at(v);
												
												
										}
										//cout << "\nIn" << flush;
										for (j = 0; j < J; j++){
											if((i != j) && ((t - ins.tau.at(j).at(i)) >= 0) ){
												if( po[subgrafo.nodeFromId(flat_ixt(ins,t,i))] > EPS3 )
													lhs += po[subgrafo.nodeFromId(flat_ixt(ins,t,i))]*X.at(j).at(i).at((t - ins.tau.at(j).at(i))).at(v);
											}
										}
									}
								}
								lhs += u.at(v);
								
								int rhs = 0;
								for(int i=0; i < I ; i++){
									for(int t=0;t < T;t++){
										rhs -= (ins.m.at(i).at(t).at(v)*po[subgrafo.nodeFromId(flat_ixt(ins,t,i))]);
									}
								}
								
								addLazy(lhs >= rhs);
								//cout << "\n\tOptimality Cut Added" << flush;
								//cout << "\n" << lhs << ">=" << rhs << flush;
								
								
							}
							
						}else{
							cerr << "\n The subproblem " << v << " is not optimal, hence there is an error. The model is " << flush;
							switch (stat_NS) {
							case NS_F::INFEASIBLE:
								cerr << "INFEASIBLE." << endl;
							break;
							case NS_F::UNBOUNDED:
								cerr << "UNBOUNDED." << endl;
							break;
							default:
								break;
							}
							exit(EXIT_FAILURE);
						}
						
						
						//high_resolution_clock::time_point tb5 = high_resolution_clock::now();
						
						//seconds_type runtime_div = duration_cast<seconds_type>(tb2 - tb1);
						//seconds_type runtime_cos = duration_cast<seconds_type>(tb3 - tb2);
						//seconds_type runtime_run = duration_cast<seconds_type>(tb4 - tb3);
						//seconds_type runtime_cut = duration_cast<seconds_type>(tb5 - tb4);
						
						//seconds_type runtime_1 = duration_cast<seconds_type>(tb3_1 - tb3);
						//seconds_type runtime_2 = duration_cast<seconds_type>(tb3_2 - tb3_1);
						//seconds_type runtime_3 = duration_cast<seconds_type>(tb3_3 - tb3_2);
						//seconds_type runtime_4 = duration_cast<seconds_type>(tb4 - tb3_3);
					   
					    //sub_runtime_div += runtime_div.count();
					    //sub_runtime_cos += runtime_cos.count();
					    //sub_runtime_run += runtime_run.count();
					    //sub_runtime_cut += runtime_cut.count();
					    
					    //sub_runtime_1 += runtime_1.count();
					    //sub_runtime_2 += runtime_2.count();
					    //sub_runtime_3 += runtime_3.count();
					    //sub_runtime_4 += runtime_4.count();
			  		
        	  }/*end_if v*/
        	  
        	  high_resolution_clock::time_point tbm = high_resolution_clock::now();
        		
        		#ifndef second_benders_BFVAP
        		/* *****************************************************************************************
        		 * 								Beginning of 2nd type subproblem
        		 * *****************************************************************************************/
        	    //while(getchar() != '\n');
        		for( i = 0; i < I;i++){
        			for (j = 0; j < J; j++){
        				if( i != j ){
        					
        					//cout << "\n\n\t --------- i=" << i << ",j=" << j  << " ----------- " << flush;//while(getchar() != '\n');
        				  
        					//vec3Dou X_sol(I,vec2Dou(J, vec1Dou(T,0.0)));//get rid of this at the end
        					//for ( v=0; v<V; v++)
        						//for (t=0; t<T; t++)
        							//if( getSolution(X.at(i).at(j).at(t).at(v)) > EPS4 ){
        								//X_sol.at(i).at(j).at(t) = getSolution(X.at(i).at(j).at(t).at(v));
        								//cout << "\n\t" << X.at(i).at(j).at(t).at(v).get(GRB_StringAttr_VarName)  << " = " << getSolution(X.at(i).at(j).at(t).at(v)) << flush;
        								//cout << setw(10) << flat_ixt(ins,t,i) << setw(10) << flat_ixt(ins,(i==j?t+1:t + ins.tau.at(i).at(j)),j) << flush;
        							//}

        					
        					
        					
							fill(diver.begin(),diver.end(),0.0);
							fill(dual.begin(),dual.end(),0.0);
							
							/* Modify objective function */
							for (t = 0; t < T; t++){ 
								coe=ins.dem.at(i).at(j).at(t);
								for ( v=0; v<V; v++ ) coe -= (getSolution(X.at(i).at(j).at(t).at(v)) > EPS4?getSolution(X.at(i).at(j).at(t).at(v)):0.0);
								//cout << "\n\tdiv[" << t << "]=" << setw(4) << std::round(coe) << flush;
								diver.at(t) = std::round(coe);
							}
							
							double soma_sup(0.0);int it;
							for(it=0; it < T; it++){
								soma_sup += diver.at(it);
								if(soma_sup < -EPS5){
									break;
								}
							}
							double ZLg(0.0);
							
							//cout << "\n\tsoma_sup: " << soma_sup << flush;
							//cout << "\n\tit: " << it << flush;
							
							if( soma_sup < -EPS5 || soma_sup > EPS5){
								if(soma_sup > EPS5){
									for(int iter(0); iter < T; iter++) dual.at(iter) = 1.0;
								}else{
									for(int iter(0); iter <= it; iter++) dual.at(iter) = -1.0;
								}
								//for(int iter(0); iter < T; iter++){
									//cout << "\n\trai[" << iter << "]=" << dual.at(iter) << flush;
								//} 
								
								// Incorporate feasibility cut.
								feas_2++;
								GRBLinExpr lhs = 0;
								for(int t=0;t<T;t++){
									for ( Veh v(0); v<V; v++ ){
										if( dual.at(t) > EPS3 || dual.at(t) < -EPS3 ){
											lhs += dual.at(t)*X.at(i).at(j).at(t).at(v);
										}
									}
								}
								
								//lhs += uu.at(i).at(j);
							
								double rhs = 0;
								for(int t=0;t < T;t++){
									if( dual.at(t) > EPS3 || dual.at(t) < -EPS3 )
										rhs += (double(ins.dem.at(i).at(j).at(t))*dual.at(t));
								}
							
								//cout << "\n\n\t" << lhs << " <= " << rhs << flush;
								//cuts << "\n\n\t" << lhs << " <= " << rhs << flush;
								addLazy(lhs >= rhs);
								//cout << "\n\tFeasibility Cut Added Inventory" << flush;
								
							}else{
								for(int iter(T-2); iter >= 0; iter--){
									dual.at(iter) = dual.at(iter+1) + ins.In.at(i).at(j).at(iter);
									ZLg += dual.at(iter)*diver.at(iter);
								}
								
								//cout << endl;
								//for(int iter(0); iter < T; iter++){
									//cout << "\n\tzol[" << iter << "]=" << dual.at(iter) << flush;
								//} 
								
								// Incorporate optimality cut.
								if( ZLg > getSolution(uu.at(i).at(j)) ){
									
									opti_2++;
									GRBLinExpr lhs = 0;
									for(int t=0;t<T;t++){
										for ( Veh v(0); v<V; v++ ){
											if( dual.at(t) > EPS4 || dual.at(t) < -EPS4){
												lhs += dual.at(t)*X.at(i).at(j).at(t).at(v);
											}
										}
									}
									
									lhs += uu.at(i).at(j);
								
									double rhs = 0;
									for(int t=0;t < T;t++){
										if( dual.at(t) > EPS4 || dual.at(t) < -EPS4)
											rhs += (double(ins.dem.at(i).at(j).at(t))*dual.at(t));
									}
								
									
									//cout << "\n\n\t" << lhs << " <= " << rhs << flush;
									//cuts << "\n\n\t" << lhs << " <= " << rhs << flush;
									addLazy(lhs >= rhs);
									//cout << "\n\tOptimality Cut Added Inventory" << flush;
									
									
								}
							}
						    
						    
        					
						}
        				

        				
        				
        			}/*end_if J*/
        		}/*end_if I*/
				#endif
        	  
        	   high_resolution_clock::time_point tbf = high_resolution_clock::now();
        	   seconds_type runtime_both = duration_cast<seconds_type>(tbf - tbi);
        	   seconds_type runtime_1st =  duration_cast<seconds_type>(tbm - tbi);
        	   seconds_type runtime_2nd =  duration_cast<seconds_type>(tbf - tbm);
        	   
			   sub_runtime_both += runtime_both.count();
			   sub_runtime_1st  += runtime_1st.count();
			   sub_runtime_2nd  += runtime_2nd.count();
        	   //cout << "\n -- Exiting Callback -- " << flush;
        	   //while(getchar() != '\n');
        	   //exit(0);
		 }
		  
	  }catch (GRBException e) {
        cout << "Error number: " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
        exit(EXIT_FAILURE);
      }catch (const exception &exc){
		cerr << "\nCallback - " << exc.what() << flush;
		exit(EXIT_FAILURE);
	  }catch (...) {
        cout << "Error during callback" << endl;
        exit(EXIT_FAILURE);
      }
    }
};

