#include "DVAP_par.hpp"


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
   
     Instance_data &ins;
    const vec4GRBVar    &X;
    const vec1GRBvar    &u;
    const vec2GRBVar    &uu;
    vec2GRBVar    &p;
    vec1GRBvar    &w;
    GRBModel &modelSubD1;
    GRBModel &modelSubD2;
    BFTBendersSub( Instance_data &ins,  const vec4GRBVar &X, const vec1GRBvar &u, const vec2GRBVar &uu, vec2GRBVar &p, vec1GRBvar &w, GRBModel &modelSubD1, GRBModel &modelSubD2)
    :ins(ins), X(X), u(u), uu(uu), p(p), w(w), modelSubD1(modelSubD1), modelSubD2(modelSubD2){
		  //cout << "\n\tI:" << ins.I << " J:" << ins.J << " T:" << ins.T << " V:" << ins.V << flush;
    }
  protected:
    void callback() {
      try {
		  
		   if (where == GRB_CB_MIPSOL /*&& where != GRB_CB_PRESOLVE*/) {
			  //cout << "\n -- Entering Callback -- " << flush;
			  
			  string var_name;  
			  int I = ins.I , J = ins.I , T = ins.T , V = ins.V,i,j,k,t;Veh v;
			  string sub;
			  
			  //cout << "\n\n\tI:" << ins.I << " J:" << ins.J << " T:" << ins.T << " V:" << ins.V << flush;
			  
			  double coe;
			  double ZLe(-DBL_MAX),ZD(-DBL_MAX);
			  
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
					 ZD = -DBL_MAX;
					 #ifdef primals
					 /* *********************************************************************************************/
						 GRBEnv envSubP = GRBEnv(true);
						 envSubP.set("LogFile", "Sub.log");
						 envSubP.start();
						 GRBModel modelSubP1 = GRBModel(envSubP);
						 modelSubP1.set(GRB_IntParam_OutputFlag, 0);
						 modelSubP1.set(GRB_StringAttr_ModelName, "P_MCFP");
						 modelSubP1.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
						 vec3GRBVar      y(I,vec2GRBVar(J,vec1GRBvar(T)));
						 vec2GRBVar      z(J,vec1GRBvar(T));
						 for(i = 0; i < I ; i++)
							for (t = 0; t < T; t++){
								var_name = "z("+ to_string(i)+ ","+to_string(t)+")";
								z.at(i).at(t) = modelSubP1.addVar(0.0, GRB_INFINITY, ins.cosf.at(i).at(t).at(v), GRB_CONTINUOUS, var_name);
								
								for(j = 0; j < J; j++){
									var_name = "y("+ to_string(i)+ ","+ to_string(j)+","+to_string(t)+")";
									y.at(i).at(j).at(t) = modelSubP1.addVar(0.0, ((ins.A.at(i).at(j).at(v)==0)?0.0:GRB_INFINITY), ins.cos.at(i).at(j).at(v), GRB_CONTINUOUS, var_name);
								}
							}
						
						 for(i = 0; i < I ; i++){
							for(t = 0; t < T; t++){
								
								coe=ins.m.at(i).at(t).at(v);
								for (j = 0; j < J; j++) if(t<T) coe -= getSolution(X.at(i).at(j).at(t).at(v)) ; 
								for (j = 0; j < J; j++) if((i != j) && ((t - ins.tau.at(j).at(i)) >= 0) ) coe += getSolution(X.at(j).at(i).at(t - ins.tau.at(j).at(i)).at(v)) ;
								
						
								GRBLinExpr sum_expr = 0;
								for (j = 0; j < J; j++) if(t<T) sum_expr += y[i][j][t] ;
								for (j = 0; j < J; j++) if((i != j) && ((t - ins.tau.at(j).at(i)) >= 0) ) sum_expr -= y[j][i][t - ins.tau.at(j).at(i)] ;
								if (t > 0) sum_expr -= y[i][i][t-1];
								sum_expr -= z[i][t];
								
								var_name = "flow("+ to_string(i)+","+to_string(t)+")";
								modelSubP1.addConstr(sum_expr, GRB_LESS_EQUAL , coe , var_name);
							}
						}

						modelSubP1.set(GRB_IntParam_Presolve, 0);
						sub = "subPI-"+to_string(v)+".lp";
						modelSubP1.write(sub);
						// Optimize model
						modelSubP1.optimize();
					 /* *********************************************************************************************/
					 #endif
					  
					  
					  sub = "subDI-"+to_string(v)+".lp";
					  modelSubD1.update();
					  //modelSubD1.write(sub);
					  //cout << "\nOK1" << flush;
					  //while(getchar() != '\n');
					  
					  #ifdef solver
					  /* Modify objective function */
					  for (t = 0; t < T; t++){
							for(i = 0; i < I ; i++){
								coe=ins.m.at(i).at(t).at(v);															
								for (j = 0; j < J; j++) if(t<T) coe -= (getSolution(X.at(i).at(j).at(t).at(v)) > EPS4?getSolution(X.at(i).at(j).at(t).at(v)):0.0) ;																											
								for (j = 0; j < J; j++) if((i != j) && ((t - ins.tau.at(j).at(i)) >= 0) ) coe += (getSolution(X.at(j).at(i).at(t - ins.tau.at(j).at(i)).at(v)) > EPS4?getSolution(X.at(j).at(i).at(t - ins.tau.at(j).at(i)).at(v)):0.0) ;
								
								p.at(i).at(t).set(GRB_DoubleAttr_Obj,coe);
								/* This is the additional constraint that accounts for the sizing decision in the dual*/
								p.at(i).at(t).set(GRB_DoubleAttr_LB, -ins.cosf.at(i).at(t).at(v));
							}
						}
						
						modelSubD1.update();
						//modelSubD1.write(sub);
						//cout << "\nOK2" << flush;
					    //while(getchar() != '\n');
						
						/* Subproblem constraints Definition */
						for(t = 0; t < T; ++t)
							for(i = 0; i < I; ++i)
								for(j = 0; j < J; ++j){
									if( ins.A.at(i).at(j).at(v) == 1 ){
										var_name = "cos("+to_string(i)+","+to_string(j)+","+to_string(t)+")";
										if( i == j ){
											modelSubD1.addConstr(((t+1)>=T?p[i][t]:p[i][t]-p[j][t+1]), GRB_LESS_EQUAL ,ins.cos.at(i).at(j).at(v), var_name);
										}else{
											modelSubD1.addConstr(((t+ins.tau[i][j])>=T?p[i][t]:p[i][t]-p[j][t+ins.tau[i][j]]), GRB_LESS_EQUAL ,ins.cos.at(i).at(j).at(v), var_name);
										
										}
									}
								}
								
						modelSubD1.update();		
								
						modelSubD1.set(GRB_IntParam_Presolve, 0);
						
						modelSubD1.write(sub);
						//cout << "\nOK3" << flush;
					    //while(getchar() != '\n');
						// Optimize model
						modelSubD1.optimize();
						
						if(modelSubD1.get(GRB_IntAttr_Status) == 2 ){
							
							//cout << "\n\tStatus Dual Solution: " << (modelSubD1.get(GRB_IntAttr_Status)==2?"Optimal":(modelSubD1.get(GRB_IntAttr_Status)==5?"Unbounded":(modelSubD1.get(GRB_IntAttr_Status)==3?"Infeasible":"none"))) << flush;
							ZD = modelSubD1.get(GRB_DoubleAttr_ObjVal);
							
							
							#ifdef primals
							cout << "\n\tStatus Primal Solution: " << (modelSubP1.get(GRB_IntAttr_Status)==2?"Optimal":(modelSubP1.get(GRB_IntAttr_Status)==5?"Unbounded":(modelSubP1.get(GRB_IntAttr_Status)==3?"Infeasible":"none"))) << flush;
							double ZP = modelSubP1.get(GRB_DoubleAttr_ObjVal);
							if( fabs(ZP-ZD) > EPS3 ){
								cerr << "\n\tError in solving optimal subproblem: " << v << " as ZP( " << ZP << " ) and ZD( " << ZD << " ) does not match " << flush;
								exit(EXIT_FAILURE);
							}
							#endif
						
							if( ZD > getSolution(u.at(v)) ){
								
								//cout << "\n\t    i  t  " << flush;
								//for(int t=0;t<T;t++)
									//for(int i=0;i<I;i++)
										//if( p[i][t].get(GRB_DoubleAttr_X) < -EPS3)
											//cout << "\n\tsol[" << i <<"][" << t << "] - (" << flat_ixt(ins,t,i ) << ") = " << p[i][t].get(GRB_DoubleAttr_X) << flush;
								
								opti++;
								GRBLinExpr lhs = 0;
								for(int i=0;i<I;i++){
									for(int t=0;t<T;t++){
										//cout << "\nOut" << flush;
										for (j = 0; j < J; j++){
											if( p[i][t].get(GRB_DoubleAttr_X) < -EPS4 ){
												lhs += p[i][t].get(GRB_DoubleAttr_X)*X.at(i).at(j).at(t).at(v);
											}
										}
										//cout << "\nIn" << flush;
										for (j = 0; j < J; j++){
											if((i != j) && ((t - ins.tau.at(j).at(i)) >= 0) ){
												if(p[i][t].get(GRB_DoubleAttr_X) < -EPS4){
													lhs -= p[i][t].get(GRB_DoubleAttr_X)*X.at(j).at(i).at((t - ins.tau.at(j).at(i))).at(v);
												}
											}
										}
									}
								}
								lhs += u.at(v);
								
								int rhs = 0;
								for(int i=0; i < I ; i++){
									for(int t=0;t < T;t++){
										rhs += (ins.m.at(i).at(t).at(v)*p[i][t].get(GRB_DoubleAttr_X));
									}
								}
								
								//cout << "\n\n\t" << lhs << " <= " << rhs << flush;
								//cuts << "\n\n\t" << lhs << " <= " << rhs << flush;
								addLazy(lhs >= rhs);
								//cout << "\n\tOptimality Cut Added Flux" << flush;
								
								
							}
						}
						if(modelSubD1.get(GRB_IntAttr_Status) == 5 ){
							
							cerr << "\nError: This should not be infeasible for this model." << endl;
							exit(EXIT_FAILURE);
							
						}
	
						
						/* Remove constraints to process next vehicle*/
						GRBConstr* array_ctr = 	modelSubD1.getConstrs ( );
						for (int k = 0; k < modelSubD1.get(GRB_IntAttr_NumConstrs); ++k)modelSubD1.remove(array_ctr[k]);
						modelSubD1.update();
						/*i=j=t=-1;
						for (int k = 0; k < modelSubD1.get(GRB_IntAttr_NumConstrs); ++k){
							cout << endl << array_ctr[k].get(GRB_StringAttr_ConstrName) << flush;
							jit_to_unflat(ins, k,t,i,j);
							array_ctr[k].set(GRB_DoubleAttr_RHS,ins.cos.at(i).at(j).at(v));
						}*/
						#else
						ListDigraph grafo;
						Node x;
						
						for( int t(0); t < ins.T+1; t++)
						  for ( int i(0); i < ins.I; i++)
							  x = grafo.addNode();
						Node supersup=grafo.addNode();
						
						double supersup_div(0.0);
						
						NodeDou divergence(grafo);
						for (t = 0; t < T; t++){
							for(i = 0; i < I ; i++){
									coe=ins.m.at(i).at(t).at(v);															
									for (j = 0; j < J; j++) if(t<T) coe -= (getSolution(X.at(i).at(j).at(t).at(v)) > EPS4?getSolution(X.at(i).at(j).at(t).at(v)):0.0) ;																											
									for (j = 0; j < J; j++) if((i != j) && ((t - ins.tau.at(j).at(i)) >= 0) ) coe += (getSolution(X.at(j).at(i).at(t - ins.tau.at(j).at(i)).at(v)) > EPS4?getSolution(X.at(j).at(i).at(t - ins.tau.at(j).at(i)).at(v)):0.0) ;
								
									//if(coe > EPS2 || coe < -EPS2 )
									//cout << "\ncoe[" << flat_ixt(ins,t,i) /*<< "]-[" << i << "][" << t*/ << "]=" << coe << flush;
									
									if(std::round(coe) < -EPS4) supersup_div += fabs(std::round(coe));
									divergence[grafo.nodeFromId(flat_ixt(ins,t,i))] = std::round(coe);
								
									if(fabs(divergence[grafo.nodeFromId(flat_ixt(ins,t,i))]-coe) > EPS3 ){
										cerr << "\nError in setting divergence in subproblem: " << v 
											 << " as LEMON( " << divergence[grafo.nodeFromId(flat_ixt(ins,t,i))] 
											 << " ) and COE( " << coe << " ) does not match " << flush;
										exit(EXIT_FAILURE);
									}
									
									//if(divergence[grafo.nodeFromId(flat_ixt(ins,t,i))] > EPS2 || divergence[grafo.nodeFromId(flat_ixt(ins,t,i))] < -EPS2 )
									//cout << "\tdiv[" << flat_ixt(ins,t,i) /*<< "]-[" << i << "][" << t*/ << "]=" << divergence[grafo.nodeFromId(flat_ixt(ins,t,i))] << flush;
											
							}
						}
						divergence[supersup] = supersup_div;
						//cout << "\nsupersup: " << supersup_div << endl;
																	
						ArcDou costo(grafo);
						Arc arco;
						int count( 0 );
						for(t = 0; t < ins.T; ++t){
							for(i = 0; i < ins.I; ++i){
								for(j = 0; j < ins.J; ++j){
									if( ins.A.at(i).at(j).at(v) == 1 ){
										//cout << "\ni: " << i << " j: " << j << " t: " << t << flush;
										if ( i == j ){
											//cout << "\t\t" << i + ins.I * t << "\t" << j + (ins.I * (t+1)) << flush;
											arco = grafo.addArc(grafo.nodeFromId(i + ins.I * t),grafo.nodeFromId(j + (ins.I * (t+1))));
											costo.set(arco, (ins.A.at(i).at(j).at(v) == 1?ins.cos.at(i).at(j).at(v):DBL_MAX) );
											
										}else if ( (t + ins.tau.at(i).at(j)) <= ins.T-1 ){  
											//cout << "\t\t" << i + ins.I * t << "\t" << j + ins.I * (t + ins.tau.at(i).at(j)) << flush;
											arco = grafo.addArc(grafo.nodeFromId(i + ins.I * t),grafo.nodeFromId(j + ins.I * (t + ins.tau.at(i).at(j))));
											costo.set(arco, (ins.A.at(i).at(j).at(v) == 1?ins.cos.at(i).at(j).at(v):DBL_MAX) );
										}
										else if ( (t + ins.tau.at(i).at(j)) > ins.T-1 ){
											//cout << "\t\t" << i + ins.I * t << "\t" << ((ins.I*(ins.T+1)+1))-1-(ins.I)+j << flush;  
											arco = grafo.addArc(grafo.nodeFromId(i + ins.I * t),grafo.nodeFromId(((ins.I*(ins.T+1)+1))-1-(ins.I)+j)) ;
											costo.set(arco, (ins.A.at(i).at(j).at(v) == 1?ins.cos.at(i).at(j).at(v):DBL_MAX) );
										}
										//cout << "\t\t" << grafo.id(arco) << flush;
									}
								}
							}
						}
						
						for(t = 0; t < ins.T; ++t){
							for(i = 0; i < ins.I; ++i){
								arco = grafo.addArc(supersup,grafo.nodeFromId(i + ins.I * t));
								costo.set(arco, ins.cosf.at(i).at(t).at(v) );
							}
						}
						
						
						/*Node head; 
						int id1,id2,ti,tj;
						for (ArcIt u(grafo); u != INVALID; ++u){
							
							id1= grafo.id(grafo.source(u)); 
							id2= grafo.id(grafo.target(u));
							
							if( id1 < (ins.I*(ins.T+1)) && id2 < (ins.I*(ins.T+1))){
								
								ixt_flat(ins,id1,ti,i );ixt_flat(ins,id2,tj,j );
								costo.set(u, (ins.A.at(i).at(j).at(v) == 1?ins.cos.at(i).at(j).at(v):DBL_MAX) );
								
								
							}else{
								
								cerr << "\nTHIS OPTION MUST NOT HAPPEN" << endl;exit(EXIT_FAILURE);
								
							}
							
							
						}*/
					
						//for (ArcIt u(grafo); u != INVALID; ++u){
							
							//cout << "\nArcId: " << setw(5) << grafo.id(u)  
												 //<< setw(5) << grafo.id(grafo.source(u))  
												 //<< setw(5) << grafo.id(grafo.target(u))  
												 //<< setw(5) << costo[u]  
												 //<< flush;
							
						//}
						
						
						
						
						
						NS sol_NS(grafo);
						sol_NS.costMap(costo).supplyMap(divergence).supplyType(NS::LEQ);
						NS::ProblemType stat_NS = sol_NS.run();
						
		

						/*
						cerr << "\nOPTIMAL COST=" << sol_NS.totalCost() << endl; 		
						for (ArcIt n(grafo); n != INVALID; ++n){
							if( flo[n] != 0 ){
									Node x = grafo.source(n);
									Node y = grafo.target(n);
									cout << "\nflow[" << grafo.id(n) << "]=" << flo[n] << flush;
									cout << "\t[" << grafo.id(x) << "][" << grafo.id(y) << "]=" << flo[n] << " Costo: " << costo[n] << flush;
			
							}
						}
						for (NodeIt n(grafo); n != INVALID; ++n){
							if( po[n] != 0  ){
								cout << "\npot[" << grafo.id(n) << "]="<< setprecision(30) << po[n] << flush;
							}
						}*/			
					
						
						if(stat_NS == NS::OPTIMAL ){
							
							ZLe = sol_NS.totalCost();
							ArcDou flo(grafo);
							NodeDou po(grafo);
							sol_NS.flowMap(flo);
							sol_NS.potentialMap(po);
							
							if( ZLe > getSolution(u.at(v)) ){
								
								opti++;
								GRBLinExpr lhs = 0;
								for(int i=0;i<I;i++){
									for(int t=0;t<T;t++){
										//cout << "\nOut" << flush;
										for (j = 0; j < J; j++){
											if( po[grafo.nodeFromId(flat_ixt(ins,t,i))] > EPS3 )	
												lhs -= po[grafo.nodeFromId(flat_ixt(ins,t,i))]*X.at(i).at(j).at(t).at(v);
												
												
										}
										//cout << "\nIn" << flush;
										for (j = 0; j < J; j++){
											if((i != j) && ((t - ins.tau.at(j).at(i)) >= 0) ){
												if( po[grafo.nodeFromId(flat_ixt(ins,t,i))] > EPS3 )
													lhs += po[grafo.nodeFromId(flat_ixt(ins,t,i))]*X.at(j).at(i).at((t - ins.tau.at(j).at(i))).at(v);
											}
										}
									}
								}
								lhs += u.at(v);
								
								int rhs = 0;
								for(int i=0; i < I ; i++){
									for(int t=0;t < T;t++){
										rhs -= (ins.m.at(i).at(t).at(v)*po[grafo.nodeFromId(flat_ixt(ins,t,i))]);
									}
								}
								
								addLazy(lhs >= rhs);
								//cout << "\n\tOptimality Cut Added" << flush;
								//cout << "\n" << lhs << ">=" << rhs << flush;
								
								
							}
							
						}else{
							cerr << "\n The subproblem " << v << " is not optimal, hence there is an error. The model is " << flush;
							switch (stat_NS) {
							case NS::INFEASIBLE:
								cerr << "INFEASIBLE." << endl;
							break;
							case NS::UNBOUNDED:
								cerr << "UNBOUNDED." << endl;
							break;
							default:
								break;
							}
							exit(EXIT_FAILURE);
						}
						
						grafo.clear();
						#endif
						
						
				}/*end_if v*/
				
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
        					#ifdef primals
        					/* ******************************************************************** */
        					GRBEnv erv = GRBEnv(true);
        					erv.set("LogFile", "VAP.log");
        					erv.start();
        					GRBModel modelSubP2 = GRBModel(erv);
        					modelSubP2.set(GRB_IntParam_OutputFlag, 0);
        					modelSubP2.set(GRB_StringAttr_ModelName, "Primal_VAP");
        					modelSubP2.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
        					vec1GRBvar      L(T);
        					for (t = 0; t < T; t++){
        						var_name = "I("+to_string(t)+")";
        						L.at(t) = modelSubP2.addVar(0.0, (T-1==0?0.0:GRB_INFINITY), ins.In.at(i).at(j).at(t), GRB_CONTINUOUS, var_name);
        					}
        					for (t = 0; t < T; t++){
        						var_name = "Dem("+to_string(t)+")";
        						coe = ins.dem.at(i).at(j).at(t);
        						for (v = 0; v < V; v++) coe -= (getSolution(X.at(i).at(j).at(t).at(v)) > EPS4?getSolution(X.at(i).at(j).at(t).at(v)):0.0);
        						
        						GRBLinExpr lhs = 0;
        						/*if(t<T-1)*/ lhs += L.at(t);
        						if(t>0)   lhs -= L.at(t-1);
        						modelSubP2.addConstr(lhs, GRB_EQUAL ,coe, var_name);
        					}
        					//cout << "\nOK2" << flush;
        					// Optimize model
        					modelSubP2.update();
        					sub = "subPII-"+to_string(i)+","+to_string(j)+".lp";
        					modelSubP2.write(sub);
        					modelSubP2.optimize();
        					//double ZIP = modelSubP2.get(GRB_DoubleAttr_ObjVal);
        					//cout << "\nOK3" << flush;
        					/* ******************************************************************** */
        					#endif
        					/* Modify objective function */
        					for (t = 0; t < T; t++){ 
        						coe=ins.dem.at(i).at(j).at(t);
        						for ( v=0; v<V; v++ ) coe -= (getSolution(X.at(i).at(j).at(t).at(v)) > EPS4?getSolution(X.at(i).at(j).at(t).at(v)):0.0);
        						w.at(t).set(GRB_DoubleAttr_Obj,std::round(coe));
        						//cout << "\n\tdiv[" << t << "]=" << setw(4) << std::round(coe) << flush;
        						
        					}
        					//cout << endl;
        					
        					modelSubD2.update();
        					
        					/* Subproblem constraints Definition */
        					for (t = 0; t < T-1; t++){ 
        						var_name = "In("+to_string(t)+")";
        						//modelSubD2.addConstr((t==T-1?-w[t]:-w[t]+w[t+1]), GRB_LESS_EQUAL , ins.In.at(i).at(j).at(t), var_name);
        						modelSubD2.addConstr(w[t]-w[t+1], GRB_LESS_EQUAL , ins.In.at(i).at(j).at(t), var_name);///////////cambio///////////
        					}
        					
        					modelSubD2.update();		
        						
        					modelSubD2.set(GRB_IntParam_Presolve, 0);
        					//sub = "subDII-"+to_string(i)+","+to_string(j)+".lp";
        					//modelSubD2.write(sub);
        					// Optimize model
        					modelSubD2.optimize();
        					
        					
        					
        					if(modelSubD2.get(GRB_IntAttr_Status) == 2 ){
        						
        						double ZD = modelSubD2.get(GRB_DoubleAttr_ObjVal);
        						
        						
        						#ifdef primals
        						double ZP = modelSubP2.get(GRB_DoubleAttr_ObjVal);
        						if( fabs(ZP-ZD) > EPS3 ){
        							cerr << "\n\tThe optimal solution of the primal(" << ZP << ") does not match the optimal solution of the dual(" << ZD << ")." << flush;
        							exit(EXIT_FAILURE);
        						}
        						#endif
        						
        						
        						if( ZD > getSolution(uu.at(i).at(j)) ){
        							
        							//cout << "\n\t    t    " << flush;
        							//for(int t=0;t<T;t++)
        								////if( w[t].get(GRB_DoubleAttr_X) > EPS4 || w[t].get(GRB_DoubleAttr_X) < -EPS4)
        									//cout << "\n\tsol[" << t << "] = " << w[t].get(GRB_DoubleAttr_X) << flush;
        								
        							
        							opti_2++;
        							GRBLinExpr lhs = 0;
        							for(int t=0;t<T;t++){
        								for ( Veh v(0); v<V; v++ ){
        									if( w[t].get(GRB_DoubleAttr_X) > EPS4 || w[t].get(GRB_DoubleAttr_X) < -EPS4){
        										lhs += w[t].get(GRB_DoubleAttr_X)*X.at(i).at(j).at(t).at(v);
        									}
        								}
        							}
        							
        							lhs += uu.at(i).at(j);
        						
        							double rhs = 0;
        							for(int t=0;t < T;t++){
        								if( w[t].get(GRB_DoubleAttr_X) > EPS4 || w[t].get(GRB_DoubleAttr_X) < -EPS4)
        									rhs += (double(ins.dem.at(i).at(j).at(t))*w[t].get(GRB_DoubleAttr_X));
        							}
        						
        							
        							//cout << "\n\n\t" << lhs << " <= " << rhs << flush;
        							//cuts << "\n\n\t" << lhs << " <= " << rhs << flush;
        							addLazy(lhs >= rhs);
        							//cout << "\n\tOptimality Cut Added Inventory" << flush;
        							
        						}
        					}
        					
        					if(modelSubD2.get(GRB_IntAttr_Status) == 5 ){
        						
        						//cout << "\n\t    t    " << flush;
        							//for(int t=0;t<T;t++)
        								//if( w[t].get(GRB_DoubleAttr_UnbdRay) > EPS3 || w[t].get(GRB_DoubleAttr_UnbdRay) < -EPS3 ){
        									//cout << "\n\tray[" << t << "] = " << w[t].get(GRB_DoubleAttr_UnbdRay) << flush;
        									////if( w[t].get(GRB_DoubleAttr_UnbdRay) > 1  ) while(getchar() != '\n');
        								//}
        						
        						feas_2++;
        						GRBLinExpr lhs = 0;
        						for(int t=0;t<T;t++){
        							for ( Veh v(0); v<V; v++ ){
        								if( w[t].get(GRB_DoubleAttr_UnbdRay) > EPS3 || w[t].get(GRB_DoubleAttr_UnbdRay) < -EPS3 ){
        									lhs += w[t].get(GRB_DoubleAttr_UnbdRay)*X.at(i).at(j).at(t).at(v);
        								}
        							}
        						}
        						
        						//lhs += uu.at(i).at(j);
        					
        						double rhs = 0;
        						for(int t=0;t < T;t++){
        							if( w[t].get(GRB_DoubleAttr_UnbdRay) > EPS3 || w[t].get(GRB_DoubleAttr_UnbdRay) < -EPS3 )
        								rhs += (double(ins.dem.at(i).at(j).at(t))*w[t].get(GRB_DoubleAttr_UnbdRay));
        						}
        					
        						//cout << "\n\n\t" << lhs << " <= " << rhs << flush;
        						//cuts << "\n\n\t" << lhs << " <= " << rhs << flush;
        						addLazy(lhs >= rhs);
        						//cout << "\n\tFeasibility Cut Added Inventory" << flush;
        						
        
        						
        					}
        					
        				}
        				
        				/* Remove constraints to process next vehicle*/
        				GRBConstr* array_ctr = 	modelSubD2.getConstrs ( );
        				for (int k = 0; k < modelSubD2.get(GRB_IntAttr_NumConstrs); ++k)modelSubD2.remove(array_ctr[k]);
        				modelSubD2.update();
        				
        				
        			}/*end_if J*/
        		}/*end_if I*/
			  
			   
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

