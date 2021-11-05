#include "BFVAP_callback.cpp"

#ifdef comparator
double BFTVAP_solveOriginal( Instance_data &ins ) 	
#else
int BFTVAP_solveOriginal( Instance_data &ins ) 	
#endif
{         
	
    try {  
		cout << "\nHello from BFTVAP_solveOriginal()\n" << endl;
		
		vec3Dou &pro     = ins.pro;
		vec3Dou &cos     = ins.cos;
		vec3Dou &cosf    = ins.cosf;
		vec2Int &tau     = ins.tau;
		vec3Int &m       = ins.m;
		vec3Int &dem     = ins.dem;
		vec3Int &A       = ins.A;
		vec3Dou &In      = ins.In;
		vec2Int &Kap     = ins.Kap;
		
		int I = ins.I , J = ins.I , T = ins.T , V = ins.V,i,j,k,t,v;
		
		string var_name;  
		var_name = "BFVAP_"+ins.name +".log";

		// Create an environment
		GRBEnv env = GRBEnv(true);
		env.set("LogFile", var_name);
		env.start();
		
		var_name.clear();
	
		
		// Create an empty model
		GRBModel model = GRBModel(env);
		model.set(GRB_StringAttr_ModelName, "BFTVAP");
		model.set(GRB_IntParam_Threads, 1);
		model.set(GRB_DoubleParam_MIPGap, 0.0);
		
		// Turn off display and heuristics
		//model.set(GRB_IntParam_OutputFlag, 0);
		//model.set(GRB_DoubleParam_Heuristics, 0.0);
		
		
		model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
		
		
		
		vec4GRBVar X(I,vec3GRBVar(J,vec2GRBVar(T,vec1GRBvar(V))));
		vec4GRBVar Y(I,vec3GRBVar(J,vec2GRBVar(T,vec1GRBvar(V))));
		vec3GRBVar Z(I,vec2GRBVar(T,vec1GRBvar(V)));
		vec3GRBVar Q(I,vec2GRBVar(J,vec1GRBvar(T)));
		cout << "\nVariable Declaration is Done " << flush;
		for(i = 0; i < I ; i++)
			for(j = 0; j < J; j++)
				for (t = 0; t < T; t++){
					var_name="q("+to_string(i)+","+to_string(j)+","+to_string(t)+")";
					//Q.at(i).at(j).at(t) = model.addVar((t==T-1?-GRB_INFINITY:0.0), (t==T-1?0.0:GRB_INFINITY), -In.at(i).at(j).at(t), GRB_CONTINUOUS, var_name);
					Q.at(i).at(j).at(t) = model.addVar(0.0, (t==T-1?0.0:GRB_INFINITY), -In.at(i).at(j).at(t), GRB_CONTINUOUS, var_name);
					for (v = 0; v < V; v++){
						var_name="x("+to_string(i)+","+to_string(j)+","+to_string(t)+","+to_string(v)+")";
						X.at(i).at(j).at(t).at(v) = model.addVar(0.0, ((A.at(i).at(j).at(v)==0)||(i==j)?0.0:GRB_INFINITY), pro.at(i).at(j).at(v), GRB_CONTINUOUS, var_name);
						var_name="y("+to_string(i)+","+to_string(j)+","+to_string(t)+","+to_string(v)+")";
						Y.at(i).at(j).at(t).at(v) = model.addVar(0.0, (A.at(i).at(j).at(v)==0?0.0:GRB_INFINITY), -cos.at(i).at(j).at(v), GRB_CONTINUOUS, var_name);
						
					}
				}
		for(i = 0; i < I ; i++)
			for (t = 0; t < T; t++)
				for (v = 0; v < V; v++){
					var_name="z("+to_string(i)+","+to_string(t)+","+to_string(v)+")";
					Z.at(i).at(t).at(v) = model.addVar(0.0, GRB_INFINITY, -cosf.at(i).at(t).at(v), GRB_CONTINUOUS, var_name);
				}
		
		
		cout << "\nVariable Definition is Done " << flush;
		
		
		
		for(i = 0; i < I ; i++)
			for (t = 0; t < T; t++)
				for (v = 0; v < V; v++){
					var_name="flux("+to_string(i)+","+to_string(t)+","+to_string(v)+")";
					GRBLinExpr sum_expr = 0;
					for (j = 0; j < J; j++) sum_expr += (X[i][j][t][v] + Y[i][j][t][v]);
					for (k = 0; k < J; k++) if((i != k) && ((t - tau.at(k).at(i)) >= 0) ) sum_expr -= (X[k][i][t - tau.at(k).at(i)][v] + Y[k][i][t - tau.at(k).at(i)][v]);
					if (t > 0) sum_expr -= Y[i][i][t-1][v];
					sum_expr -= Z[i][t][v];
					//model.addConstr(sum_expr, GRB_EQUAL ,m.at(i).at(t).at(v), var_name);
					model.addConstr(sum_expr, GRB_LESS_EQUAL ,m.at(i).at(t).at(v), var_name);
					
				}
		
		cout << "\nFlow Constraint Definition is Done " << flush;
		
		for (i = 0; i < I; i++)

			for (t = 0; t < T; t++){
				var_name="cap("+to_string(i)+","+to_string(t)+")";
				GRBLinExpr sum_expr = 0;
				for (v = 0; v < V; v++)
					for (j = 0; j < J; j++)
						 if((i != j) && ((t - tau.at(j).at(i)) >= 0) )
							 sum_expr += X[j][i][t - tau.at(j).at(i)][v] ;
				model.addConstr(sum_expr, GRB_LESS_EQUAL ,Kap.at(i).at(t), var_name);
			}
		
		cout << "\nKapacity Constraint Definition is Done " << flush;
		
				
		for(i = 0; i < I ; i++)
			for(j = 0; j < J; j++)
				for (t = 0; t < T; t++){
					    var_name="backl("+to_string(i)+","+to_string(j)+","+to_string(t)+")";
					    GRBLinExpr sum_expr = 0;	 
						sum_expr += Q[i][j][t];
						if(t > 0) sum_expr -= Q[i][j][t-1];
						for(v = 0; v < V; v++) sum_expr += X[i][j][t][v];
						model.addConstr(sum_expr, GRB_EQUAL ,dem.at(i).at(j).at(t), var_name);
						//model.addConstr(sum_expr, GRB_LESS_EQUAL ,dem.at(i).at(j).at(t), var_name);
						//model.addConstr(sum_expr, GRB_GREATER_EQUAL ,dem.at(i).at(j).at(t), var_name);
					}
		
		cout << "\nBacklog COnstraints is done " << flush;

					
			 
		  
		//model.write("BFTVAP.lp");
				
		high_resolution_clock::time_point tbi = high_resolution_clock::now();
		// Optimize model
		model.optimize();
		#ifdef comparator
		if( model.get(GRB_IntAttr_Status) == 3 || model.get(GRB_IntAttr_Status) == 4  ){
			cerr << "\nSolution is infeasible or unbounded " << flush;
			exit(EXIT_SUCCESS);
		}
		#endif
		if( model.get(GRB_IntAttr_Status) == 4 ){
			model.set(GRB_IntParam_DualReductions, 0);
			model.set(GRB_IntParam_Presolve, 0);
			model.update();
			model.optimize();
			print_status( model );
		}
		
		high_resolution_clock::time_point tbf = high_resolution_clock::now();
		seconds_type solver_runtime = duration_cast<seconds_type>(tbf - tbi);
		if( model.get(GRB_IntAttr_Status) == 5 ){
			for (v = 0; v < V; v++){
				cout << endl;
				for (t = 0; t < T; t++){
					for(i = 0; i < I ; i++){
						
						if( Z.at(i).at(t).at(v).get(GRB_DoubleAttr_UnbdRay) > EPS4 ){
							cout << endl << Z.at(i).at(t).at(v).get(GRB_StringAttr_VarName) << " = " << Z.at(i).at(t).at(v).get(GRB_DoubleAttr_UnbdRay) << setw(10) << flat_ixt(ins,t,i) << flush;
							////cout << setw(15) << ins.cosf.at(i).at(t).at(v) << flush;
							////cout << setw(15) << ins.cosf.at(i).at(t).at(v)*( Z.at(i).at(t).at(v).get(GRB_DoubleAttr_UnbdRay) - ins.m.at(i).at(t).at(v)) << flush;
							
						}
						
						for(j = 0; j < J; j++){	
							if( X.at(i).at(j).at(t).at(v).get(GRB_DoubleAttr_UnbdRay) > EPS4 ){
								cout << endl << X.at(i).at(j).at(t).at(v).get(GRB_StringAttr_VarName) << " = " << X.at(i).at(j).at(t).at(v).get(GRB_DoubleAttr_UnbdRay) << flush;
								cout << setw(15) << flat_ixt(ins,t,i) << setw(10) << flat_ixt(ins,(i==j?t+1:t + ins.tau.at(i).at(j)),j) << flush;
								//cout << setw(15) << pro.at(i).at(j).at(v) << flush;
								//cout << setw(15) << pro.at(i).at(j).at(v)*X.at(i).at(j).at(t).at(v).get(GRB_DoubleAttr_UnbdRay) << flush;
								
								
							}
							if( Y.at(i).at(j).at(t).at(v).get(GRB_DoubleAttr_UnbdRay) > EPS4 ){
								cout << endl << Y.at(i).at(j).at(t).at(v).get(GRB_StringAttr_VarName) << " = " << Y.at(i).at(j).at(t).at(v).get(GRB_DoubleAttr_UnbdRay) << flush;
								cout << setw(15) << flat_ixt(ins,t,i) << setw(10) << flat_ixt(ins,(i==j?t+1:t + ins.tau.at(i).at(j)),j) << flush;
								//cout << setw(15) << -cos.at(i).at(j).at(v) << flush;
								//cout << setw(15) << -cos.at(i).at(j).at(v)*Y.at(i).at(j).at(t).at(v).get(GRB_DoubleAttr_UnbdRay) << flush;
								
								
							}
							
						}
					}
				}
				
				
				
			}
			cout << endl;
			for (t = 0; t < T; t++){
				for(i = 0; i < I ; i++){
					for(j = 0; j < J; j++){	
						if( Q.at(i).at(j).at(t).get(GRB_DoubleAttr_UnbdRay) > EPS4 ){
							cout << endl << Q.at(i).at(j).at(t).get(GRB_StringAttr_VarName) << " = " << Q.at(i).at(j).at(t).get(GRB_DoubleAttr_UnbdRay) <<  flush;
							cout << setw(15) << flat_ixt(ins,t,i) << setw(10) << flat_ixt(ins,(i==j?t+1:t + ins.tau.at(i).at(j)),j) << flush;
							//cout << setw(15) << ins.In.at(i).at(j).at(t) << flush;
							//cout << setw(15) << ins.In.at(i).at(j).at(t)*Q.at(i).at(j).at(t).get(GRB_DoubleAttr_UnbdRay)  << flush;
							
						}
					}
				}
			}
								
						
		}
		if( model.get(GRB_IntAttr_Status) == 2 ){
				double ZLP = model.get(GRB_DoubleAttr_ObjVal);
				op_solV_LP << model.get(GRB_IntAttr_Status) << setw(10) << ZLP << setw(10) << solver_runtime.count() << flush;
				op_sol_LP  << model.get(GRB_IntAttr_Status) << setw(10) << ZLP << setw(10) << solver_runtime.count() << flush;
				
				for (v = 0; v < V; v++){
					op_sol_LP << endl;
					for (t = 0; t < T; t++){
						for(i = 0; i < I ; i++){
							
							if( Z.at(i).at(t).at(v).get(GRB_DoubleAttr_X) > EPS4 ){
								op_sol_LP << endl << Z.at(i).at(t).at(v).get(GRB_StringAttr_VarName) << " = " << Z.at(i).at(t).at(v).get(GRB_DoubleAttr_X) << setw(10) << flat_ixt(ins,t,i) << flush;
								op_sol_LP << setw(15) << ins.cosf.at(i).at(t).at(v) << flush;
								op_sol_LP << setw(15) << ins.cosf.at(i).at(t).at(v)*( Z.at(i).at(t).at(v).get(GRB_DoubleAttr_X) - ins.m.at(i).at(t).at(v)) << flush;
								
							}
							
							for(j = 0; j < J; j++){	
								if( X.at(i).at(j).at(t).at(v).get(GRB_DoubleAttr_X) > EPS4 ){
									op_sol_LP << endl << X.at(i).at(j).at(t).at(v).get(GRB_StringAttr_VarName) << " = " << X.at(i).at(j).at(t).at(v).get(GRB_DoubleAttr_X) << flush;
									op_sol_LP << setw(15) << flat_ixt(ins,t,i) << setw(10) << flat_ixt(ins,(i==j?t+1:t + ins.tau.at(i).at(j)),j) << flush;
									op_sol_LP << setw(15) << pro.at(i).at(j).at(v) << flush;
									op_sol_LP << setw(15) << pro.at(i).at(j).at(v)*X.at(i).at(j).at(t).at(v).get(GRB_DoubleAttr_X) << flush;
									
									
								}
								if( Y.at(i).at(j).at(t).at(v).get(GRB_DoubleAttr_X) > EPS4 ){
									op_sol_LP << endl << Y.at(i).at(j).at(t).at(v).get(GRB_StringAttr_VarName) << " = " << Y.at(i).at(j).at(t).at(v).get(GRB_DoubleAttr_X) << flush;
									op_sol_LP << setw(15) << flat_ixt(ins,t,i) << setw(10) << flat_ixt(ins,(i==j?t+1:t + ins.tau.at(i).at(j)),j) << flush;
									op_sol_LP << setw(15) << -cos.at(i).at(j).at(v) << flush;
									op_sol_LP << setw(15) << -cos.at(i).at(j).at(v)*Y.at(i).at(j).at(t).at(v).get(GRB_DoubleAttr_X) << flush;
									
									
								}
								
							}
						}
					}
					
					
					
				}
				op_sol_LP << endl;
				for (t = 0; t < T; t++){
					for(i = 0; i < I ; i++){
						for(j = 0; j < J; j++){	
							if( Q.at(i).at(j).at(t).get(GRB_DoubleAttr_X) > EPS4 ){
								op_sol_LP << endl << Q.at(i).at(j).at(t).get(GRB_StringAttr_VarName) << " = " << Q.at(i).at(j).at(t).get(GRB_DoubleAttr_X) <<  flush;
								op_sol_LP << setw(15) << flat_ixt(ins,t,i) << setw(10) << flat_ixt(ins,(i==j?t+1:t + ins.tau.at(i).at(j)),j) << flush;
								op_sol_LP << setw(15) << ins.In.at(i).at(j).at(t) << flush;
								op_sol_LP << setw(15) << ins.In.at(i).at(j).at(t)*Q.at(i).at(j).at(t).get(GRB_DoubleAttr_X)  << flush;
								
							}
						}
					}
				}
								
							

				
				#ifndef comparator
				resultados_OV << "For_LPBF:"  << setprecision(6) << setw(15) << ZLP << "\t" << fixed << showpoint << setprecision(6) << solver_runtime.count() <<  noshowpoint << setprecision(0) <<  flush;
				#endif
				
				
				/* **************** Solve Integer Problem*********************** */
				
				for(i = 0; i < I ; i++)
					for(j = 0; j < J; j++)
						for (t = 0; t < T; t++){
							Q.at(i).at(j).at(t).set(GRB_CharAttr_VType, GRB_INTEGER);
							for (v = 0; v < V; v++){
								X.at(i).at(j).at(t).at(v).set(GRB_CharAttr_VType, GRB_INTEGER);
								Y.at(i).at(j).at(t).at(v).set(GRB_CharAttr_VType, GRB_INTEGER);
							}
						}
							
				for(i = 0; i < I ; i++)
					for (t = 0; t < T; t++)
						for (v = 0; v < V; v++){
							Z.at(i).at(t).at(v).set(GRB_CharAttr_VType, GRB_INTEGER);
						}
				
				double start_cplex2 = clock();
				// Optimize model
				model.optimize();
				double ZIP = model.get(GRB_DoubleAttr_ObjVal);
				double stop_cplex2 = clock();
				double solver_runtime2 = (stop_cplex2-start_cplex2)/double(CLOCKS_PER_SEC);
				
				op_solV_IP << model.get(GRB_IntAttr_Status) << setw(10) << ZIP << setw(10) << solver_runtime.count() << flush;
				op_sol_IP  << model.get(GRB_IntAttr_Status) << setw(10) << ZIP << setw(10) << solver_runtime.count() << flush;
				
				for (v = 0; v < V; v++){
					op_sol_IP << endl;
					for (t = 0; t < T; t++){
						for(i = 0; i < I ; i++){
							
							if( Z.at(i).at(t).at(v).get(GRB_DoubleAttr_X) > EPS4 ){
								op_sol_IP << endl << Z.at(i).at(t).at(v).get(GRB_StringAttr_VarName) << " = " << Z.at(i).at(t).at(v).get(GRB_DoubleAttr_X) << setw(10) << flat_ixt(ins,t,i) << flush;
								op_sol_IP << setw(15) << ins.cosf.at(i).at(t).at(v) << flush;
								op_sol_IP << setw(15) << ins.cosf.at(i).at(t).at(v)*( Z.at(i).at(t).at(v).get(GRB_DoubleAttr_X) - ins.m.at(i).at(t).at(v)) << flush;
								
							}
							
							for(j = 0; j < J; j++){	
								if( X.at(i).at(j).at(t).at(v).get(GRB_DoubleAttr_X) > EPS4 ){
									op_sol_IP << endl << X.at(i).at(j).at(t).at(v).get(GRB_StringAttr_VarName) << " = " << X.at(i).at(j).at(t).at(v).get(GRB_DoubleAttr_X) << flush;
									op_sol_IP << setw(15) << flat_ixt(ins,t,i) << setw(10) << flat_ixt(ins,(i==j?t+1:t + ins.tau.at(i).at(j)),j) << flush;
									op_sol_IP << setw(15) << pro.at(i).at(j).at(v) << flush;
									op_sol_IP << setw(15) << pro.at(i).at(j).at(v)*X.at(i).at(j).at(t).at(v).get(GRB_DoubleAttr_X) << flush;
									
									
								}
								if( Y.at(i).at(j).at(t).at(v).get(GRB_DoubleAttr_X) > EPS4 ){
									op_sol_IP << endl << Y.at(i).at(j).at(t).at(v).get(GRB_StringAttr_VarName) << " = " << Y.at(i).at(j).at(t).at(v).get(GRB_DoubleAttr_X) << flush;
									op_sol_IP << setw(15) << flat_ixt(ins,t,i) << setw(10) << flat_ixt(ins,(i==j?t+1:t + ins.tau.at(i).at(j)),j) << flush;
									op_sol_IP << setw(15) << -cos.at(i).at(j).at(v) << flush;
									op_sol_IP << setw(15) << -cos.at(i).at(j).at(v)*Y.at(i).at(j).at(t).at(v).get(GRB_DoubleAttr_X) << flush;
									
									
								}
								
							}
						}
					}
					
					
					
				}
				op_sol_IP << endl;
				for (t = 0; t < T; t++){
					for(i = 0; i < I ; i++){
						for(j = 0; j < J; j++){	
							if( Q.at(i).at(j).at(t).get(GRB_DoubleAttr_X) > EPS4 ){
								op_sol_IP << endl << Q.at(i).at(j).at(t).get(GRB_StringAttr_VarName) << " = " << Q.at(i).at(j).at(t).get(GRB_DoubleAttr_X) <<  flush;
								op_sol_IP << setw(15) << flat_ixt(ins,t,i) << setw(10) << flat_ixt(ins,(i==j?t+1:t + ins.tau.at(i).at(j)),j) << flush;
								op_sol_IP << setw(15) << ins.In.at(i).at(j).at(t) << flush;
								op_sol_IP << setw(15) << ins.In.at(i).at(j).at(t)*Q.at(i).at(j).at(t).get(GRB_DoubleAttr_X)  << flush;
								
							}
						}
					}
				}
							
			
				#ifndef comparator
				resultados_OV << "\tFor_IPBF:" << setw(15) << setprecision(6) << ZIP << "\t" << fixed << showpoint << setprecision(6) << solver_runtime2 <<  noshowpoint << setprecision(0) <<  flush;
				#endif
				
				
				cout << "\nZIP  = " << ZIP  << endl;
				cout << "ZLP  = " << ZLP  << endl;
				double gap;
				gap = ((ZLP - ZIP)/ZIP)*100;
				cout << "GAP  = " << gap << " %" << endl;
				
				#ifndef comparator
				resultados_OV << "\t" << setw(10) << setprecision(6) << gap << flush;
				#endif
				// Solve the model with different values of Method
				//tune_parameters_VAP(model);
				
				cout << "\nBye from BFTVAP_solveOriginal()\n" << endl;
				
				#ifdef comparator
					resultados_OV << "OrgBF:" << "\t" << setw(8) << ZIP << "\t" << setw(10) << solver_runtime2 << flush;
					return ZIP;
				#endif
		}
		return 0;
		
		
		
  } catch(GRBException e) {
    cout << "\n\nError code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
    exit(EXIT_FAILURE);
  } catch(...) {
    cout << "\n\nException during optimization" << endl;
    exit(EXIT_FAILURE);
  }
  
  
	
}

#ifdef comparator
double BFTVAP_solveBenders( Instance_data &ins ) 	
#else
int BFTVAP_solveBenders( Instance_data &ins ) 	
#endif
{   
       
	
    try {  
		cout << "\nHello from BFTVAP_solveBenders()\n" << endl;
		
		vec3Dou &pro     = ins.pro;
		vec3Dou &cos     = ins.cos;
		vec3Dou &cosf    = ins.cosf;
		vec2Int &tau     = ins.tau;
		vec3Int &m       = ins.m;
		vec3Int &dem     = ins.dem;
		vec3Int &A       = ins.A;
		vec3Dou &In      = ins.In;
		vec2Int &Kap     = ins.Kap;
		
		int I = ins.I , J = ins.I , T = ins.T , V = ins.V,i,j,k,t,v;
		
		string var_name;  
		
		string sol_name = "solution"+ins.name+".dat";
		ofstream sol( sol_name, ios::out);
		
		int total_dem(0);
		for (t = 0; t < T; t++)
			for(i = 0; i < I ; i++)
				for(j = 0; j < J; j++)
					if( ins.dem.at(i).at(j).at(t) > 0 )
						total_dem += ins.dem.at(i).at(j).at(t);

		var_name = "BFVAP_Bender_"+ins.name +".log";
		
		// Create an environment
		GRBEnv env = GRBEnv(true);
		env.set("LogFile", var_name);
		env.start();
		var_name.clear();
		
		// Create an empty model
		GRBModel model = GRBModel(env);
		model.set(GRB_StringAttr_ModelName, "BFTVAP_Benders");
		
		// Turn off display and heuristics
		//model.set(GRB_IntParam_OutputFlag, 0);
		//model.set(GRB_DoubleParam_Heuristics, 0.0);
		model.set(GRB_IntParam_Threads, 1);
		model.set(GRB_DoubleParam_MIPGap, 0.0);
		model.set(GRB_DoubleParam_TimeLimit, 1000);
		
		model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
		
		// Must set LazyConstraints parameter when using lazy constraints
		model.set(GRB_IntParam_LazyConstraints, 1);
		
		/* ----------------------Subproblem Structure -------------------------*/
		
		GRBModel modelSubD1 = GRBModel(env);
		modelSubD1.set(GRB_StringAttr_ModelName, "BFTVAP_sub1");
		modelSubD1.set(GRB_IntParam_OutputFlag, 0);
		modelSubD1.set(GRB_StringAttr_ModelName, "D_MCFP1");
		modelSubD1.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
		vec2GRBVar      p(I,vec1GRBvar(T));
		for (t = 0; t < T; t++)
			for(i = 0; i < I ; i++){
				var_name = "p("+to_string(i)+","+to_string(t)+")";
				p.at(i).at(t) = modelSubD1.addVar(-GRB_INFINITY,0.0, 0.0, GRB_CONTINUOUS, var_name);
			}
			
			
		GRBModel modelSubD2 = GRBModel(env);
		modelSubD2.set(GRB_StringAttr_ModelName, "BFTVAP_sub2");
		modelSubD2.set(GRB_IntParam_OutputFlag, 0);
		modelSubD2.set(GRB_StringAttr_ModelName, "D_MCFP2");
		modelSubD2.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);///////////cambio///////////
		vec1GRBvar      w(T);
		for (t = 0; t < T; t++){
			var_name = "w("+to_string(t)+")";
			w.at(t) = modelSubD2.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, var_name);
		}
		
	
	
		
		
		/* ----------------------Subproblem Structure -------------------------*/
		
		cout << "\nSubproblem structure has been set " << endl;
		vec4GRBVar X(I,vec3GRBVar(J,vec2GRBVar(T,vec1GRBvar(V))));
		vec1GRBvar u(V);
		vec2GRBVar uu(I,vec1GRBvar(J));
		
		
		
		
		for(i = 0; i < I ; i++)
			for(j = 0; j < J; j++)
				for (t = 0; t < T; t++)
					for (v = 0; v < V; v++){
						var_name = "x("+to_string(i)+","+to_string(j)+","+to_string(t)+","+to_string(v)+")";
						X.at(i).at(j).at(t).at(v) = model.addVar(0.0, ((A.at(i).at(j).at(v)==0)||(i==j)?0.0:(t+tau.at(i).at(j)>T-1?total_dem:GRB_INFINITY)), pro.at(i).at(j).at(v), GRB_INTEGER, var_name);
					}
		
		
		cout << "\nVariables X has been set " << endl;			
		for (v = 0; v < V; v++){
			var_name = "u("+to_string(v)+")";
			u.at(v) = model.addVar(0.0, GRB_INFINITY, -1.0, GRB_INTEGER, var_name);
		}
		cout << "\nVariables U has been set " << endl;	
		
	
		for(i = 0; i < I ; i++)
			for(j = 0; j < J; j++){
				var_name = "uu("+to_string(i)+","+to_string(j)+")";
				uu.at(i).at(j) = model.addVar(0.0, GRB_INFINITY, -1.0, GRB_CONTINUOUS, var_name);///////////cambio///////////
		    }
		    
		cout << "\nVariables UU has been set " << endl;	
		
		
		for (i = 0; i < I; i++)
			for (t = 0; t < T; t++){
				var_name="cap("+to_string(i)+","+to_string(t)+")";
				GRBLinExpr sum_expr = 0;
				for (v = 0; v < V; v++)
					for (j = 0; j < J; j++)
						 if((i != j) && ((t - tau.at(j).at(i)) >= 0) )
							 sum_expr += X[j][i][t - tau.at(j).at(i)][v] ;
				model.addConstr(sum_expr, GRB_LESS_EQUAL ,Kap.at(i).at(t), var_name);
			}
		cout << "\nTerminal Capacity has been set " << endl;	
		model.update();
		
		cout << "\nSize of constraints: " << model.get(GRB_IntAttr_NumConstrs) << flush;
		//model.write("Benders_BFTVAP.lp");
		
		
		BFTBendersSub cb = BFTBendersSub(ins,X,u,uu,p,w,modelSubD1,modelSubD2);
		
		model.setCallback(&cb);
		
		vec4Dou X_sol(I,vec3Dou(J,vec2Dou(T,vec1Dou(V,0.0))));   
        
        
        
        
        high_resolution_clock::time_point tbi = high_resolution_clock::now();
		// Optimize model
		model.optimize();
		double ZIP = model.get(GRB_DoubleAttr_ObjVal);
		high_resolution_clock::time_point tbf = high_resolution_clock::now();
		seconds_type solver_runtime = duration_cast<seconds_type>(tbf - tbi);
		
		
		cout << (model.get(GRB_IntAttr_Status)==2?"\nOptimal Solution":(model.get(GRB_IntAttr_Status)==5?"\nUnbounded Solution":"\nNeither")) << flush; 
		/*
		for (v = 0; v < V; v++){
			cout << endl;
			double obj_veh( 0.0 );
			for (t = 0; t < T; t++){
				for(i = 0; i < I ; i++){
					for(j = 0; j < J; j++){	
						if( X.at(i).at(j).at(t).at(v).get(GRB_DoubleAttr_X) > EPS4 ){
							X_sol.at(i).at(j).at(t).at(v) = X.at(i).at(j).at(t).at(v).get(GRB_DoubleAttr_X);
							sol << "\n" << setw(5) << i << setw(5) << j << setw(5) << t << setw(5) << v << setw(5) << X.at(i).at(j).at(t).at(v).get(GRB_DoubleAttr_X) << flush;
							cout << "\n" << X.at(i).at(j).at(t).at(v).get(GRB_StringAttr_VarName) << " = " << X.at(i).at(j).at(t).at(v).get(GRB_DoubleAttr_X) << flush;
							cout << "\t" << pro.at(i).at(j).at(v) << flush;
							cout << setw(10) << flat_ixt(ins,t,i) << setw(10) << flat_ixt(ins,(i==j?t+1:t + ins.tau.at(i).at(j)),j) << flush;
							obj_veh += X.at(i).at(j).at(t).at(v).get(GRB_DoubleAttr_X)*pro.at(i).at(j).at(v);
						}
					}
				}
			}
			cout << "\nvehicle_cost[" << v << "]=" << obj_veh << flush;
		}
		//cout << "\nHALO1" << flush;
		for (v = 0; v < V; v++){
			if( u.at(v).get(GRB_DoubleAttr_X) > EPS4 ){
				cout << "\n" << u.at(v).get(GRB_StringAttr_VarName) << " = " << u.at(v).get(GRB_DoubleAttr_X) << flush;
			}
		}
		//cout << "\nHALO2" << flush;
		#ifndef second_benders_BFVAP
		for(i = 0; i < I ; i++){
			for(j = 0; j < J; j++){
				if( uu.at(i).at(j).get(GRB_DoubleAttr_X) < -EPS4 ){
					cout << "\n" << uu.at(i).at(j).get(GRB_StringAttr_VarName) << " = " << uu.at(i).at(j).get(GRB_DoubleAttr_X) << flush;
				}
			}
		}
		#endif
		//cout << "\nHALO2" << flush;
		*/
		
		cout << "\nZIP  = " << ZIP  << endl;
		cout << "\nN_opti: " << opti << " N_feas: " << feas << flush;
		cout << "\tN_opti_2: " << opti_2 << " N_feas_2: " << feas_2 << flush;
		cout << "\tsub_rt_both: " << sub_runtime_both << " sub_rt_1st: " << sub_runtime_1st<< " sub_rt_2nd: " << sub_runtime_2nd << flush;
		//cout << "\nsub_runtime_div: " << sub_runtime_div << " sub_runtime_cos: " << sub_runtime_cos<< " sub_runtime_run: " << sub_runtime_run << " sub_runtime_cut: " << sub_runtime_cut << flush;
		//cout << "\nsub_runtime_1: " << sub_runtime_1 << " sub_runtime_2: " << sub_runtime_2<< " sub_runtime_3: " << sub_runtime_3 << " sub_runtime_4: " << sub_runtime_4 << flush;
		//model.write("Final_Benders.lp");
		
		resultados_OV << "\t" << "BenFleet:" << "\t" << setprecision(2) << setw(8) << ZIP  << setw(10) << solver_runtime.count() << flush;
		resultados_OV << "\tN_opti: "   << setw(8) << opti   << " N_feas: "   << setw(8) << feas << flush;
		resultados_OV << "\tN_opti_2: " << setw(8) << opti_2 << " N_feas_2: " << setw(8) << feas_2 << flush;
		resultados_OV << setprecision(4);
		resultados_OV << "\tsub_rt_both: "<< setw(8) << sub_runtime_both << " sub_rt_1st: " << setw(8) << sub_runtime_1st<< " sub_rt_2nd: " << setw(8) << sub_runtime_2nd << flush;
		//resultados_OV << "\tsub_runtime_div: " << setw(8) << sub_runtime_div << " sub_runtime_cos: " << setw(8) << sub_runtime_cos<< " sub_runtime_run: " << setw(8) << sub_runtime_run << " sub_runtime_cut: " << setw(8) << sub_runtime_cut << flush;
		//resultados_OV << "\tsub_runtime_1: " << setw(8) << sub_runtime_1 << " sub_runtime_2: " << setw(8) << sub_runtime_2<< " sub_runtime_3: " << setw(8) << sub_runtime_3 << " sub_runtime_4: " << setw(8) << sub_runtime_4 << flush;
		
		/*{
			cout << endl << flush;
			
			GRBModel modelSubY = GRBModel(env);
			modelSubY.set(GRB_IntParam_OutputFlag, 0);
			modelSubY.set(GRB_StringAttr_ModelName, "Y_FINAl");
			modelSubY.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
			
			
			vec4GRBVar Y(I,vec3GRBVar(J,vec2GRBVar(T,vec1GRBvar(V))));
			vec3GRBVar Z(I,vec2GRBVar(T,vec1GRBvar(V)));
			vec3GRBVar Q(I,vec2GRBVar(J,vec1GRBvar(T)));
			
			for(i = 0; i < I ; i++)
				for(j = 0; j < J; j++)
					for (t = 0; t < T; t++){
						var_name="q("+to_string(i)+","+to_string(j)+","+to_string(t)+")";
						Q.at(i).at(j).at(t) = modelSubY.addVar(0.0, (t==T-1?0.0:GRB_INFINITY), -In.at(i).at(j).at(t), GRB_INTEGER, var_name);
					}
				
			
			for (t = 0; t < T; t++)
				for (v = 0; v < V; v++)
					for(i = 0; i < I ; i++){
						var_name="z("+to_string(i)+","+to_string(t)+","+to_string(v)+")";
						Z.at(i).at(t).at(v) = modelSubY.addVar(0.0, GRB_INFINITY, -cosf.at(i).at(t).at(v), GRB_INTEGER, var_name);
						for(j = 0; j < J; j++){
							var_name = "y("+to_string(i)+","+to_string(j)+","+to_string(t)+","+to_string(v)+")";
							Y.at(i).at(j).at(t).at(v) = modelSubY.addVar(0.0, (A.at(i).at(j).at(v)==0?0.0:GRB_INFINITY), -cos.at(i).at(j).at(v), GRB_INTEGER, var_name);
						}
					}
				
			
							
			
			 for(i = 0; i < I ; i++){
				for(t = 0; t < T; t++){
					for (v = 0; v < V; v++){
						
						
						double coe(ins.m.at(i).at(t).at(v));
						for (j = 0; j < J; j++) if(t<T) coe -= X_sol.at(i).at(j).at(t).at(v); 
						for (j = 0; j < J; j++) if((i != j) && ((t - ins.tau.at(j).at(i)) >= 0) ) coe += X_sol.at(j).at(i).at(t - ins.tau.at(j).at(i)).at(v);
						
						////for (j = 0; j < J; j++) if(t<T) coe -= X.at(i).at(j).at(t).at(v).get(GRB_DoubleAttr_X); 
						////for (j = 0; j < J; j++) if((i != j) && ((t - ins.tau.at(j).at(i)) >= 0) ) coe += X.at(j).at(i).at(t - ins.tau.at(j).at(i)).at(v).get(GRB_DoubleAttr_X);
						//cout << "\t*\t" << flush;
						
						//if(coe != 0) cout << "\n\tnode: " << flat_ixt(ins,t,i) << " v: " << v << "\t = " << coe << flush;
						
				
						var_name="flux("+to_string(i)+","+to_string(t)+","+to_string(v)+")";
						GRBLinExpr sum_expr = 0;
						for (j = 0; j < J; j++) sum_expr +=  Y[i][j][t][v];
						for (j = 0; j < J; j++) if((i != j) && ((t - tau.at(j).at(i)) >= 0) ) sum_expr -= (Y[j][i][t - tau.at(j).at(i)][v]);
						if (t > 0) sum_expr -= Y[i][i][t-1][v];
						sum_expr -= Z[i][t][v];
						modelSubY.addConstr(sum_expr, GRB_EQUAL ,coe, var_name);
					}
				}
			}
			
			for (i = 0; i < I; i++)
				for (t = 0; t < T; t++){
					
					double coe(0.0);
					
					var_name="cap("+to_string(i)+","+to_string(t)+")";

					for (v = 0; v < V; v++)
						for (j = 0; j < J; j++)
							 if((i != j) && ((t - tau.at(j).at(i)) >= 0) )
								 coe += X_sol.at(j).at(i).at(t - tau.at(j).at(i)).at(v);
					if( coe > Kap.at(i).at(t) ){
						cerr << "\nError in constraint " << var_name << " --- > " << setprecision(10) << coe << " > " << setprecision(10) << Kap.at(i).at(t) << flush;
						exit(EXIT_FAILURE);
					}
					
					
				}
			
				
			for(i = 0; i < I ; i++)
				for(j = 0; j < J; j++)
					for (t = 0; t < T; t++){
							var_name="backl("+to_string(i)+","+to_string(j)+","+to_string(t)+")";
							double coe(dem.at(i).at(j).at(t));
							GRBLinExpr sum_expr = 0;	 
							sum_expr += Q[i][j][t];
							if(t > 0) sum_expr -= Q[i][j][t-1];
							for(v = 0; v < V; v++) coe -= X_sol.at(i).at(j).at(t).at(v);
							modelSubY.addConstr(sum_expr, GRB_EQUAL , coe , var_name);
						}
			
			
			
			cout << endl;
			
			modelSubY.set(GRB_IntParam_Presolve, 0);
			modelSubY.write("Benders_Xfixed.lp");
			modelSubY.optimize();
			if(modelSubY.get(GRB_IntAttr_Status) == 2 ){
				for (v = 0; v < V; v++)
					for (t = 0; t < T; t++)
						for(i = 0; i < I ; i++)
							for(j = 0; j < J; j++)
							if( Y.at(i).at(j).at(t).at(v).get(GRB_DoubleAttr_X) > EPS4){
								cout << endl << Y.at(i).at(j).at(t).at(v).get(GRB_StringAttr_VarName) << " = " << Y.at(i).at(j).at(t).at(v).get(GRB_DoubleAttr_X) << flush;
								cout << setw(15) << flat_ixt(ins,t,i) << setw(10) << flat_ixt(ins,(i==j?t+1:t + ins.tau.at(i).at(j)),j) << flush;
								cout << setw(15) << -cos.at(i).at(j).at(v) << flush;
								cout << setw(15) << -cos.at(i).at(j).at(v)*Y.at(i).at(j).at(t).at(v).get(GRB_DoubleAttr_X) << flush;
							}
				cout << endl;			
				for (v = 0; v < V; v++)
					for (t = 0; t < T; t++)
						for(i = 0; i < I ; i++)
							if( Z.at(i).at(t).at(v).get(GRB_DoubleAttr_X) > EPS4 ){
								cout << endl << Z.at(i).at(t).at(v).get(GRB_StringAttr_VarName) << " = " << Z.at(i).at(t).at(v).get(GRB_DoubleAttr_X) << setw(10) << flat_ixt(ins,t,i) << flush;
								cout << setw(15) << ins.cosf.at(i).at(t).at(v) << flush;
								cout << setw(15) << ins.cosf.at(i).at(t).at(v)*( Z.at(i).at(t).at(v).get(GRB_DoubleAttr_X) - ins.m.at(i).at(t).at(v)) << flush;
								
							}
				cout << endl;
				for (t = 0; t < T; t++)
					for(i = 0; i < I ; i++)
						for(j = 0; j < J; j++)
							if( Q.at(i).at(j).at(t).get(GRB_DoubleAttr_X) > EPS4 ){
								cout << endl << Q.at(i).at(j).at(t).get(GRB_StringAttr_VarName) << " = " << Q.at(i).at(j).at(t).get(GRB_DoubleAttr_X) <<  flush;
								cout << setw(15) << flat_ixt(ins,t,i) << setw(10) << flat_ixt(ins,(i==j?t+1:t + ins.tau.at(i).at(j)),j) << flush;
								cout << setw(15) << ins.In.at(i).at(j).at(t) << flush;
								cout << setw(15) << ins.In.at(i).at(j).at(t)*Q.at(i).at(j).at(t).get(GRB_DoubleAttr_X)  << flush;
								
							}
						
							
			}
			if( modelSubY.get(GRB_IntAttr_Status) == 3 ){
				cerr << "\nSolution from the Master Problem is infeasible " << flush;
			}
	    }*/
		
		sol.close();
		cout << "\nBye from BFTVAP_solveBenders()\n" << endl;
		//cout << "\n" << ins.name << flush;
		
		#ifdef comparator
			return ZIP;
		#endif
		return 0;
		
		
		
  } catch(GRBException e) {
    cerr << "\n\nError code = " << e.getErrorCode() << "\t";
    cerr << flush;
    cerr << e.getMessage() << endl;
    exit(EXIT_FAILURE);
  } catch(...) {
    cerr << "\n\nException during optimization" << endl;
    cerr << flush;
    exit(EXIT_FAILURE);
  }
  
  
	
}
