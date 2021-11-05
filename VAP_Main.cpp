#include "Functions/functions.cpp"
#include "Functions/BFVAP.cpp"
using namespace std;



int main(int   argc, char *argv[])
{

	Instance_data ins;
   
    
    
    #ifdef comparator
		
		
		resultados_OV << "\n" << fixed << showpoint << setprecision(4) << flush;
		
		
		while( !Generate_Instance(argc, argv, ins) )
		{
			cerr << "\nThe input files could not be read.\n" << flush;
			return 1;
		}
		double z  = BFTVAP_solveOriginal( ins );
		double za  = BFTVAP_solveBenders( ins );
		
		cout << "\nSemilla: " << setw(10) << semilla << flush;
		resultados_OV << fixed << showpoint << "\t\t" << flush;

		if ( abs(z-za) > EPS3 ){
			cout << "\nDifferent Results\n" << flush;
			/* Both EXIT_FAILURE and the value one indicate Unsuccessful program execution status */
			return EXIT_FAILURE;
			//return 1;
		}else{
			/* Both EXIT_SUCCESS and the value zero indicate successful program execution status */
			cout << "\nEqual Results\n" << flush;
			return EXIT_SUCCESS;
			//return 0;
		}
		
    #else
		
		
		#ifdef Primera_Vez
		Create_costs_for_existing_instances(argc, argv, ins);
		#else
		resultados_OV << "\n" << fixed << showpoint << setprecision(6) << flush;
		if ( argc == 12 ){
			if(DVAP_read_instance(argc, argv, ins) < 1)
			{
				cerr << "\nThe input files could not be read.\n" << flush;
				/* Both EXIT_SUCCESS and the value zero indicate successful program execution status */
				exit(EXIT_FAILURE);
			}
		}else{
				cerr << "\nIncorrect number of arguments.\n" << flush;
				/* Both EXIT_FAILURE and the value one indicate Unsuccessful program execution status */
				exit(EXIT_FAILURE);
		}   
		
		/* Solve BacklogFleet VAP */
		if( argc == 12 ){
			BFTVAP_solveOriginal( ins ); 
			////while(getchar() != '\n' );	
			BFTVAP_solveBenders( ins ); 
		}else{
			cerr << "\n -- Incorrect number of arguments for reading instance --" << endl;
			return EXIT_FAILURE;
		}
		

		return EXIT_SUCCESS; 
		#endif
		
		
	#endif
	
	
	
	
	
	
	
}




//// Set the TuneResults parameter to 1
//model.set(GRB_IntParam_TuneResults, 1);
//// Tune the model
//model.tune();
//// Get the number of tuning results
//int resultcount = model.get(GRB_IntAttr_TuneResultCount);
//if (resultcount > 0) {
  //// Load the tuned parameters into the model's environment
  //model.getTuneResult(0);
  //// Write tuned parameters to a file
  //model.write("tune.prm");
  //// Solve the model using the tuned parameters
  //model.optimize();
//}

//{
	    //// Create an environment
		//GRBEnv env = GRBEnv(true);
		//env.set("LogFile", "VAP.log");
		//env.start();
		
		//// Create an empty model
		//GRBModel model = GRBModel(env);
		//model.set(GRB_StringAttr_ModelName, "VAP");
		
		//// Turn off display and heuristics
		//model.set(GRB_IntParam_OutputFlag, 0);
		////model.set(GRB_DoubleParam_Heuristics, 0.0);
		
		//model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
		
		//GRBVar p1 = model.addVar(0.0, GRB_INFINITY,  1,  GRB_CONTINUOUS, "p1");
		//GRBVar p2 = model.addVar(0.0, GRB_INFINITY,  0,  GRB_CONTINUOUS, "p2");
		//GRBVar p3 = model.addVar(0.0, GRB_INFINITY,  0,  GRB_CONTINUOUS, "p3");
		//GRBVar p4 = model.addVar(0.0, GRB_INFINITY, -2,  GRB_CONTINUOUS, "p4");
		//GRBVar p5 = model.addVar(0.0, GRB_INFINITY,  0,  GRB_CONTINUOUS, "p5");
		//GRBVar p6 = model.addVar(0.0, GRB_INFINITY,  0,  GRB_CONTINUOUS, "p6");
		//GRBVar p7 = model.addVar(0.0, GRB_INFINITY,  1,  GRB_CONTINUOUS, "p7");
		//GRBVar p8 = model.addVar(0.0, GRB_INFINITY,  0,  GRB_CONTINUOUS, "p8");
		//GRBVar p9 = model.addVar(0.0, GRB_INFINITY,  0,  GRB_CONTINUOUS, "p9");
		//GRBVar pn = model.addVar(0.0, GRB_INFINITY, -2, GRB_CONTINUOUS,  "pn");
		//model.addConstr(pn-p9, GRB_LESS_EQUAL ,2, "c1");
		//model.addConstr(p9-p8, GRB_LESS_EQUAL ,2, "c2");
		//model.addConstr(p8-p7, GRB_LESS_EQUAL ,2, "c3");
		//model.addConstr(p7-p6, GRB_LESS_EQUAL ,2, "c4");
		//model.addConstr(p6-p5, GRB_LESS_EQUAL ,2, "c5");
		//model.addConstr(p5-p4, GRB_LESS_EQUAL ,2, "c6");
		//model.addConstr(p4-p3, GRB_LESS_EQUAL ,2, "c7");
		//model.addConstr(p3-p2, GRB_LESS_EQUAL ,2, "c8");
		//model.addConstr(p2-p1, GRB_LESS_EQUAL ,2, "c9");
		//model.set(GRB_IntParam_Presolve, 0);
		//model.optimize();
		//model.write("small.lp");
		
		//if(model.get(GRB_IntAttr_Status) == 2) {
			
			//cout << "\nOF: " << model.get(GRB_DoubleAttr_ObjVal) << flush;
			
			//cout << "\n" << p1.get(GRB_StringAttr_VarName) << " = " << p1.get(GRB_DoubleAttr_X) << flush;
			//cout << "\n" << p2.get(GRB_StringAttr_VarName) << " = " << p2.get(GRB_DoubleAttr_X) << flush;
			//cout << "\n" << p3.get(GRB_StringAttr_VarName) << " = " << p3.get(GRB_DoubleAttr_X) << flush;
			//cout << "\n" << p4.get(GRB_StringAttr_VarName) << " = " << p4.get(GRB_DoubleAttr_X) << flush;
			//cout << "\n" << p5.get(GRB_StringAttr_VarName) << " = " << p5.get(GRB_DoubleAttr_X) << flush;
			//cout << "\n" << p6.get(GRB_StringAttr_VarName) << " = " << p6.get(GRB_DoubleAttr_X) << flush;
			//cout << "\n" << p7.get(GRB_StringAttr_VarName) << " = " << p7.get(GRB_DoubleAttr_X) << flush;
			//cout << "\n" << p8.get(GRB_StringAttr_VarName) << " = " << p8.get(GRB_DoubleAttr_X) << flush;
			//cout << "\n" << p9.get(GRB_StringAttr_VarName) << " = " << p9.get(GRB_DoubleAttr_X) << flush;
			//cout << "\n" << pn.get(GRB_StringAttr_VarName) << " = " << pn.get(GRB_DoubleAttr_X) << flush;
			
		//}
		//if(model.get(GRB_IntAttr_Status) == 5) {
			
			//cout << "\n" << p1.get(GRB_StringAttr_VarName) << " = " << p1.get(GRB_DoubleAttr_UnbdRay) << flush;
			//cout << "\n" << p2.get(GRB_StringAttr_VarName) << " = " << p2.get(GRB_DoubleAttr_UnbdRay) << flush;
			//cout << "\n" << p3.get(GRB_StringAttr_VarName) << " = " << p3.get(GRB_DoubleAttr_UnbdRay) << flush;
			//cout << "\n" << p4.get(GRB_StringAttr_VarName) << " = " << p4.get(GRB_DoubleAttr_UnbdRay) << flush;
			//cout << "\n" << p5.get(GRB_StringAttr_VarName) << " = " << p5.get(GRB_DoubleAttr_UnbdRay) << flush;
			//cout << "\n" << p6.get(GRB_StringAttr_VarName) << " = " << p6.get(GRB_DoubleAttr_UnbdRay) << flush;
			//cout << "\n" << p7.get(GRB_StringAttr_VarName) << " = " << p7.get(GRB_DoubleAttr_UnbdRay) << flush;
			//cout << "\n" << p8.get(GRB_StringAttr_VarName) << " = " << p8.get(GRB_DoubleAttr_UnbdRay) << flush;
			//cout << "\n" << p9.get(GRB_StringAttr_VarName) << " = " << p9.get(GRB_DoubleAttr_UnbdRay) << flush;
			//cout << "\n" << pn.get(GRB_StringAttr_VarName) << " = " << pn.get(GRB_DoubleAttr_UnbdRay) << flush;
			
		//}
		//if(model.get(GRB_IntAttr_Status) == 3) {
				//cout << "\nInfeasible" << flush;
		//}
						
	    //exit(0);
	//}

//{
		//// Create an environment
		//GRBEnv env = GRBEnv(true);
		//env.set("LogFile", "mip1.log");
		//env.start();
		
		//// Create an empty model
		//GRBModel model = GRBModel(env);

		//// Create variables
		//GRBVar u0 = model.addVar(-GRB_INFINITY, GRB_INFINITY, -1.0, GRB_CONTINUOUS, "u0");
		//GRBVar u1 = model.addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "u1");
		//GRBVar u3 = model.addVar(-GRB_INFINITY, GRB_INFINITY, -1.0, GRB_CONTINUOUS, "u3");
		//GRBVar u4 = model.addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "u4");
		
		//model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
		
		//model.addConstr(- u0 <= 0, "c0");
		//model.addConstr(- u0 <= 2, "c1");
		//model.addConstr(- u1 <= 2, "c2");
		//model.addConstr(- u1 <= 0, "c3");
		//model.addConstr(- u3 <= 0, "c4");
		//model.addConstr(- u4 <= 0, "c5");
		
		//model.set(GRB_IntParam_Presolve, 0);
		////model.set(GRB_IntParam_OutputFlag, 0);
		
		//// Optimize model
        //model.optimize();
        
        //model.write("ex.lp");
        
        //cout  << "\nStatus: " << model.get(GRB_IntAttr_Status) << flush;
        
        //if(model.get(GRB_IntAttr_Status) == 2 ){
			//cout << "\nu0 =" << u0.get(GRB_DoubleAttr_X) << flush;
			//cout << "\nu1 =" << u1.get(GRB_DoubleAttr_X) << flush;
		//}
		//if(model.get(GRB_IntAttr_Status) == 5 ){
			//cout << "\nu0 =" << u0.get(GRB_DoubleAttr_UnbdRay) << flush;
			//cout << "\nu1 =" << u1.get(GRB_DoubleAttr_UnbdRay) << flush;
			//cout << "\nu3 =" << u3.get(GRB_DoubleAttr_UnbdRay) << flush;
			//cout << "\nu4 =" << u4.get(GRB_DoubleAttr_UnbdRay) << flush;
		//}
		
		
	//}
