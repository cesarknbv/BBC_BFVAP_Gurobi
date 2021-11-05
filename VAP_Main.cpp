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


