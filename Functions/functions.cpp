#include "DVAP_par.hpp"



bool Create_costs_for_existing_instances( int argc, char *argv[],  Instance_data &ins ){
	cout << "\nHello from Create_costs_for_existing_instances()\n" << endl;
	 
	ifstream Datos(argv[10], ios::in);
	if( !Datos ){
		cerr << "File of Data Size could not be opened" << endl;
		exit(EXIT_FAILURE);
	}
	
	ifstream Cijv(argv[2], ios::in);
	if( !Cijv ){
		cerr << "File of Costs could not be opened" << endl;
		exit(EXIT_FAILURE);
	}
	int x, y, z;
	
	#ifdef Rej
		
	aux_name = argv[1];
	
	//cout << "Hola: " << aux_name << flush;
	size_t first_pos = aux_name.find_first_of("p");
	size_t sec_pos = aux_name.find_first_of("r");
	aux_name = aux_name.substr (first_pos,(sec_pos-first_pos));
	//cout << "Hola: " << aux_name << flush;
	
	#endif

 
    ////////////////////////////////////////////////////////////////////////////
    
   

    int number_of_demand, restrictions_per_vehicle, num_veh;
    Datos >> ins.I;
    Datos >> ins.T;
    Datos >> ins.V;
    Datos >> num_veh;
    Datos >> number_of_demand;
    Datos >> restrictions_per_vehicle;
    
    
    ins.name += to_string(ins.I);ins.name += "x"; ins.name += to_string(ins.T); ins.name += "x"; ins.name += to_string(ins.V); ins.name += "x"; ins.name += to_string(num_veh); 
    ins.name += "x"; ins.name += to_string(number_of_demand); ins.name += "x"; ins.name += to_string(restrictions_per_vehicle);
    

    
    Datos.close();
    #ifdef Rej
	ins.name = ins.name.substr(0,ins.name.size()-1);
	ins.name = ins.name + aux_name;
	#endif
	
    cout << "\nInstance: " << ins.name << flush;
    //exit(0);
    int I = ins.I , J = ins.J = ins.I , T = ins.T , V = ins.V;
    
    
    cout << "\n\tI:" << ins.I << " J:" << ins.J << " T:" << ins.T << " V:" << ins.V << flush;
    ins.cos = vec3Dou( I , vec2Dou( J , vec1Dou( V )));
   
   	for (y = 0; y < J; y++) 
   		for (x = 0; x < I; x++)
   			Cijv >> ins.cos.at(x).at(y).at(0);
   		
   	    
   	
   	for (y = 0; y < J; y++) 
   		for (x = 0; x < I; x++) 
			for ( z = 1; z < V; z++)
				ins.cos.at(x).at(y).at(z) = ins.cos.at(x).at(y).at(0);
			
		
   	Cijv.close();	
    

    string scosf( (ins.name + "cosf.dat" ) );
    string sH(    (ins.name + "H.dat"   ) ); 
	string sK(    (ins.name + "Kap.dat"   ) ); 
	cout << "\nOK1" << endl;
    /* --------------------------------------------------------------- 
	 * 					Fleet Sizing Costs
	 * --------------------------------------------------------------- */
	ofstream instanceData( scosf.c_str() , ios::out );
	if( !instanceData )
    {
       cerr << "File could not be opened" << endl;
       return false;
    } 
	
	ins.cosf = vec3Dou( I, vec2Dou( T, vec1Dou(V , 0.0 ) ) );
	for (int i = 0; i < I; i++){ 
		for (int t = 0; t < T; t++){ 
			for (int v = 0; v < V; v++) {
	
				
				instanceData << i+1 << ";";
				instanceData << t+1 << ";";
				instanceData << v+1 << ";";
				double cushto(0.0);
				for (int j = 0; j < I; j++){ cushto += ins.cos.at(i).at(j).at(v); }
				cushto = cushto/I;
				cushto = ceil(((rand() / (RAND_MAX + 1.))*2)*(cushto));
				if ( cushto < EPS2 && cushto > -EPS2 ){
					
					cushto = 1+rand()%50;//cout << "\n" << cushto<< flush;
				}
				instanceData << cushto << endl;
				
				ins.cosf.at(i).at(t).at(v) = cushto;
		
	}}}
	
	instanceData.close();
	cout << "\nOK2" << endl;
	/* --------------------------------------------------------------- 
	 * 					Backlog Costs
	 * --------------------------------------------------------------- */
	instanceData.open( sH.c_str() , ios::out );
	if( !instanceData )
    {
       cerr << "File could not be opened" << endl;
       return false;
    } 
	
	ins.In = vec3Dou( I, vec2Dou( J, vec1Dou(T , 0.0 ) ) );
	for (int i = 0; i < I; i++){ 
		for (int j = 0; j < I; j++){ 
			for (int t = 0; t < T; t++){ 
	
				if( i != j ){
					instanceData << i+1 << ";";
					instanceData << j+1 << ";";
					instanceData << t+1 << ";";
					double backo( 0.5 + (rand()%3) );
					instanceData << backo << endl;
				
				    ins.In.at(i).at(j).at(t) = backo;
				}
	}}}
	
	instanceData.close();
    cout << "\nOK3" << endl;
	/* --------------------------------------------------------------- 
	 * 					Terminals Capacity
	 * --------------------------------------------------------------- */
	ins.Kap = vec2Int( I , vec1Int( T , 0) );
	instanceData.open( sK.c_str() , ios::out );
	if( !instanceData )
    {
       cerr << "File could not be opened" << endl;
       return false;
    } 
	
	
	for (int i = 0; i < I; i++){ 
		for (int t = 0; t < T; t++){ 

	
				
				instanceData << i+1 << ";";
				instanceData << t+1 << ";";
				int Kapo( 5 + (rand()%5) );
				
				instanceData << Kapo << endl;
				
				ins.Kap.at(i).at(t) = Kapo;
		
	}}
	
	instanceData.close();
	 
	cout << "\nBye from Create_costs_for_existing_instances()\n" << endl;
	return true;
	
}

int DVAP_read_instance(int argc, char *argv[],  Instance_data &ins)
{
	
    cout << "\nHello from VAP_ReadInstnace()\n" << endl;
    
   
    
    ifstream Datos(argv[10], ios::in);
	if( !Datos ){
		cerr << "File of Data Size could not be opened" << endl;
		exit(EXIT_FAILURE);
	}
    ifstream Tauij(argv[1], ios::in);
	if( !Tauij ){
		cerr << "File of Distances could not be opened" << endl;
		exit(EXIT_FAILURE);
	}
	ifstream Cijv(argv[2], ios::in);
	if( !Cijv ){
		cerr << "File of Costs could not be opened" << endl;
		exit(EXIT_FAILURE);
	}
	ifstream Cv(argv[3], ios::in);
	if( !Cv ){
		cerr << "File of CostFleet could not be opened" << endl;
		exit(EXIT_FAILURE);
	}
	ifstream Pijv(argv[4], ios::in);
	if( !Pijv ){
		cerr << "File of Profits could not be opened" << endl;
		exit(EXIT_FAILURE);
	}
	ifstream Demijt(argv[5], ios::in);
	if( !Demijt ){
		cerr << "File of Demands could not be opened" << endl;
		exit(EXIT_FAILURE);
	}
	ifstream Mitv(argv[6], ios::in);
	if( !Mitv ){
		cerr << "File of Supplies could not be opened" << endl;
		exit(EXIT_FAILURE);
	}
	ifstream Aijv(argv[7], ios::in);
	if( !Aijv ){
		cerr << "File of Movement Restrictions could not be opened" << endl;
		exit(EXIT_FAILURE);
	}
	ifstream Hijt(argv[8], ios::in);
	if( !Hijt ){
		cerr << "File of Backlog Costs could not be opened" << endl;
		exit(EXIT_FAILURE);
	}
	ifstream Kit(argv[9], ios::in);
	if( !Kit ){
		cerr << "File of Terminals Capacity could not be opened" << endl;
		exit(EXIT_FAILURE);
	}
	
	int x, y, z;
	
	#ifdef Rej
		
	aux_name = argv[1];
	
	//cout << "Hola: " << aux_name << flush;
	size_t first_pos = aux_name.find_first_of("p");
	size_t sec_pos = aux_name.find_first_of("r");
	aux_name = aux_name.substr (first_pos,(sec_pos-first_pos));
	//cout << "Hola: " << aux_name << flush;
	
	#endif

 
    ////////////////////////////////////////////////////////////////////////////
    
   

    int number_of_demand, restrictions_per_vehicle, num_veh;
    Datos >> ins.I;
    Datos >> ins.T;
    Datos >> ins.V;
    Datos >> num_veh;
    Datos >> number_of_demand;
    Datos >> restrictions_per_vehicle;
    
    
    ins.name += to_string(ins.I);ins.name += "x"; ins.name += to_string(ins.T); ins.name += "x"; ins.name += to_string(ins.V); ins.name += "x"; ins.name += to_string(num_veh); 
    ins.name += "x"; ins.name += to_string(number_of_demand); ins.name += "x"; ins.name += to_string(restrictions_per_vehicle);
    

    
    Datos.close();
    #ifdef Rej
	ins.name = ins.name.substr(0,ins.name.size()-1);
	ins.name = ins.name + aux_name;
	#endif
	
	op_sol_LP.open( "Data/"+ins.name+"_LP_sol", ios_base::out);
    op_sol_IP.open( "Data/"+ins.name+"_IP_sol", ios_base::out);
    op_solV_LP.open( "Data/"+ins.name+"_LP_veh", ios_base::out);
    op_solV_IP.open( "Data/"+ins.name+"_IP_veh", ios_base::out);
    
    
    
    
    
    cout << "\nInstance: " << ins.name << flush;
    //exit(0);
    int I = ins.I , J = ins.J = ins.I , T = ins.T , V = ins.V;
    
    
    cout << "\n\tI:" << ins.I << " J:" << ins.J << " T:" << ins.T << " V:" << ins.V << flush;
    
    
    	  
    
    resultados_OV << setw(30) << ins.name << "\t\t" << flush;
   
     ins.tau = vec2Int( I , vec1Int( J ));
     for (y = 0; y < I; y++)
       for (x = 0; x < J; x++) 
         Tauij >> ins.tau.at(x).at(y);
   
			
   
   
      Tauij.close();
    
     cout << "\nOKTau" << flush;

    ins.cos = vec3Dou( I , vec2Dou( J , vec1Dou( V )));
   
   	for (y = 0; y < J; y++){ 
   		for (x = 0; x < I; x++){ 
   			Cijv >> ins.cos.at(x).at(y).at(0);
   		}
   	}    
   	
   	for (y = 0; y < J; y++){ 
   		for (x = 0; x < I; x++){ 
			for ( z = 1; z < V; z++){
				ins.cos.at(x).at(y).at(z) = ins.cos.at(x).at(y).at(0);
			}
		}
	}
	
	
   

   
   	Cijv.close();	
   	cout << "\nOKCos" << flush;
   
   ins.cosf = vec3Dou( I , vec2Dou( T , vec1Dou( V )));
   y=-1,x=-1,z=-1; int itt(0);char vall;
   while ( !Cv.eof() ){
		if (itt%2 == 0){
			//if (Cv.eof() ) break;
			if( itt == 0 ){
				Cv >> y;
			}else if(itt == 2){
				Cv >> x;
			}else if(itt == 4){
				Cv >> z;
			}else{
				Cv >> ins.cosf.at(y-1).at(x-1).at(z-1);
			}
			//cout << num;
		}else{
			Cv >> vall;
			//if (Demijt.eof() ) break;
			//cout << val;
		}
		
		
		if ( ++itt == 7 ){ 
			itt = 0; //cout << endl;
			y=x=z=-1;
		}
	}
  
	  
   	
   	
   

   
   	Cv.close();	
    cout << "\nOKCof" << flush;
   
    

    ins.pro = vec3Dou( I , vec2Dou( J , vec1Dou( V )));
   
   
   	for (y = 0; y < J; y++){ 
   		for (x = 0; x < I; x++){ 
   			Pijv >> ins.pro.at(x).at(y).at(0);
   		}
   	}
    
    for (y = 0; y < J; y++){ 
   		for (x = 0; x < I; x++){ 
			for ( z = 1; z < V; z++){
				ins.pro.at(x).at(y).at(z) = ins.pro.at(x).at(y).at(0);
			}
		}
	}
	
	
   

   
  
   	Pijv.close();
   	cout << "\nOKPro" << flush;
   
    ins.dem = vec3Int( I , vec2Int( J , vec1Int( T, 0 )));
 

	int it = 0; char val; int num;
	
	
	int i,j,t,dem;
	while ( !Demijt.eof() ){
		if (it%2 == 0){
			Demijt >> num;
			//if (Demijt.eof() ) break;
			if( it == 0 ){
				i = num;
			}else if(it == 2){
				j = num;
			}else if(it == 4){
				t = num;
			}else{
				dem = num;
			}
			//cout << num;
		}else{
			Demijt >> val;
			//if (Demijt.eof() ) break;
			//cout << val;
		}
		
		
		if ( ++it == 7 ){ 
			it = 0; //cout << endl;
			ins.dem.at(i-1).at(j-1).at(t-1) += dem;
			i = j = t = dem =0;
		}
	}
	Demijt.close();
    cout << "\nOKDem" << flush;
    ins.m = vec3Int( I , vec2Int( T , vec1Int( V, 0 )));


    it = 0; num = 0;
    

    i = t = 0; 
    int v, m;
	 while ( !Mitv.eof() ){
		 if (it%2 == 0){
			 Mitv >> num;
			 //if (Mitv.eof() ) break;
			 if( it == 0 ){
				 i = num;
			 }else if(it == 2){
				 t = num;
			 }else if(it == 4){
				 v = num;
			 }else{
				 m = num;
			 }
			 //cout << num;
		 }else{
			 Mitv >> val;
			 //if (Mitv.eof() ) break;
			 //cout << val;
		 }
		 
		
		 if ( ++it == 7 ){ 
			 it = 0; //cout << endl;
			 (ins.m.at(i-1).at(t-1).at(m-1))++;
			 i = t = v = m = 0;
		 }
		 
	 }

	Mitv.close();
    cout << "\nOKM" << flush;


   //////////////////////////////////////////////////////////////////////////////////////////////
   ins.In = vec3Dou( I , vec2Dou( J , vec1Dou( T )));
   y=-1,x=-1,z=-1; itt=0;
   while ( !Hijt.eof() ){
		if (itt%2 == 0){
			//if (Hijt.eof() ) break;
			if( itt == 0 ){
				Hijt >> y;
			}else if(itt == 2){
				Hijt >> x;
			}else if(itt == 4){
				Hijt >> z;
			}else{
				Hijt >> ins.In.at(y-1).at(x-1).at(z-1);
			}
			//cout << num;
		}else{
			Hijt >> vall;
			//if (Demijt.eof() ) break;
			//cout << val;
		}
		
		
		if ( ++itt == 7 ){ 
			itt = 0; //cout << endl;
			y=x=z=-1;
		}
	}
  
	  
   	Hijt.close();	
    cout << "\nOKIn" << flush;
   //////////////////////////////////////////////////////////////////////////////////////////////
     
   ins.Kap = vec2Int( I , vec1Int( T ));
   y=-1,x=-1; itt=0;
   while ( !Kit.eof() ){
		if (itt%2 == 0){
			//if (Kit.eof() ) break;
			if( itt == 0 ){
				Kit >> y;
			}else if(itt == 2){
				Kit >> x;
			}else{
				Kit >> ins.Kap.at(y-1).at(x-1);
			}
			//cout << num;
		}else{
			Kit >> vall;
			//if (Demijt.eof() ) break;
			//cout << val;
		}
		
		
		if ( ++itt == 5 ){ 
			itt = 0; //cout << endl;
			y=x=-1;
		}
	}
  
	  
   	Kit.close();	
    cout << "\nOKKap" << flush;
   
     
   
    //  /////////////////////////////////////////////////////////////////////////////
    ins.A = vec3Int( I , vec2Int( J , vec1Int( V, 1 )));


		
	
	it = 0; num = 0;
	
	i = j = v = 0; 
	int a;
	#ifndef Rej
	while ( !Aijv.eof() ){
		if (it%2 == 0){
			Aijv >> num;
			//if (Aijv.eof() ) break;
			if( it == 0 ){
				i = num;
			}else if(it == 2){
				j = num;
			}else if(it == 4){
				v = num;
			}else{
				a = num;
			}
			//cout << num;
		}else{
			Aijv >> val;
			//if (Aijv.eof() ) break;
			//cout << val;
		}
		
		
		if ( ++it == 7 ){ 
			it = 0; //cout << endl;
			ins.A.at(i-1).at(j-1).at(v-1) = a;
			i = j = v = a = 0;
		}
	}
	#else
	/* ********************************RejaneInstances**************************************** */
		while ( !Aijv.eof() ){
			if (it%2 == 0){
				Aijv >> num;
				//if (Aijv.eof() ) break;
				if( it == 0 ){
					v = num;
				}else if(it == 2){
					i = num;
				}else if(it == 4){
					j = num;
				}else{
					a = num;
				}
				//cout << num;
			}else{
				Aijv >> val;
				//if (Aijv.eof() ) break;
				//cout << val;
			}
			
			
			if ( ++it == 7 ){ 
				it = 0; //cout << endl;
				ins.A.at(i-1).at(j-1).at(v-1) = a;
				i = j = v = a = 0;
			}
		}
		/* ********************************RejaneInstances**************************************** */
		#endif
		Aijv.close();
		cout << "\nOKA" << flush;
		
		
		ofstream output ("Parameter_VAP.txt", ios_base::out);
		if( !output ){
			cerr << "File of output could not be created" << endl;
			exit(1);
		}
		 for (y = 0; y < I; y++){
		   output << endl;
		   for (x = 0; x < J; x++){ 
				output <<  "tau["<< x << "][" << y << "]=" << setprecision(2) << ins.tau.at(x).at(y) << "  ";
			}
		 }
		 output << endl << endl;
		 for (z = 0; z < V; z++) {
			output << endl;
			for(y = 0; y < I; ++y){
			output << endl;
				for(x = 0; x < I; ++x){	   
					output << "pro["<< x << "][" << y << "][" << z << "]=" << setprecision(2) << ins.pro.at(x).at(y).at(z) << "  ";
				}
			}
		}
		output << endl << endl;
		for (z = 0; z < V; z++) {
			output << endl;
			for(y = 0; y < I; ++y){
			output << endl;
				for(x = 0; x < I; ++x){	   
					output << "cos["<< x << "][" << y << "][" << z << "]=" << setprecision(2) << ins.cos.at(x).at(y).at(z) << "  ";
				}
			}
		}
		output << endl << endl;
		for (z = 0; z < V; z++) {
			output << endl;
			for (y = 0; y < I; y++){ 
				for (x = 0; x < T; x++){ 
					output << "\ncosf["<< y << "]["<< x << "][" << z << "]=" << setw(5) << ins.cosf.at(y).at(x).at(z)  << setw(10) << flat_ixt(ins,x,y ) << flush;
				}
			}
		}
        output << endl << endl;
		for (v = 0; v < V; v++) {
			for (i = 0; i < I; i++){ 
				for (t = 0; t < T; t++) {  
					if ( ins.m.at(i).at(t).at(v) != 0 ){
					output << "\nm["<< i << "][" << t << "][" << v << "]=" 
						<< setw(5) << setprecision(2) << ins.m.at(i).at(t).at(v) << setw(15) << flat_ixt(ins,t,i) << flush;
					}
				}
			}
		 }
		 output << endl << endl;
		 for (t = 0; t < T; t++) {
			for (i = 0; i < I; i++){ 
				for (j = 0; j < J; j++){ 
	   				output << "\nh["<< i << "][" << j << "][" << t << "]=" << setw(5) << fixed << showpoint << setprecision(1) << ins.In.at(i).at(j).at(t) 
	   				<< setw(15) << flat_ixt(ins,t,i) << setw(10) << flat_ixt(ins,(i==j?t+1:t+ins.tau.at(i).at(j)),j);
	   			}
	   		 }
	   	 }
	   	 output << endl << endl;
	   	  for (i = 0; i < I; i++){ 
			//output << endl;
			for (t = 0; t < T; t++) {  
				output << "\nK["<< i << "][" << t << "]=" << setw(5) << fixed << showpoint << setprecision(1) << ins.Kap.at(i).at(t) 
				<< setw(15) << flat_ixt(ins,t,i);
			}	
		}
		output << endl << endl;
		 for (t = 0; t < T; t++) {
			for (i = 0; i < I; i++){ 
				for (j = 0; j < J; j++){ 
					if ( ins.dem.at(i).at(j).at(t) != 0 ){
						output << "\ndem["<< i << "][" << j << "][" << t << "]=" 
							<< setw(5) << ins.dem.at(i).at(j).at(t) << setw(15) << flat_ixt(ins,t,i) << setw(10) << flat_ixt(ins,(i==j?t+1:t+ins.tau.at(i).at(j)),j);
					}
				}
			}
		 }
         output << endl << endl;
		 for (z = 0; z < V; z++) {
		 	output << endl;
		 	for (y = 0; y < I; y++){ 
		 		for (x = 0; x < J; x++){ 
		 			if ( ins.A.at(y).at(x).at(z) == 0 ){
		 			output << "\nA["<< y << "][" << x << "][" << z << "]=" 
		 				<< setw(5) << ins.A.at(y).at(x).at(z);
		 			}
		 		}
		 	}
		 }
	     output << endl << endl;
	cout << "\nChao from VAP_ReadInstnace()\n" << endl;
   	

	return 1;
}

bool Generate_Instance(int argc, char *argv[], Instance_data &instance ){
	
	
	cout << "\nHello from Generate_Instance\n" << flush;
	
	vec3Dou &pro     = instance.pro;
	vec3Dou &cos     = instance.cos;
	vec3Dou &cosf    = instance.cosf;
	vec2Int &tau     = instance.tau;
	vec3Int &m       = instance.m;
	vec3Int &dem     = instance.dem;
	vec3Int &A       = instance.A;
	vec3Dou &In      = instance.In;
	vec2Int &Kap     = instance.Kap;

    //int I = instance.I , J = instance.I , T = instance.T , V = instance.V;
  
	int terminals(25), periods(25), vehicles(25), num_vehicles(130), num_cargas(130), p(0);
	
	
	if ( (argc >= 7) ){
		terminals = atoi(argv[1]); periods = atoi(argv[2]); vehicles = atoi(argv[3]); num_vehicles = atoi(argv[4]); num_cargas = atoi(argv[5]); p = atoi(argv[6]);
	}else if ( argc == 6 ){
		terminals = atoi(argv[1]); periods = atoi(argv[2]); vehicles = atoi(argv[3]); num_vehicles = vehicles; 
		p = (terminals*terminals)/10;
	}else if(argc == 5){
		terminals = atoi(argv[1]); periods = atoi(argv[2]); vehicles = atoi(argv[3]); num_vehicles = vehicles; 
		num_cargas = (terminals*terminals*periods)/4;
		p = (terminals*terminals)/10;
	}else{
		
	}
	
	
	
	instance.name = to_string(terminals) + "x" + to_string(periods) + "x" + to_string(vehicles) + "x" + to_string(num_vehicles) + "x" + to_string(num_cargas) + "x" + to_string(p);
	
	
	
	op_sol_LP.open( "Data/"+instance.name+"_LP_sol", ios_base::out);
    op_sol_IP.open( "Data/"+instance.name+"_IP_sol", ios_base::out);
    op_solV_LP.open( "Data/"+instance.name+"_LP_veh", ios_base::out);
    op_solV_IP.open( "Data/"+instance.name+"_IP_veh", ios_base::out);
	
	string stau(  (instance.name + "dist.dat") ); 
	string  scos( (instance.name + "cos.dat" ) ); 
	string  scosf( (instance.name + "cosf.dat" ) ); 
	string spro(  (instance.name + "pro.dat" ) ); 
	string sm(    (instance.name + "m.dat"   ) ); 
	string sdem(  (instance.name + "dem.dat" ) ); 
	string sA(    (instance.name + "A.dat"   ) ); 
	string sH(    (instance.name + "H.dat"   ) ); 
	string sK(    (instance.name + "Kap.dat"   ) ); 
	string ss(    (instance.name + ".dat"    ) );
	
	cout << "\nTerminals = "		    << terminals
		 << "\nPeriods = " 		   		<< periods
		 << "\nVehicles = " 			<< vehicles 
		 << "\nNumber of Vehicles = "   << num_vehicles
		 << "\nNumber of Demands = "    << num_cargas
		 << "\nP = " 					<< p
		 << flush;
		 
		 
    if ( (terminals*terminals) < p ){
		cerr << "\nError in logic of input data, since the number of trips prohibited " << p << " can not exceed the nuumber of pair (terminalxterminal) " << (terminals*terminals) << endl;
		exit(EXIT_FAILURE);
	}
	
	if ( (terminals*terminals*periods) < num_cargas ){
		cerr << "\nError in logic of input data, since the number of freights " << num_cargas << " can not exceed the nuumber of pair (terminalxterminalxperiods) " << (terminals*terminals*periods) << endl;
		exit(EXIT_FAILURE);
	}
	if ( num_vehicles < vehicles ){
		cerr << "\nError in logic of input data, since the number of types of vehicles " << vehicles << " can not exceed the nuumber of vehicles " << num_vehicles << endl;
		exit(EXIT_FAILURE);
	}
		 
	resultados_OV << setw(3) << terminals << " " << setw(4) << periods << " " << setw(4) << vehicles << " " << setw(4) << num_vehicles << " " << setw(4) << num_cargas << " " << setw(4) << p <<"\t\t";
	
	instance.I = instance.J = terminals;
	instance.T = periods;
	instance.V = vehicles;
	
	cout << "\nI: " << instance.I << " J: " << instance.J << " T: " << instance.T << " V: " << instance.V << flush;
	
	int I(instance.I), J(instance.J),T(instance.T),V(instance.V);
	
	if( argc == 8 ){
		semilla = atoi(argv[7]);
	}else{
		semilla = time(0);
    }

	//semilla = 1616398604;
	
	srand(semilla);
	resultados_OV << setw(10) << semilla << "\t" << flush;
	cout << "\nSemilla: " << setw(10) << semilla << flush;
	
	
	
    
    	
	/* --------------------------------------------------------------- 
	 * 						Restriction of Movements
	 * --------------------------------------------------------------- */
	cout << "\nOK1" << flush;
	
	
	ofstream instanceData( sA.c_str() , ios::out );
	if( !instanceData )
    {
       cerr << "File could not be opened" << endl;
       return false;
    } 
    
    A = vec3Int( terminals, vec2Int( terminals, vec1Int(vehicles , 1 ) ) );
    for ( int vv = 0; vv < vehicles; vv++  ){
		
		for ( int jj = 0; jj < p; jj++ ){
			int k;
			int i = 1 + rand()%terminals;
			instanceData << i << ";";
			do{
				k = 1 + rand()%terminals;
			}while ( i == k );
			instanceData << k << ";";
			instanceData << vv+1 << ";";
			instanceData << 0 << endl;
			
			A.at(i-1).at(k-1).at(vv) = 0;
			
		}
		
	}

    instanceData.close();
    
    /* ****************************************************************************** */
    /* ****************************************************************************** */
    /* ****************************************************************************** */
    /* ****************************************************************************** */
    //for (int z = 0; z < V; z++) {
		//for (int y = 0; y < I; y++){ 
			//A.at(y).at(y).at(z) = 0;
		//}
	//}
    /* ****************************************************************************** */
    /* ****************************************************************************** */
    /* ****************************************************************************** */
    /* ****************************************************************************** */
    
    
    
    
    
   /* --------------------------------------------------------------- 
	 * 					Distances
	 * --------------------------------------------------------------- */
    cout << "\tOK2" << flush;
    instanceData.open( stau.c_str(), ios::out );
	if( !instanceData )
    {
       cerr << "File could not be opened" << endl;
       return false;
    } 

	tau = vec2Int( terminals , vec1Int( terminals , 0) );
	vec1Int x_cor( terminals );
	vec1Int y_cor( terminals );
	
	
	bool flag;
    do{
		flag = true;
		for( int k = 0; k < terminals; k++ ){
			x_cor.at(k) = 1 + rand()%(periods);
			y_cor.at(k) = 1 + rand()%(periods);
			//cout << "\n" << x[k] << "   " << y[k];
		}
	
	
		for(int i = 0; i < terminals; i++ ){
			for(int j = i; j < terminals; j++ ){
				if( i == j ){
					tau.at(i).at(j) = 0;
				}else{
					double val = (x_cor.at(i)^2)+(y_cor.at(j)^2);
					tau.at(i).at(j) = 1/*sqrt(val)*/ ;
					tau.at(j).at(i) = 1/*tau.at(i).at(j)*/;
				}
			}
		}
		
		/* Check triangular inequality for distance matrix*/
		//for (int co = 0; co < terminals; co++ ){
			//for(int i = 0; i < terminals; i++ ){
				//for(int j = i; j < terminals; j++ ){
					
					//if( tau.at(i).at(j) > ( tau.at(i).at(co) + tau.at(co).at(j) ) ) flag = false;
					
				//}
			//}
		//}
	}while( flag == false );
	
	


   for(int j = 0; j < terminals; ++j){
		if ( j > 0 ) {instanceData << endl; }
		for(int i = 0; i < terminals; ++i){	   
			if ( (i != j)  && (tau.at(i).at(j) == 0)  ){
				tau.at(i).at(j) = 1;
			}
			if( tau.at(i).at(j) < 10 ){ instanceData << tau.at(i).at(j) << "   ";}else{
				instanceData << tau.at(i).at(j) << "  ";}
		}
	}
	
	
	instanceData.close();
	
	/* --------------------------------------------------------------- 
	 * 					Empty vehicle Costs
	 * --------------------------------------------------------------- */
	 cout << "\tOK3" << flush;
	instanceData.open( scos.c_str(), ios::out );
	if( !instanceData )
    {
       cerr << "File could not be opened" << endl;
       return false;
    } 

	
	cos = vec3Dou( terminals, vec2Dou( terminals, vec1Dou(vehicles , 0.0 ) ) );
	vec1Dou maximo(vehicles,0.0);
    
	for(int vv = 0; vv < vehicles; vv++ ){
		
		int mult( 1+ rand()%20 );
		for(int i = 0; i < terminals; i++ ){
			for(int j = i; j < terminals; j++ ){
				
				
				cos.at(i).at(j).at(vv) = tau.at(i).at(j)*mult;
				cos.at(j).at(i).at(vv) = cos.at(i).at(j).at(vv);
				
				
				if ( cos.at(i).at(j).at(vv) > maximo.at(vv) ) maximo.at(vv) = cos.at(i).at(j).at(vv);
				
			}
		}
	}


	for(int vv = 0; vv < vehicles; vv++ ){
		if ( vv > 0 ) {instanceData << endl;}
		for(int j = 0; j < terminals; ++j){
			if ( j > 0 ) {instanceData << endl; }
			for(int i = 0; i < terminals; ++i){	  
				if( cos.at(i).at(j).at(vv) < 10 ){ 
					instanceData << fixed << showpoint << setprecision(1) << cos.at(i).at(j).at(vv) << "   ";
				}
				else{
					instanceData << fixed << showpoint << setprecision(1) << cos.at(i).at(j).at(vv) << "  ";
				}
			}
		}
	}
	
	
	
	
	instanceData.close();
	
	/* --------------------------------------------------------------- 
	 * 					Fleet Sizing Costs
	 * --------------------------------------------------------------- */
	cout << "\tOK4" << flush;
	instanceData.open( scosf.c_str() , ios::out );
	if( !instanceData )
    {
       cerr << "File could not be opened" << endl;
       return false;
    } 
	
	cosf = vec3Dou( I, vec2Dou( T, vec1Dou(V , 0.0 ) ) );
	for (int i = 0; i < I; i++){ 
		for (int t = 0; t < T; t++){ 
			for (int v = 0; v < V; v++) {
	
				
				instanceData << i+1 << ";";
				instanceData << t+1 << ";";
				instanceData << v+1 << ";";
				double cushto(0.0);
				for (int j = 0; j < I; j++){ cushto += instance.cos.at(i).at(j).at(v); }
				cushto = cushto/I;
				cushto = ceil(((rand() / (RAND_MAX + 1.))*2)*(cushto));
				if ( cushto < EPS2 && cushto > -EPS2 ){
					
					cushto = 1+rand()%50;//cout << "\n" << cushto<< flush;
				}
				instanceData << cushto << endl;
				
				instance.cosf.at(i).at(t).at(v) = cushto;
		
	}}}
	
	instanceData.close();
	/* --------------------------------------------------------------- 
	 * 					Backlog Costs
	 * --------------------------------------------------------------- */
	cout << "\tOK5" << flush;
	instanceData.open( sH.c_str() , ios::out );
	if( !instanceData )
    {
       cerr << "File could not be opened" << endl;
       return false;
    } 
	
	In = vec3Dou( I, vec2Dou( J, vec1Dou(T , 0.0 ) ) );
	for (int i = 0; i < I; i++){ 
		for (int j = 0; j < I; j++){ 
			for (int t = 0; t < T; t++){ 
	
				if( i != j ){
					instanceData << i+1 << ";";
					instanceData << j+1 << ";";
					instanceData << t+1 << ";";
					double backo( 0.5 + (rand()%3) );
					instanceData << backo << endl;
				
				    instance.In.at(i).at(j).at(t) = backo;
				}
	}}}
	
	instanceData.close();

	/* --------------------------------------------------------------- 
	 * 					Terminals Capacity
	 * --------------------------------------------------------------- */
	cout << "\tOK6" << flush;
	Kap = vec2Int( I , vec1Int( T , 0) );
	instanceData.open( sK.c_str() , ios::out );
	if( !instanceData )
    {
       cerr << "File could not be opened" << endl;
       return false;
    } 
	
	
	for (int i = 0; i < I; i++){ 
		for (int t = 0; t < T; t++){ 

	
				
				instanceData << i+1 << ";";
				instanceData << t+1 << ";";
				int Kapo( 5 + (rand()%5) );
				
				instanceData << Kapo << endl;
				
				instance.Kap.at(i).at(t) = Kapo;
		
	}}
	
	instanceData.close();

	/* --------------------------------------------------------------- 
	 * 					Profits
	 * --------------------------------------------------------------- */
	cout << "\tOK7" << flush;
	instanceData.open( spro.c_str(), ios::out );
	if( !instanceData )
    {
       cerr << "File could not be opened" << endl;
       return false;
    } 

	
   pro = vec3Dou( terminals, vec2Dou( terminals, vec1Dou(vehicles , 0.0 ) ) );

   for(int vv = 0; vv < vehicles; vv++ ){
		for(int i = 0; i < terminals; i++ ){
			for(int j = i; j < terminals; j++ ){
				if( i == j ){
					pro.at(i).at(j).at(vv) = 0;
				}else{
					pro.at(i).at(j).at(vv) = maximo.at(vv) + rand()%20;
					pro.at(j).at(i).at(vv) = pro.at(i).at(j).at(vv);
				}			    
				
			}
		}
	}


	for(int vv = 0; vv < vehicles; vv++ ){
		if ( vv > 0 ) {instanceData << endl;}
		for(int j = 0; j < terminals; ++j){
			if ( j > 0 ) {instanceData << endl; }
			for(int i = 0; i < terminals; ++i){	   
				if( pro.at(i).at(j).at(vv) < 10 ){ 
					instanceData << fixed << showpoint << setprecision(1) << pro.at(i).at(j).at(vv) << "   ";
				}else{
					instanceData << fixed << showpoint << setprecision(1) << pro.at(i).at(j).at(vv) << "  ";
				}
			}
		}
	}
	
	
	



     
    instanceData.close();
     /* --------------------------------------------------------------- 
	 * 					Supply of Vehicles
	 * --------------------------------------------------------------- */
    cout << "\tOK8" << flush;
    m = vec3Int( terminals, vec2Int( periods, vec1Int(vehicles , 0 ) ) );
    
    instanceData.open( sm.c_str(), ios::out );
    if( !instanceData )
    {
       cerr << "File could not be opened" << endl;
       return false;
    } 
    

    for ( int vv = 0; vv < vehicles; vv++  ){
		
		int i = 1 + rand()%terminals;
		instanceData << i << ";";
		int t = 1 + rand()%(periods==1?periods:periods-1);
		instanceData << t << ";";
		instanceData << vv+1 << ";";
		instanceData << vv+1 << endl;
		
		m.at(i-1).at(t-1).at(vv)++;
		
	}    
    for ( int vv = vehicles; vv < num_vehicles; vv++  ){
		
		int i = 1 + rand()%terminals;
		instanceData << i << ";";
		int t = 1 + rand()%(periods-1);
		instanceData << t << ";";
		instanceData << vv+1 << ";";
		int g = 1 + rand()%vehicles;
		instanceData << g << endl;
		
		m.at(i-1).at(t-1).at(g-1)++;
		
	}

	instanceData.close();
	
	/* --------------------------------------------------------------- 
	 * 					Demand for Vehicles
	 * --------------------------------------------------------------- */
    
	cout << "\tOK9" << flush;
	dem = vec3Int( terminals, vec2Int( terminals, vec1Int(periods , 0 ) ) );
    
    //cout << "\OK1,1" << flush; while( getchar() != '\n' );
    instanceData.open( sdem.c_str(), ios::out );
    if( !instanceData )
    {
       cerr << "File could not be opened" << endl;
       return false;
    } 
    //cout << "\OK1,2" << flush; while( getchar() != '\n' );
    for ( int n = 0; n < num_cargas; n++  ){
		int k;
		int i = 1 + rand()%terminals;
		instanceData << i << ";";
		do{
			k = 1 + rand()%terminals;
		}while ( i == k );
		instanceData << k << ";";
		int g = 1 + rand()%(periods==1?periods:periods-1);
		instanceData << g << ";";
		int h = 1 + rand()%10;
		instanceData << h << endl;
		
		//cout << "\n" << i << " " << k << " " << g << " " << h;
		
		dem.at(i-1).at(k-1).at(g-1) += h;
		
	}	
	instanceData.close();
	//cout << "\OK1,3" << flush; while( getchar() != '\n' );
	
     
    

   
  
     
   	/* --------------------------------------------------------------- 
	 * 					Data Summary
	 * --------------------------------------------------------------- */
    cout << "\tOK10" << flush;
    instanceData.open( ss.c_str(), ios::out );
    instanceData << terminals << "\t" << periods << "\t" << vehicles << "\n" << num_vehicles << "\t" << num_cargas << "\t" << p;
    instanceData.close();
    
    ofstream output ("Parameter_VAP_Generated.txt", ios_base::out);
    if( !output ){
    	cerr << "File of output could not be created" << endl;
    	return false;
    }    
    output << endl << endl;

    for(int x = 0; x < I; ++x){	   
	    output << endl;
		for(int y = 0; y < I; ++y){
			output << "tau["<< x+1 << "][" << y+1 << "]=" << fixed << showpoint << setprecision(1) << tau[x][y] << "  ";
		}
    }
    
    
 
    for (int z = 0; z < V; z++) {
	    output << endl;
		for(int y = 0; y < I; ++y){
		output << endl;
			for(int x = 0; x < I; ++x){	   
				output << "pro["<< x+1 << "][" << y+1 << "][" << z+1 << "]=" << fixed << showpoint << setprecision(1) << pro[x][y][z] << "  ";
			}
		}
    }
    output << endl << endl;

    for (int z = 0; z < V; z++) {
	    output << endl;
		for(int y = 0; y < I; ++y){
		output << endl;
			for(int x = 0; x < I; ++x){	   
				output << "cos["<< x+1 << "][" << y+1 << "][" << z+1 << "]=" << fixed << showpoint << setprecision(1) << cos[x][y][z] << "  ";
			}
		}
    }
    output << endl << endl;

    for (int z = 0; z < V; z++) {
	    output << endl;
		for(int y = 0; y < I; ++y){
		    //output << endl;
			for (int x = 0; x < T; x++){   
				output << "\ncosf["<< y+1 << "][" << x+1 << "][" << z+1 << "]=" << fixed << showpoint << setprecision(1) << cosf[y][x][z] << "  ";
			}
		}
    }
    output << endl << endl;

    for (int z = 0; z < T; z++){   
	    output << endl;
		for(int x = 0; x < I; ++x){
		//output << endl;
			for(int y = 0; y < J; ++y){
				output << "\nh["<< x+1 << "][" << y+1 << "][" << z+1 << "]=" << fixed << showpoint << setprecision(1) << In[x][y][z] << "  ";
			}
		}
    }
    output << endl << endl;

    for(int x = 0; x < I; ++x){	
	    //output << endl;
		for (int z = 0; z < T; z++){   
			output << "\nK["<< x+1 << "][" << z+1 << "]=" << fixed << showpoint << setprecision(1) << Kap[x][z] << setw(10) << flat_ixt(instance,z,x);
		}
    }
    output << endl << endl;

    for (int z = 0; z < V; z++) {
   	 //output << endl;
		for (int y = 0; y < I; y++){ 
		//output << endl;
			for (int x = 0; x < T; x++){ 
				if ( instance.m[y][x][z] != 0 ){
				output << "\nm["<< y << "][" << x << "][" << z << "]=" 
					<< setprecision(2) << instance.m[y][x][z] << setw(10) << flat_ixt(instance,x,y) << flush;
				}
			}
		}
     }
     output << endl << endl;

     for (int z = 0; z < T; z++) {
	 	//output << endl;
	 	for (int y = 0; y < I; y++){ 
	 	//output << endl;
	 		for (int x = 0; x < J; x++){ 
	 			if ( instance.dem[y][x][z] != 0 ){
	 				output << "\ndem["<< y << "][" << x << "][" << z << "]=" 
	 					<< instance.dem[y][x][z] << setw(10) << flat_ixt(instance,z,y) << setw(10) << flat_ixt(instance,(y==x?z+1:z + instance.tau.at(y).at(x)),x) << flush;
	 			}
	 		}
	 	}
      }
      output << endl << endl;
      cout << "\nOK18" << flush;
      for (int z = 0; z < V; z++) {
	  	output << endl;
	  	for (int y = 0; y < I; y++){ 
	  		output << endl;
	  		for (int x = 0; x < J; x++){ 
	  			if ( instance.A[y][x][z] == 0 ){
	  			output << "\nA["<< y << "][" << x << "][" << z << "]=" << 
	  				setprecision(2) << instance.A[y][x][z] << endl;
	  			}
	  		}
	  	}
	  }
      output << endl << endl;
    
    
    cout << endl << flush;
    
    
    

     
     cout << "\nBye from Generate_Instance\n" << flush;
     
     return true;
     
	
}


void tune_parameters_VAP( GRBModel &model ){
		// Solve the model with different values of Method
		int    bestMethod = -1;
		double bestTime = model.get(GRB_DoubleParam_TimeLimit);
		for (int i = 0; i <= 2; ++i) {
		  model.reset();
		  model.set(GRB_IntParam_Method, i);
		  model.optimize();
		  if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
			bestTime = model.get(GRB_DoubleAttr_Runtime);
			bestMethod = i;
			// Reduce the TimeLimit parameter to save time
			// with other methods
			model.set(GRB_DoubleParam_TimeLimit, bestTime);
		  }
		}

		// Report which method was fastest
		if (bestMethod == -1) {
		  cout << "Unable to solve this model" << endl;
		} else {
		  cout << "Solved in " << bestTime
			<< " seconds with Method: " << bestMethod << endl;
		}
	}



