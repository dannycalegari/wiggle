/* roots.cc */

void plot_roots(int degree, bool plot_all){

	vector<polynomial> P_list,PP_list,PP_new_list,left_list,right_list;
	polynomial P,Q,Q1,Q2,R;
	cvec root_list;
	int i,j,k,N;
	cpx z;
	double t;
	
	/* recursively generate finite list of polynomials
	polynomial is string of coefficients, constant term first
	
	(seed) 			0 1 is on the list
	(recursion)		if P is on the list then 
					z - zP is on the list; and
					1 + (z-1)P is on the list.
	*/
	//	P(z) = z
	P.clear();
	P.push_back(1.0);
	P.push_back(0.0);
	
	P_list.clear();
	PP_list.clear();
	P_list.push_back(P);
	PP_list.push_back(P);
		
	for(i=0;i<degree-1;i++){
		cout << i << " iterates \n";
		cout.flush();
		PP_new_list.clear();	// refresh new list each time
		for(j=0;j<PP_list.size();j++){
			Q=PP_list[j];
			Q1=Q;
			// add Q1 and Q2 to PP_new_list
			Q1[Q1.size()-1]=Q1[Q1.size()-1]-1.0;
			Q1.push_back(0.0);
			for(k=0;k<Q1.size();k++){
				Q1[k]=-1.0*Q1[k];
			};
			PP_new_list.push_back(Q1);
			Q2=Q;
			Q2.push_back(1.0);
			for(k=1;k<Q2.size();k++){
				Q2[k]=Q2[k]-Q[k-1];
			};
			PP_new_list.push_back(Q2);
		};
		PP_list=PP_new_list;
		// concatenate
		P_list.insert(P_list.end(),PP_list.begin(),PP_list.end());
	};
	
	left_list=PP_list;
	right_list=PP_list;
	for(j=0;j<left_list.size();j++){
		Q=left_list[j];
		Q1=Q;
		Q1[Q1.size()-1]=Q1[Q1.size()-1]-1.0;
		Q1.push_back(0.0);
		for(k=0;k<Q1.size();k++){
			Q1[k]=-1.0*Q1[k];
		};
		left_list[j]=Q1;
	};
	for(j=0;j<right_list.size();j++){
		Q=right_list[j];
		Q2=Q;
		Q2.push_back(1.0);
		for(k=1;k<Q2.size();k++){
			Q2[k]=Q2[k]-Q[k-1];
		};
		right_list[j]=Q2;	
	};
	
	
	setup_graphics();

	if(plot_all){
		cout << PP_list.size()-2 << " polynomials \n";
		for(i=0;i<left_list.size()-1;i++){
			P=left_list[i];
			cout << i << " ";
			cout.flush();

			for(j=0;j<right_list.size();j++){
				Q=right_list[j];
				if(P.size()!=Q.size()){
					cout << "error! comparing polynomials of different degree! \n";
				};
					
				R.clear();
				t=0.0;
				for(k=0;k<P.size();k++){
					R.push_back(P[k]-Q[k]);
					t=t+abs(P[k]-Q[k]);
				};
	//				cout << "polynomial " << i << "," << j << " has degree " << R.size()-1 << "\n";
				if(t>0.5){	// degenerate P=Q
					root_list=roots(R);
					/*
					
					*/
					for(k=0;k<root_list.size();k++){
						z=root_list[k];
			//			cout << z << " ";
						draw_point(1.5*(z.real()-0.5),-1.5*z.imag(),0x000000);
					};
				};
			};
		};	
	} else {	
		N=left_list.size();
		while(1){
			i=rand() % N;
			j=rand() % N;
			if(i!=j){
				P=left_list[i];
				Q=right_list[j];
				R.clear();
				for(k=0;k<P.size();k++){
					R.push_back(P[k]-Q[k]);
				};		
				root_list=roots(R);
				for(k=0;k<root_list.size();k++){
					z=root_list[k];
			//			cout << z << " ";
					draw_point(1.5*(z.real()-0.5),-1.5*z.imag(),0x000000);
				};		
			};
		};
	};

};