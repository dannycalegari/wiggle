/* polynomial.cc */

/* 
	polynomial P is a vector of cpx
   	degree d = P.size()-1;
   	P corresponds to the polynomial  z -> P[0]z^d + P[1]z^{d-1} + . .  + P[d]
*/

double pow(double base, int exponent){
	if(exponent==0){
		return(1.0);
	} else {
		return(base*pow(base,exponent-1));
	};
};

int degree(polynomial P){
	return(P.size()-1);
};

double norm(polynomial P){
	int i;
	double t;
	t=0.0;
	for(i=0;i<P.size();i++){
		t=t+abs(P[i]);
	};
	return(t);	
};

cpx eval(polynomial P, cpx z){
	cpx w;
	int i;
	
	w=0.0;
	for(i=0;i<P.size();i++){
		w=(z*w) + P[i];
		// w_0 = P[0], w_1 = zP[0]+P[1], w_2 = z^2P[0]+zP[1]+P[2] etc
	};
	
	return(w);
};

cpx eval_iterate(polynomial P, cpx z, int n){
	cpx w;
	int i;
	w=z;
	for(i=0;i<n;i++){
		w=eval(P,w);
	};
	
	return(w);
};

polynomial derivative(polynomial P){
	polynomial Q;
	int i;
	double t;
	Q.clear();
	for(i=0;i<P.size()-1;i++){
		t=(double) P.size()-i-1;
		Q.push_back(P[i]*t);
	};
	return(Q);
};

polynomial quotient(polynomial A, polynomial B){
	// returns the unique polynomial Q with A = QB+R and deg(R) < deg(B)
	polynomial Q,R;	
	Q.clear();	// initialize
	cpx w;
	int i;
	R=A;	// initialize
	
	if(B[0]==0.0){ // make sure B doesn't start with zero!
		cout << "leading coefficient of nonzero polynomial is zero! \n";
		cout.flush();
		assert(1==0);
	};
	
	
	while(R.size()>=B.size()){
		w=R[0]/B[0];
		Q.push_back(w);	// add R[0]/B[0] as coefficient to Q
		for(i=0;i<B.size();i++){	// subtract B*R[0]*z^{difference in degrees} from R
			R[i]=R[i]-B[i]*w;
		};
		// leading coefficient of R is zero; we need to remove it.
		R.erase(R.begin());
	};
	
	return(Q);
};

polynomial remainder(polynomial A, polynomial B){
	// returns the unique polynomial R with A = QB+R and deg(R) < deg(B)
	polynomial Q,R;	
	Q.clear();	// initialize
	cpx w;
	int i;
	R=A;	// initialize
	
	if(B[0]==0.0){ // make sure B doesn't start with zero!
		cout << "leading coefficient of nonzero polynomial is zero! \n";
		cout.flush();
		assert(1==0);
	};
	
	
	while(R.size()>=B.size()){
		w=R[0]/B[0];
		Q.push_back(w);	// add R[0]/B[0] as coefficient to Q
		for(i=0;i<B.size();i++){	// subtract B*R[0]*z^{difference in degrees} from R
			R[i]=R[i]-B[i]*w;
		};
		// leading coefficient of R is zero; we need to remove it.
		R.erase(R.begin());
	};
	
	return(R);
};

cvec roots(polynomial P){
	polynomial PP, Q, L;
	cvec R;
	cpx z;
	int step;
	bool verbose_mode;
	
	PP=P; // initialize;
	verbose_mode=false;
	
	while(PP[0]==0.0){
//		cout << "erasing initial 0 coefficient \n";
		PP.erase(PP.begin());
//		verbose_mode=true;
	};
	
	R.clear();	// find each root by Newton's method then iteratively divide
	step=0;
	while(PP.size()-1>2){
		if(verbose_mode){
			cout << "polynomial has degree " << PP.size()-1 << "\n";
		};
	// find next root
		z=7.227+2.141*I;
		Q=derivative(PP);
		step=0;

		while(abs(eval(PP,z))>0.00000000001){
			if(step>100 || abs(eval(Q,z))<0.00000000001){
				step=0;
				z=((double)rand()/(double)(RAND_MAX/3.0))*I+(double)rand()/(double)(RAND_MAX/3.0);
				if(verbose_mode){
					cout << "new seed\n";
					cout << "z= " << z << "\n";
				};
			};
			if(eval(Q,z)==0.0){
				cout << "division by zero!! \n";
			};
			z=z-(eval(PP,z)/eval(Q,z));	// what if eval(Q,z)=0?
			step++;
			if(verbose_mode){
				cout << "z " << z << " step " << step << "\n";
			};
		};
		if(verbose_mode){
			cout << "found an approximate root " << z << "\n";
		};
		R.push_back(z);
		// create L the linear polynomial z - z_0 where z_0 is the root we found
		L.clear();
		L.push_back(1.0);
		L.push_back(-z);	
		// divide PP by L and throw away the remainder
		PP=quotient(PP,L);
	};
	if(PP.size()-1==2){	// quadratic formula
		z=(-1.0*PP[1]+sqrt(PP[1]*PP[1]-4.0*PP[0]*PP[2]))/(2.0*PP[0]);
		R.push_back(z);
		z=(-1.0*PP[1]-sqrt(PP[1]*PP[1]-4.0*PP[0]*PP[2]))/(2.0*PP[0]);
		R.push_back(z);
	};
	if(PP.size()-1==1){	// linear formula
		z=-1.0*PP[1]/PP[0];
		R.push_back(z);
	};
	
	return(R);
};

cvec critical_points(polynomial P){
	return roots(derivative(P));
};

