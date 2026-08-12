/* prune.cc */

struct sim{
	cpx alpha,beta;
	// transformation w -> alpha w + beta
};

sim mul(sim A, sim B){
	sim C;
	C.alpha = A.alpha*B.alpha;
	C.beta = A.alpha*B.beta + A.beta;
	return(C);
};

sim inv(sim A){
	sim C;
	C.alpha = 1.0/A.alpha;
	C.beta = -1.0*A.beta/A.alpha;
	return(C);
};

sim make_sim(cpx alpha, cpx beta){
	sim A;
	A.alpha=alpha;
	A.beta=beta;
	return(A);
};

double maxradius(cpx z){
	double T1,T2;
	
	T1=abs(z/2.0 - 0.5)/abs(1.0-abs(z));
	T2=abs((z+1.0)/2.0 - 0.5)/abs(1.0-abs(1.0-z));
	if(T1>T2){
		return(T1);
	} else {
		return(T2);
	};
};

bool overlap(sim A, cpx z){
	// might cross; false means definitely doesn't cross
	double R;
	R=maxradius(z);
	if(R*(1.0+abs(A.alpha)) >= abs(A.beta+A.alpha/2.0 - 1.0/2.0)){
		return(true);
	} else {
		return(false);
	};
};

bool in_envelope(cpx z){
	// tests whether z is in the polygonal loop with vertices at envelope[i] for i=0,..,10
	
	cpx envelope[11] = {{0.35,0.5}, {0.2,0.3}, {0.05,0.3}, {-0.05,0.20}, {-0.05,-0.05},
		{0.4,0.03}, {0.55,-0.2}, {1.05,-0.1}, {1.05,0.05}, {0.8,0.15}, {0.8,0.35}};

	if(z.real()<-0.05 || z.real()>1.05 || z.imag()<-0.2 || z.imag()>0.5){
		return(false);
	} else {
		return(true);
	};
};


bool cross(sim A, cpx z){
	// true means definitely does cross
	/* note: this algorithm is only valid for 
		|z-(0.34095+0.434843)|<0.00005
	returns true iff the following 3 criteria are all satisfied:
	1. the line from 0 to 1 crosses the line from beta to alpha+beta;
	2. -beta/alpha and (1-beta)/alpha are not contained in envelope
	3. beta and alpha+beta are not contained in envelope
	*/

	cpx alpha,beta;
	alpha=A.alpha;
	beta=A.beta;
	double t,s;
	
	// step 1: does [0,1] cross [beta,alpha+beta]?
	
	if(beta.imag() * (alpha+beta).imag() < 0.0){	// beta, alpha+beta on opposite sides of x axis
		t=-beta.imag()/alpha.imag();
		s=beta.real()+t*alpha.real();
		if((0.0<s) && (s<1.0)){
			
		} else {
			return(false);
		};
	} else {
		return(false);
	};
	
	// step 2: are -beta/alpha and (1-beta)/alpha in envelope
	if(in_envelope(-1.0*beta/alpha)==false && in_envelope((1.0-beta)/alpha)==false && in_envelope(beta)==false && in_envelope(alpha+beta)==false){
		return(true);
	} else {
		return(false);
	};
};

int prune(cpx z, int depth){
	vector<sim> viable_list, new_list;
	sim A,AA,B,BB,C,D;
	int i,j,imax;
	
	viable_list.clear();
	new_list.clear();
	
	A=make_sim(-1.0*z,z);
//	cout << "A : " << A.alpha << " , " << A.beta << "\n";
//	AA=make_sim(-1.0/z,1.0);
	AA=inv(A);
	B=make_sim(z-1.0,1.0);
//	cout << "B : " << B.alpha << " , " << B.beta << "\n";
//	BB=make_sim(1.0/(z-1.0),-1.0/(z-1.0));
	BB=inv(B);
	C=make_sim((1.0-z)/z,(z-1.0)/z);	// seed
//	cout << "seed : " << C.alpha << " , " << C.beta << "\n";
//	cout << "R: " << maxradius(z) << "\n";
	
	// generate viable list to depth 1

	D=mul(mul(AA,C),B);	// 00 11; do this first because it *will* have overlap
	if(overlap(D,z)){
		viable_list.push_back(D);
	};	
	
	D=mul(mul(AA,C),A);	// 00 10
	if(overlap(D,z)){
		viable_list.push_back(D);
	};	

	D=mul(mul(BB,C),A);	// 01 10
	if(overlap(D,z)){
		viable_list.push_back(D);
	};
	D=mul(mul(BB,C),B);	// 01 11
	if(overlap(D,z)){
		viable_list.push_back(D);
	};
//	cout << "generated to depth 1; viable list is: \n";
//	for(i=0;i<viable_list.size();i++){
//		cout << i << ": " << viable_list[i].alpha << " , " << viable_list[i].beta << "\n";
//	};
	new_list=viable_list;
	viable_list.clear();
	
	// generate viable list to depth 2 *pruning special case*
	
	for(i=0;i<new_list.size();i++){
		C=new_list[i];
		if(i==0){
			D=mul(mul(AA,C),A);
			if(overlap(D,z)){
				viable_list.push_back(D);
			};
			D=mul(mul(AA,C),B);
			if(overlap(D,z)){
				viable_list.push_back(D);
			};
			
			// omit this one!	001 110
		//	D=mul(mul(BB,C),A);
		//	if(overlap(D,z)){
		//		viable_list.push_back(D);
		//	};
			
			D=mul(mul(BB,C),B);
			if(overlap(D,z)){
				viable_list.push_back(D);
			};	
		} else {
			D=mul(mul(AA,C),A);
			if(overlap(D,z)){
				viable_list.push_back(D);
			};
			D=mul(mul(AA,C),B);
			if(overlap(D,z)){
				viable_list.push_back(D);
			};
			D=mul(mul(BB,C),A);
			if(overlap(D,z)){
				viable_list.push_back(D);
			};	
			D=mul(mul(BB,C),B);
			if(overlap(D,z)){
				viable_list.push_back(D);
			};			
		};
	};
//	cout << "generated to depth 2; viable list is: \n";
//	for(i=0;i<viable_list.size();i++){
//		cout << i << ": " << viable_list[i].alpha << " , " << viable_list[i].beta << "\n";
//	};
	if(viable_list.size()==0){
		return(0);
	};
	imax=0;
	
	// loop
	for(j=2;j<depth;j++){
		new_list.clear();
		if(viable_list.size()==0){
			return(0);
		//	return(256*imax);
		};
		for(i=0;i<viable_list.size();i++){
			C=viable_list[i];
			if(abs(C.alpha)>1){	// subdivide right side
				D=mul(C,A);
				if(cross(D,z)){
		//			cout << C.alpha << " " << C.beta << " " << z << " cross!\n";
					return(viable_list.size());
				};
				if(overlap(D,z)){
					new_list.push_back(D);
				};
				D=mul(C,B);
				if(cross(D,z)){
					return(viable_list.size());
				};
				if(overlap(D,z)){
					new_list.push_back(D);
				};
			} else {	// subdivide left side
				D=mul(AA,C);
				if(cross(D,z)){
					return(viable_list.size());
				};
				if(overlap(D,z)){
					new_list.push_back(D);
				};
				D=mul(BB,C);
				if(cross(D,z)){
					return(viable_list.size());
				};
				if(overlap(D,z)){
					new_list.push_back(D);
				};			
			};
			if(new_list.size()>5000){
				return(5000);
			};				
		};
		viable_list=new_list;
		if(imax<viable_list.size()){
			imax=viable_list.size();
		};
	};
	
	// if we get this far there's still overlap
	
	return(viable_list.size());
};

void draw_prune_map(int depth, int step, double xl, double xr, double yl, double yr){
	double x,y,dstep;
	int XX,YY;
	dstep=(double) step;
	cpx z;
	int i,j;
	x=xl;
	for(XX=0;XX<=800;XX=XX+step){
		x=x+(xr-xl)*dstep/800.0;
		y=yl;
		for(YY=800;YY>=0;YY=YY-step){
			y=y+(yr-yl)*dstep/800.0;
			if((x-0.5)*(x-0.5)+(y*y)<=0.25){	// experiment!
				z=x+I*y;
				i=prune(z,depth);
				
				if(i>0){
					j=i/4;
					if(j>255){
						j=255;
					};
	//				cout << z << " " << i << "\n";
					fill_rectangle(XX,YY,step,step,j);
	//				draw_point(XX,YY,i/4);
	//				draw_point(5.0*(z.real()-0.5),-5.0*z.imag(),i);			
				};
			
			};	
		};
		XFlush(display);
	};
//	XFlush(display);
};