/* figures.cc

	Non-interactive figure generation; added for the 2026 revision.

	Writes an image of the roots of the differences p_u - p_v, where p_u denotes
	the affine map u applied to the point z (a polynomial in z of degree |u|+1),
	u runs over words beginning with f and v over words beginning with g.
	See the subsection "Roots and R".

	Polynomials are vector<cpx> with the HIGHEST degree coefficient first,
	as in polynomial.cc.

	Usage:	wiggle -fig [maxlen] [rmax] [N] [outfile.ppm]
	e.g.	wiggle -fig 10 0.55 1000 roots_deg10.ppm
*/

polynomial poly_apply_f(polynomial P){
	// P -> z - zP		(prepend the letter f)
	polynomial Q=P;
	int k;
	Q[Q.size()-1] = Q[Q.size()-1] - 1.0;
	Q.push_back(0.0);
	for(k=0;k<(int)Q.size();k++){
		Q[k] = -1.0*Q[k];
	};
	return(Q);
};

polynomial poly_apply_g(polynomial P){
	// P -> 1 + (z-1)P		(prepend the letter g)
	polynomial Q=P;
	int k;
	Q.push_back(1.0);
	for(k=1;k<(int)Q.size();k++){
		Q[k] = Q[k] - P[k-1];
	};
	return(Q);
};

polynomial poly_diff(polynomial A, polynomial B){
	// A - B, aligned at the constant term, with leading zeros stripped
	int i,j,n;
	n = (A.size()>B.size() ? A.size() : B.size());
	polynomial C(n,cpx(0.0,0.0));
	for(i=0;i<(int)A.size();i++){
		C[n-A.size()+i] = C[n-A.size()+i] + A[i];
	};
	for(i=0;i<(int)B.size();i++){
		C[n-B.size()+i] = C[n-B.size()+i] - B[i];
	};
	j=0;
	while(j<(int)C.size()-1 && abs(C[j])<1.0e-14){
		j++;
	};
	C.erase(C.begin(),C.begin()+j);
	return(C);
};

cvec roots_dk(polynomial P){
	/* Durand-Kerner (Weierstrass) iteration, followed by a Newton polish
	   against the original polynomial.  Unlike roots() in polynomial.cc this
	   is deterministic: no random restarts, so figures are reproducible. */
	cvec R;
	polynomial M,D;
	int i,j,it,n;
	cpx d,step,dv;
	double delta;

	n = (int)P.size()-1;
	if(n<1){
		return(R);
	};
	M=P;
	for(i=0;i<=n;i++){
		M[i]=M[i]/P[0];			// make monic
	};
	cpx seed(0.4,0.9);			// classical spiral of initial guesses
	cpx w(1.0,0.0);
	for(i=0;i<n;i++){
		R.push_back(w);
		w=w*seed;
	};
	for(it=0;it<500;it++){
		delta=0.0;
		for(i=0;i<n;i++){
			d=cpx(1.0,0.0);
			for(j=0;j<n;j++){
				if(j!=i){
					d=d*(R[i]-R[j]);
				};
			};
			if(abs(d)<1.0e-300){
				continue;
			};
			step=eval(M,R[i])/d;
			R[i]=R[i]-step;
			if(abs(step)>delta){
				delta=abs(step);
			};
		};
		if(delta<1.0e-14){
			break;
		};
	};
	D=derivative(P);
	for(i=0;i<n;i++){
		for(j=0;j<3;j++){
			dv=eval(D,R[i]);
			if(abs(dv)>1.0e-300){
				R[i]=R[i]-eval(P,R[i])/dv;
			};
		};
	};
	return(R);
};

void write_roots_image(int maxlen, double rmax, int N, const char *filename){
	vector<polynomial> F,G,level,next;
	polynomial seed,A,B,C;
	cvec R;
	int n,i,j,k,px,py;
	long plotted,rejected_value,rejected_residual;
	double t,res,scale;
	cpx r;

	// generate p_u for all words u of length 1..maxlen, sorted by first letter
	seed.clear();
	seed.push_back(1.0);
	seed.push_back(0.0);		// p(z) = z
	level.clear();
	level.push_back(seed);
	for(n=1;n<=maxlen;n++){
		next.clear();
		for(i=0;i<(int)level.size();i++){
			A=poly_apply_f(level[i]);
			B=poly_apply_g(level[i]);
			F.push_back(A);
			G.push_back(B);
			next.push_back(A);
			next.push_back(B);
		};
		level=next;
	};
	cout << "generated " << F.size() << " words beginning with f and "
		<< G.size() << " beginning with g (lengths 1 to " << maxlen << ")\n";
	cout << "comparing " << (long)F.size()*(long)G.size() << " pairs\n";
	cout.flush();

	// white background
	vector<unsigned char> img(3*(long)N*(long)N,255);

	plotted=0;
	rejected_value=0;
	rejected_residual=0;
	scale=(double)N/(2.0*rmax);
	for(i=0;i<(int)F.size();i++){
		if(i%100==0){
			cout << "\r" << i << "/" << F.size();
			cout.flush();
		};
		A=F[i];
		for(j=0;j<(int)G.size();j++){
			B=G[j];
			C=poly_diff(A,B);
			if(C.size()<2){
				continue;		// identically zero, or a nonzero constant
			};
			R=roots_dk(C);
			for(k=0;k<(int)R.size();k++){
				r=R[k];
				if(abs(r-0.5)>=rmax){
					continue;
				};
				// discard anything that is not really a root
				res=abs(eval(C,r));
				if(res>1.0e-8*norm(C)){
					rejected_residual++;
					continue;
				};
				// the common value p_u(r) = p_v(r) must not be r itself
				if(abs(eval(A,r)-r)<1.0e-6){
					rejected_value++;
					continue;
				};
				px=(int)(scale*(r.real()-(0.5-rmax)));
				py=(int)(scale*(rmax-r.imag()));
				if(px>=0 && px<N && py>=0 && py<N){
					img[3*((long)py*(long)N+(long)px)+0]=0;
					img[3*((long)py*(long)N+(long)px)+1]=0;
					img[3*((long)py*(long)N+(long)px)+2]=0;
					plotted++;
				};
			};
		};
	};
	cout << "\rplotted " << plotted << " roots in |z-1/2| < " << rmax << "\n";
	cout << "discarded " << rejected_residual << " for large residual and "
		<< rejected_value << " with common value z\n";

	// the circle |z-1/2| = 1/2, in red
	for(k=0;k<40000;k++){
		t=TWOPI*(double)k/40000.0;
		r=cpx(0.5+0.5*cos(t),0.5*sin(t));
		px=(int)(scale*(r.real()-(0.5-rmax)));
		py=(int)(scale*(rmax-r.imag()));
		for(i=-1;i<=1;i++){
			for(j=-1;j<=1;j++){
				if(px+i>=0 && px+i<N && py+j>=0 && py+j<N){
					img[3*((long)(py+j)*(long)N+(long)(px+i))+0]=220;
					img[3*((long)(py+j)*(long)N+(long)(px+i))+1]=0;
					img[3*((long)(py+j)*(long)N+(long)(px+i))+2]=0;
				};
			};
		};
	};

	ofstream out(filename,ios::binary);
	out << "P6\n" << N << " " << N << "\n255\n";
	out.write((const char *)&img[0],img.size());
	out.close();
	cout << "wrote " << filename << " (" << N << "x" << N << ")\n";
};
