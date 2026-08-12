/* certify.cc

	Certified radii for stable crossings; added for the 2026 revision.

	Given z and a pair of words (u,v), let A = u^{-1}v be the similarity
	w -> alpha w + beta comparing gamma with A gamma.  The stable crossing test
	of prune.cc is the conjunction of eight strict inequalities on the four
	points

		P0 = A(0) = beta            Q0 = A^{-1}(0) = -beta/alpha
		P1 = A(1) = alpha+beta      Q1 = A^{-1}(1) = (1-beta)/alpha

	namely that Im(P0), Im(P1) have opposite signs, that Im(Q0), Im(Q1) have
	opposite signs (together these say the segment [0,1] crosses the segment
	[P0,P1]), and that all four points lie outside the envelope E.

	For each of the eight we compute a margin m_i at z, and a Lipschitz constant
	L_i valid on a disk about z, and certify the crossing on the disk of radius

		eps = min_i m_i / L_i .

	The Lipschitz constants are obtained by carrying d/dz along the composition
	and evaluating in disk ("ball") arithmetic over the disk, so they are genuine
	upper bounds for the derivative there, not derivatives at the centre.

	CAVEAT: the arithmetic below inflates radii by a fixed relative factor to
	absorb floating point error.  That is not a substitute for directed rounding;
	see the discussion in the paper.

	Usage:	wiggle -cert x y [depth]
*/

#define ROUND_FUDGE 1.0e-14		// relative inflation per ball operation

// envelope: the bounding box used by in_envelope() in prune.cc
#define ENV_XMIN (-0.05)
#define ENV_XMAX ( 1.05)
#define ENV_YMIN (-0.20)
#define ENV_YMAX ( 0.50)

/* The envelope.  prune.cc declares an 11 vertex polygon and then tests only
   its axis aligned bounding box; the polygon is the region Gamma of the paper.
   Since the polygon is contained in the box, "outside the polygon" is a weaker
   condition than "outside the box" and the clearance to the polygon is larger,
   so using the polygon certifies more pairs with bigger margins.  It does rely
   on the stronger containment gamma(z) subset polygon, which is checked by
   verify_gamma_in_envelope() below.				*/

#define ENV_N 11
cpx env_vertex(int i){
	static const double ex[ENV_N]={0.35,0.2,0.05,-0.05,-0.05,0.4,0.55,1.05,1.05,0.8,0.8};
	static const double ey[ENV_N]={0.5,0.3,0.3,0.20,-0.05,0.03,-0.2,-0.1,0.05,0.15,0.35};
	return(cpx(ex[i%ENV_N],ey[i%ENV_N]));
};

bool env_use_polygon = true;		// set false to recover prune.cc's box

double seg_distance(cpx w, cpx a, cpx b){
	cpx d;
	double t,n2;
	d=b-a;
	n2=norm(d);
	if(n2<=0.0){
		return(abs(w-a));
	};
	t=((w.real()-a.real())*d.real()+(w.imag()-a.imag())*d.imag())/n2;
	if(t<0.0){ t=0.0; };
	if(t>1.0){ t=1.0; };
	return(abs(w-(a+t*d)));
};

bool env_inside_polygon(cpx w){
	// ray casting along the positive real direction
	int i,cross;
	cpx a,b;
	double x;
	cross=0;
	for(i=0;i<ENV_N;i++){
		a=env_vertex(i);
		b=env_vertex(i+1);
		if((a.imag()>w.imag()) != (b.imag()>w.imag())){
			x = a.real() + (w.imag()-a.imag())*(b.real()-a.real())
				/(b.imag()-a.imag());
			if(x > w.real()){
				cross++;
			};
		};
	};
	return((cross%2)==1);
};

double env_boundary_distance(cpx w){
	int i;
	double d,t;
	d=1.0e300;
	for(i=0;i<ENV_N;i++){
		t=seg_distance(w,env_vertex(i),env_vertex(i+1));
		if(t<d){
			d=t;
		};
	};
	return(d);
};

double box_clearance(cpx w){
	// 0 if w is in the bounding box, otherwise the Euclidean distance to it
	double dx,dy;
	dx=0.0;
	dy=0.0;
	if(w.real()<-0.05){ dx=-0.05-w.real(); };
	if(w.real()> 1.05){ dx=w.real()-1.05; };
	if(w.imag()<-0.20){ dy=-0.20-w.imag(); };
	if(w.imag()> 0.50){ dy=w.imag()-0.50; };
	return(sqrt(dx*dx+dy*dy));
};

double dist_to_envelope(cpx w){
	// clearance of w from the envelope; 0 if w is inside it
	if(env_use_polygon==false){
		return(box_clearance(w));
	};
	if(env_inside_polygon(w)){
		return(0.0);
	};
	return(env_boundary_distance(w));
};

/* ---------------- disk arithmetic: { w : |w-c| <= r } ---------------- */

struct ball{
	cpx c;
	double r;
};

ball bmake(cpx c, double r){
	ball a;
	a.c=c;
	a.r=r;
	return(a);
};
ball bpt(cpx c){
	return(bmake(c,0.0));
};
double bsup(ball a){			// upper bound for |w|, w in a
	return(abs(a.c)+a.r);
};
double binf(ball a){			// lower bound for |w|, w in a
	double t;
	t=abs(a.c)-a.r;
	if(t<0.0){
		t=0.0;
	};
	return(t);
};
ball binflate(ball a){
	a.r = a.r*(1.0+ROUND_FUDGE) + ROUND_FUDGE*(abs(a.c)+1.0)*1.0e-2;
	return(a);
};
ball badd(ball a, ball b){
	return(binflate(bmake(a.c+b.c, a.r+b.r)));
};
ball bsub(ball a, ball b){
	return(binflate(bmake(a.c-b.c, a.r+b.r)));
};
ball bneg(ball a){
	return(bmake(-1.0*a.c, a.r));
};
ball bmul(ball a, ball b){
	return(binflate(bmake(a.c*b.c,
		abs(a.c)*b.r + abs(b.c)*a.r + a.r*b.r)));
};
ball binv(ball a){
	// requires a.r < |a.c|; returns a ball containing 1/w for all w in a
	double lo;
	lo=binf(a);
	if(lo<=0.0){
		return(bmake(cpx(0.0,0.0),1.0e300));	// no useful bound
	};
	return(binflate(bmake(1.0/a.c, a.r/(lo*abs(a.c)))));
};
ball bdiv(ball a, ball b){
	return(bmul(a,binv(b)));
};

/* ------- similarity carrying its derivative, in disk arithmetic ------- */

struct csim{
	ball alpha,beta;		// w -> alpha w + beta
	ball dalpha,dbeta;		// d/dz of the above
};

csim csim_f(ball z){		// f: w -> -zw + z
	csim A;
	A.alpha=bneg(z);
	A.beta=z;
	A.dalpha=bpt(cpx(-1.0,0.0));
	A.dbeta=bpt(cpx(1.0,0.0));
	return(A);
};
csim csim_g(ball z){		// g: w -> (z-1)w + 1
	csim A;
	A.alpha=bsub(z,bpt(cpx(1.0,0.0)));
	A.beta=bpt(cpx(1.0,0.0));
	A.dalpha=bpt(cpx(1.0,0.0));
	A.dbeta=bpt(cpx(0.0,0.0));
	return(A);
};
csim csim_id(){
	csim A;
	A.alpha=bpt(cpx(1.0,0.0));
	A.beta=bpt(cpx(0.0,0.0));
	A.dalpha=bpt(cpx(0.0,0.0));
	A.dbeta=bpt(cpx(0.0,0.0));
	return(A);
};
csim csim_mul(csim A, csim B){		// A o B
	csim C;
	C.alpha = bmul(A.alpha,B.alpha);
	C.beta  = badd(bmul(A.alpha,B.beta),A.beta);
	C.dalpha= badd(bmul(A.dalpha,B.alpha),bmul(A.alpha,B.dalpha));
	C.dbeta = badd(badd(bmul(A.dalpha,B.beta),bmul(A.alpha,B.dbeta)),A.dbeta);
	return(C);
};
csim csim_inv(csim A){
	csim C;
	ball ia,ia2;
	ia = binv(A.alpha);
	ia2 = bmul(ia,ia);
	C.alpha = ia;
	C.beta  = bneg(bmul(A.beta,ia));
	C.dalpha= bneg(bmul(A.dalpha,ia2));
	C.dbeta = badd(bneg(bmul(A.dbeta,ia)),bmul(bmul(A.beta,A.dalpha),ia2));
	return(C);
};
csim csim_word(string w, ball z){
	csim A;
	int i;
	A=csim_id();
	for(i=0;i<(int)w.size();i++){
		A = csim_mul(A, (w[i]=='f' ? csim_f(z) : csim_g(z)));
	};
	return(A);
};
csim csim_pair(string u, string v, ball z){		// A = u^{-1} v
	return(csim_mul(csim_inv(csim_word(u,z)),csim_word(v,z)));
};

/* ---------------------- the certified radius ---------------------- */

struct certificate{
	bool ok;
	string u,v;
	double eps;			// certified radius about z
	double worst_margin;		// the binding margin
	double worst_lipschitz;		// and the Lipschitz constant that bound it
	int binding;			// which of the eight conditions bound eps
};

certificate certified_radius(string u, string v, cpx z, double eps0){
	certificate cert;
	csim Ac,Ab;
	cpx P0,P1,Q0,Q1;
	double m[8],L[8],Lp0,Lp1,Lq0,Lq1,alo,bhi,chi,e,err;
	int i;

	cert.ok=false;
	cert.u=u;
	cert.v=v;
	cert.eps=0.0;
	cert.worst_margin=0.0;
	cert.binding=-1;

	// values at the centre
	Ac = csim_pair(u,v,bpt(z));
	if(binf(Ac.alpha)<=0.0){
		return(cert);
	};
	P0 = Ac.beta.c;
	P1 = Ac.alpha.c+Ac.beta.c;
	Q0 = -1.0*Ac.beta.c/Ac.alpha.c;
	Q1 = (1.0-Ac.beta.c)/Ac.alpha.c;
	err = Ac.alpha.r + Ac.beta.r;		// crude bound on the centre error

	// the two sign conditions must hold at the centre
	if(!(P0.imag()*P1.imag() < 0.0)){
		return(cert);
	};
	if(!(Q0.imag()*Q1.imag() < 0.0)){
		return(cert);
	};

	m[0]=fabs(P0.imag())-err;	m[1]=fabs(P1.imag())-err;
	m[2]=fabs(Q0.imag())-err;	m[3]=fabs(Q1.imag())-err;
	m[4]=dist_to_envelope(P0)-err;	m[5]=dist_to_envelope(P1)-err;
	m[6]=dist_to_envelope(Q0)-err;	m[7]=dist_to_envelope(Q1)-err;
	for(i=0;i<8;i++){
		if(m[i]<=0.0){
			return(cert);
		};
	};

	// Lipschitz constants, as sup bounds over the disk of radius eps0
	Ab = csim_pair(u,v,bmake(z,eps0));
	alo = binf(Ab.alpha);
	if(alo<=0.0){
		return(cert);
	};
	bhi = bsup(Ab.beta);
	chi = bsup(bsub(bpt(cpx(1.0,0.0)),Ab.beta));
	Lp0 = bsup(Ab.dbeta);
	Lp1 = bsup(Ab.dalpha)+bsup(Ab.dbeta);
	Lq0 = bsup(Ab.dbeta)/alo + bhi*bsup(Ab.dalpha)/(alo*alo);
	Lq1 = bsup(Ab.dbeta)/alo + chi*bsup(Ab.dalpha)/(alo*alo);
	L[0]=Lp0; L[1]=Lp1; L[2]=Lq0; L[3]=Lq1;
	L[4]=Lp0; L[5]=Lp1; L[6]=Lq0; L[7]=Lq1;

	e=eps0;
	for(i=0;i<8;i++){
		if(L[i]<=0.0){
			continue;		// this quantity is constant in z
		};
		if(m[i]/L[i] < e){
			e = m[i]/L[i];
			cert.binding=i;
			cert.worst_margin=m[i];
			cert.worst_lipschitz=L[i];
		};
	};
	cert.eps=e;
	cert.ok=true;
	return(cert);
};

/* The crossing test of prune.cc, in the symmetric determinant form and using
   whichever envelope is selected.  Equivalent to cross() when the box is used;
   with the polygon it accepts strictly more pairs, and so tends to certify with
   shorter words -- which matters because the Lipschitz constants, and hence the
   certified radius, degrade exponentially in the word length.	*/

bool cross_env(sim A, cpx z){
	cpx P0,P1,Q0,Q1;
	if(abs(A.alpha)<=0.0){
		return(false);
	};
	P0 = A.beta;
	P1 = A.alpha+A.beta;
	Q0 = -1.0*A.beta/A.alpha;
	Q1 = (1.0-A.beta)/A.alpha;
	if(!(P0.imag()*P1.imag() < 0.0)){
		return(false);
	};
	if(!(Q0.imag()*Q1.imag() < 0.0)){
		return(false);
	};
	if(dist_to_envelope(P0)<=0.0 || dist_to_envelope(P1)<=0.0){
		return(false);
	};
	if(dist_to_envelope(Q0)<=0.0 || dist_to_envelope(Q1)<=0.0){
		return(false);
	};
	return(true);
};

bool cross_rays(sim A, cpx z, int K);
struct certificate;
certificate certified_radius_rays(string u, string v, cpx z, double eps0);

/* ------------- search for a certifying pair, as in prune.cc ------------- */

struct wpair{
	string u,v;
	cpx alpha,beta;
};

cpx wp_alpha(cpx z, string w){
	cpx a;
	int i;
	a=cpx(1.0,0.0);
	for(i=0;i<(int)w.size();i++){
		a = a*(w[i]=='f' ? -1.0*z : z-1.0);
	};
	return(a);
};

sim sim_word(cpx z, string w){
	sim A;
	int i;
	A=make_sim(cpx(1.0,0.0),cpx(0.0,0.0));
	for(i=0;i<(int)w.size();i++){
		A = mul(A, (w[i]=='f' ? make_sim(-1.0*z,z) : make_sim(z-1.0,1.0)));
	};
	return(A);
};

wpair make_wpair(cpx z, string u, string v){
	wpair P;
	sim A;
	P.u=u;
	P.v=v;
	A = mul(inv(sim_word(z,u)),sim_word(z,v));
	P.alpha=A.alpha;
	P.beta=A.beta;
	return(P);
};

/* find_crossing returns
	 1 : a stable crossing was found (u_out,v_out); z is in R
	 0 : the viable list emptied; z is in the interior of W
	-1 : inconclusive (list exceeded the cap, or depth ran out)
   The old prune() conflated 0 with 1 and -1 in its return value, so its
   colour map could not distinguish a certificate from a failure.	*/
int find_crossing(cpx z, int depth, string &u_out, string &v_out, long cap,
		vector<wpair> *hits){
	/* the enumeration of prune.cc, carrying the words: four pairs at level 1,
	   four children each at level 2 with (ffg,ggf) omitted, then adaptive
	   two-way subdivision of the larger side. */
	vector<wpair> cur,next;
	wpair P;
	sim A;
	int j,i,k,window;
	string uu,vv;
	const char *L="fg";

	window=0;
	cur.clear();
	for(i=0;i<2;i++){
		for(j=0;j<2;j++){
			P = make_wpair(z,string("f")+L[i],string("g")+L[j]);
			A = make_sim(P.alpha,P.beta);
			if(overlap(A,z)){
				cur.push_back(P);
			};
		};
	};
	next.clear();				// one four-way level, as in prune.cc
	for(i=0;i<(int)cur.size();i++){
		for(j=0;j<4;j++){
			uu = cur[i].u + L[j/2];
			vv = cur[i].v + L[j%2];
			if(uu=="ffg" && vv=="ggf"){
				continue;		// the excluded pair
			};
			P = make_wpair(z,uu,vv);
			A = make_sim(P.alpha,P.beta);
			if(cross_rays(A,z,0)){
				hits->push_back(P);
			};
			if(overlap(A,z)){
				next.push_back(P);
			};
		};
	};
	cur=next;
	if(hits->size()>0){
		window=1;
	};
	for(j=2;j<depth;j++){
		next.clear();
		if(cur.size()==0){
			return(0);
		};
		for(i=0;i<(int)cur.size();i++){
			for(k=0;k<2;k++){
				if(abs(cur[i].alpha)>1.0){
					uu=cur[i].u;
					vv=cur[i].v+L[k];
				} else {
					uu=cur[i].u+L[k];
					vv=cur[i].v;
				};
				P = make_wpair(z,uu,vv);
				A = make_sim(P.alpha,P.beta);
				if(cross_rays(A,z,0)){
					hits->push_back(P);
				};
				if(overlap(A,z)){
					next.push_back(P);
				};
			};
			if((long)next.size()>cap){
				return(-1);
			};
		};
		cur=next;
		if(hits->size()>0){
			window++;
			if(window>=4 || hits->size()>4000){
				return(1);	// first productive level plus a few more
			};
		};
	};
	if(hits->size()>0){
		return(1);
	};
	return(-1);			// depth ran out
};

/* ------------------------ driver and self check ------------------------ */

/* Since eps ~ margin / L and L grows like |alpha|^{-1}, i.e. exponentially in
   the word length, the FIRST crossing found is usually not the best one.  We
   therefore collect every crossing available at the first level that has any,
   certify each, and keep the largest radius.			*/

double cross_min_margin(sim A, cpx z){
	cpx P0,P1,Q0,Q1;
	double m,t;
	if(abs(A.alpha)<=0.0){
		return(-1.0);
	};
	P0=A.beta; P1=A.alpha+A.beta;
	Q0=-1.0*A.beta/A.alpha; Q1=(1.0-A.beta)/A.alpha;
	if(!(P0.imag()*P1.imag()<0.0) || !(Q0.imag()*Q1.imag()<0.0)){
		return(-1.0);
	};
	m=fabs(P0.imag());
	t=fabs(P1.imag()); if(t<m){ m=t; };
	t=fabs(Q0.imag()); if(t<m){ m=t; };
	t=fabs(Q1.imag()); if(t<m){ m=t; };
	t=dist_to_envelope(P0); if(t<m){ m=t; };
	t=dist_to_envelope(P1); if(t<m){ m=t; };
	t=dist_to_envelope(Q0); if(t<m){ m=t; };
	t=dist_to_envelope(Q1); if(t<m){ m=t; };
	return(m);
};

certificate best_of_candidates(vector<wpair> &cand, cpx z, int howmany){
	certificate best,cert;
	vector<int> ord;
	vector<double> mar;
	sim A;
	int i,j,k,n;
	double e0;
	best.ok=false;
	best.eps=0.0;
	n=(int)cand.size();
	for(i=0;i<n;i++){
		A=make_sim(cand[i].alpha,cand[i].beta);
		mar.push_back(cross_min_margin(A,z)*abs(cand[i].alpha));
		ord.push_back(i);
	};
	// partial selection sort on the margin, largest first
	for(i=0;i<howmany && i<n;i++){
		k=i;
		for(j=i+1;j<n;j++){
			if(mar[ord[j]]>mar[ord[k]]){
				k=j;
			};
		};
		j=ord[i]; ord[i]=ord[k]; ord[k]=j;
		e0=1.0e-3;
		while(e0>1.0e-14){
			cert=certified_radius_rays(cand[ord[i]].u,cand[ord[i]].v,z,e0);
			if(cert.ok && cert.eps<e0){
				break;
			};
			e0=e0*0.5;
		};
		if(cert.ok && cert.eps>best.eps){
			best=cert;
		};
	};
	return(best);
};

certificate certify_point(cpx z, int depth, long cap, int &status){
	vector<wpair> hits;
	string u,v;
	certificate c;
	c.ok=false;
	c.eps=0.0;
	c.binding=-1;
	status = find_crossing(z,depth,u,v,cap,&hits);
	if(status!=1){
		return(c);
	};
	c = best_of_candidates(hits,z,300);
	if(c.ok==false){
		status=-1;
	};
	return(c);
};

bool crossing_holds_on_circle(string u, string v, cpx z, double r, int n){
	int i;
	double t;
	cpx zz;
	sim A;
	for(i=0;i<n;i++){
		t=TWOPI*(double)i/(double)n;
		zz = z + r*cpx(cos(t),sin(t));
		A = mul(inv(sim_word(zz,u)),sim_word(zz,v));
		if(cross_env(A,zz)==false){
			return(false);
		};
	};
	return(true);
};

double empirical_radius(string u, string v, cpx z){
	/* indicative only: the largest r (to within a factor 1.01) for which the
	   double precision test still holds all round the circle of radius r */
	double lo,hi,mid;
	int i;
	lo=0.0;
	hi=1.0e-8;
	for(i=0;i<60 && crossing_holds_on_circle(u,v,z,hi,400);i++){
		lo=hi;
		hi=hi*2.0;
	};
	for(i=0;i<40;i++){
		mid=0.5*(lo+hi);
		if(crossing_holds_on_circle(u,v,z,mid,400)){
			lo=mid;
		} else {
			hi=mid;
		};
	};
	return(lo);
};

void report_certificate(cpx z, int depth, long cap){
	string u,v;
	certificate cert;
	int i,n,bad;
	double t,rr,e0;
	cpx zz;
	sim A;
	const char *names[8]={"Im A(0)","Im A(1)","Im A^-1(0)","Im A^-1(1)",
		"env A(0)","env A(1)","env A^-1(0)","env A^-1(1)"};

	cout << "z = " << z << "\n";
	cert = certify_point(z,depth,cap,i);
	u=cert.u;
	v=cert.v;
	if(i!=1){
		cout << "  no stable crossing: status " << i
			<< (i==0 ? " (certified in the interior of W)"
				: " (inconclusive: list cap or depth)") << "\n";
		return;
	};
	cout << "  crossing found: u = " << u << " , v = " << v
		<< "  (lengths " << u.size() << "," << v.size() << ")\n";

	cout << "  certified radius eps = " << cert.eps << "\n";
	cout << "  binding condition   : "
		<< (cert.binding>=0 ? names[cert.binding] : "none (eps capped)")
		<< " , margin " << cert.worst_margin
		<< " , Lipschitz " << cert.worst_lipschitz << "\n";
	cout << "  empirical radius    : " << empirical_radius(u,v,z)
		<< "  (indicative; ratio to certified "
		<< empirical_radius(u,v,z)/cert.eps << ")\n";

	// self check: the double precision test should hold everywhere in the disk
	n=20000;
	bad=0;
	for(i=0;i<n;i++){
		t = TWOPI*(double)i/(double)n;
		rr = cert.eps*(0.999);
		zz = z + rr*cpx(cos(t),sin(t));
		A = mul(inv(sim_word(zz,u)),sim_word(zz,v));
		if(cross_env(A,zz)==false){
			bad++;
		};
	};
	cout << "  self check on the boundary circle: " << bad << "/" << n
		<< " failures\n";
};

/* ------------------- map of certificates and radii ------------------- */

void scan_certificates(cpx c, double hw, int N, int depth, const char *filename){
	/* Colour map on the square of half width hw about c:
		 white  : certified in the interior of W (viable list emptied)
		 grey   : inconclusive (list cap or depth exhausted)
		 blue->red by log10 of the certified radius : certified in R
	   This is the map Figure 7 should have been: the old prune() returned a
	   single int that could not distinguish these three outcomes.	*/
	vector<unsigned char> img(3*(long)N*(long)N,0);
	string u,v;
	certificate cert;
	cpx z;
	int i,j,st,ncert,nW,nbad;
	double x,y,t,e0;
	long k;

	ncert=0;
	nW=0;
	nbad=0;
	for(j=0;j<N;j++){
		for(i=0;i<N;i++){
			x = c.real() - hw + 2.0*hw*((double)i+0.5)/(double)N;
			y = c.imag() + hw - 2.0*hw*((double)j+0.5)/(double)N;
			z = cpx(x,y);
			k = 3*((long)j*(long)N+(long)i);
			cert = certify_point(z,depth,5000,st);
			if(st==0){
				nW++;
				img[k]=255; img[k+1]=255; img[k+2]=255;
			} else if(st<0){
				nbad++;
				img[k]=130; img[k+1]=130; img[k+2]=130;
			} else {
				if(cert.ok==false || cert.eps<=0.0){
					nbad++;
					img[k]=60; img[k+1]=60; img[k+2]=60;
				} else {
					ncert++;
					t = (log10(cert.eps)+10.0)/5.0;	// 1e-10 .. 1e-5
					if(t<0.0){ t=0.0; };
					if(t>1.0){ t=1.0; };
					img[k]  =(unsigned char)(255.0*t);
					img[k+1]=0;
					img[k+2]=(unsigned char)(255.0*(1.0-t));
				};
			};
		};
		cout << "\r" << j+1 << "/" << N;
		cout.flush();
	};
	cout << "\n  certified in R : " << ncert << "\n  certified in W : " << nW
		<< "\n  inconclusive   : " << nbad << "\n";
	ofstream out(filename,ios::binary);
	out << "P6\n" << N << " " << N << "\n255\n";
	out.write((const char *)&img[0],img.size());
	out.close();
	cout << "  wrote " << filename << "\n";
};

/* ------------------------- the certified loop -------------------------

   We cover the circle of radius rad about c by certified balls.  A ball
   D(z(theta_i),eps_i) covers the arc |theta-theta_i| < 2 asin(eps_i/(2 rad)),
   so each certificate certifies an interval of angles.  Rather than walking
   the circle step by step we

	phase 1: sample a coarse equally spaced set of angles, which normally
		 certifies most of the circle cheaply; then
	phase 2: repeatedly bisect whatever gaps remain, spending effort only
		 where phase 1 failed, and retrying resistant angles at
		 greater depth.

   Progress is reported as it goes, and certificates are appended to the
   output file as they are found, so a long run can be watched (and salvaged
   if it is interrupted).						*/

struct arc{
	double lo,hi;
};

bool arc_less(arc a, arc b){
	return(a.lo < b.lo);
};

double covered_measure(vector<arc> &v){
	// v must be sorted; returns total length of the union
	double t,hi;
	int i;
	t=0.0;
	hi=-1.0e300;
	for(i=0;i<(int)v.size();i++){
		if(v[i].lo > hi){
			t = t + (v[i].hi - v[i].lo);
			hi = v[i].hi;
		} else if(v[i].hi > hi){
			t = t + (v[i].hi - hi);
			hi = v[i].hi;
		};
	};
	return(t);
};

void find_gaps(vector<arc> &v, double lo, double hi, vector<arc> &gaps){
	// v sorted; collect the components of [lo,hi] not covered by the union
	double at;
	int i;
	arc g;
	gaps.clear();
	at=lo;
	for(i=0;i<(int)v.size();i++){
		if(v[i].lo > at){
			g.lo=at;
			g.hi=(v[i].lo<hi ? v[i].lo : hi);
			if(g.hi>g.lo){
				gaps.push_back(g);
			};
		};
		if(v[i].hi > at){
			at=v[i].hi;
		};
		if(at>=hi){
			return;
		};
	};
	if(at<hi){
		g.lo=at;
		g.hi=hi;
		gaps.push_back(g);
	};
};

struct loopcert{
	int status;			// as find_crossing
	certificate cert;
};

loopcert certify_angle(cpx c, double rad, double theta, int depth, long cap){
	loopcert L;
	string u,v;
	cpx z;
	double e0;
	z = c + rad*cpx(cos(theta),sin(theta));
	L.cert = certify_point(z,depth,cap,L.status);
	return(L);
};

void certified_loop(cpx c, double rad, int depth, int coarse, const char *filename){
	vector<arc> cov,gaps;
	arc a;
	loopcert L;
	cpx z;
	double theta,delta,frac,prev;
	int i,k,nc,nW,ninc,ndeep,round;
	long cap;
	long budget;

	nc=0; nW=0; ninc=0; ndeep=0;
	budget=0;
	cov.clear();

	ofstream out(filename);
	out << std::setprecision(17);
	out << "# certified balls covering the circle of radius " << rad
		<< " about " << c.real() << " " << c.imag() << "\n";
	out << "# columns: theta Re(z) Im(z) eps u v\n";
	out.flush();

	cout << "certifying the circle of radius " << rad << " about " << c << "\n";
	cout << "phase 1: " << coarse << " equally spaced angles at depth "
		<< depth << "\n";
	cout.flush();

	for(i=0;i<coarse;i++){
		theta = TWOPI*(double)i/(double)coarse;
		L = certify_angle(c,rad,theta,depth,5000);
		budget++;
		if(L.status==0){
			nW++;
		} else if(L.status<0){
			ninc++;
		};
		if(L.cert.ok && L.cert.eps>0.0){
			nc++;
			delta = 2.0*asin(L.cert.eps/(2.0*rad) < 1.0
				? L.cert.eps/(2.0*rad) : 1.0);
			a.lo=theta-delta;
			a.hi=theta+delta;
			cov.push_back(a);
			z = c + rad*cpx(cos(theta),sin(theta));
			out << theta << " " << z.real() << " " << z.imag() << " "
				<< L.cert.eps << " " << L.cert.u << " " << L.cert.v << "\n";
			out.flush();
		};
		if(i%(coarse/20 > 0 ? coarse/20 : 1)==0){
			sort(cov.begin(),cov.end(),arc_less);
			cout << "\r  phase 1: " << i+1 << "/" << coarse
				<< "  certified " << nc << "  W " << nW << "  inconclusive "
				<< ninc << "  coverage "
				<< 100.0*covered_measure(cov)/TWOPI << "%      ";
			cout.flush();
		};
	};
	sort(cov.begin(),cov.end(),arc_less);
	frac = covered_measure(cov)/TWOPI;
	cout << "\r  phase 1 done: " << nc << " certificates, coverage "
		<< 100.0*frac << "%                       \n";
	if(nW>0){
		cout << "  WARNING: " << nW << " angles are certified to lie in W;"
			<< " this radius does not work\n";
	};
	cout.flush();

	// phase 2: bisect the gaps, retrying resistant angles at greater depth
	for(round=0;round<60;round++){
		find_gaps(cov,0.0,TWOPI,gaps);
		if(gaps.size()==0){
			break;
		};
		prev = covered_measure(cov);
		cout << "  phase 2 round " << round+1 << ": " << gaps.size()
			<< " gaps, uncovered " << TWOPI-prev << " rad\n";
		cout.flush();
		for(k=0;k<(int)gaps.size();k++){
			if(gaps[k].hi-gaps[k].lo < 1.0e-13){
				continue;		// below resolution; give up on this one
			};
			theta = 0.5*(gaps[k].lo+gaps[k].hi);
			cap=5000;
			for(i=0;i<3;i++){		// spend more depth if resistant
				L = certify_angle(c,rad,theta,depth,cap);
				budget++;
				if(L.cert.ok && L.cert.eps>0.0){
					break;
				};
				if(L.status==0){
					break;		// really in W; more depth will not help
				};
				cap = cap*8;
				ndeep++;
			};
			if(L.cert.ok && L.cert.eps>0.0){
				nc++;
				delta = 2.0*asin(L.cert.eps/(2.0*rad) < 1.0
					? L.cert.eps/(2.0*rad) : 1.0);
				a.lo=theta-delta;
				a.hi=theta+delta;
				cov.push_back(a);
				z = c + rad*cpx(cos(theta),sin(theta));
				out << theta << " " << z.real() << " " << z.imag() << " "
					<< L.cert.eps << " " << L.cert.u << " " << L.cert.v << "\n";
				out.flush();
			} else if(L.status==0){
				cout << "    theta = " << theta << " is in W\n";
			};
			if(budget>2000000){
				cout << "    budget exhausted\n";
				break;
			};
		};
		sort(cov.begin(),cov.end(),arc_less);
		cout << "    coverage now " << 100.0*covered_measure(cov)/TWOPI
			<< "%  (" << nc << " certificates, " << ndeep
			<< " deep retries)\n";
		cout.flush();
		if(covered_measure(cov) <= prev + 1.0e-15 && gaps.size()>0){
			cout << "    no progress this round\n";
		};
		if(budget>2000000){
			break;
		};
	};

	find_gaps(cov,0.0,TWOPI,gaps);
	out.close();
	cout << "RESULT: " << nc << " certificates, coverage "
		<< 100.0*covered_measure(cov)/TWOPI << "%, "
		<< gaps.size() << " gaps remaining\n";
	if(gaps.size()==0){
		cout << "  the circle of radius " << rad << " about " << c
			<< " is covered by balls each contained in R,\n"
			<< "  hence lies in the interior of R and separates the island"
			<< " from infinity.\n";
	} else {
		for(k=0;k<(int)gaps.size() && k<10;k++){
			cout << "  gap " << k << " : [" << gaps[k].lo << ","
				<< gaps[k].hi << "] width " << gaps[k].hi-gaps[k].lo << "\n";
		};
	};
	cout << "  wrote " << filename << "\n";
};

void locate_map(cpx c, double hw, int N, int depth, long cap, const char *filename){
	/* status only, no certification: white = int(W), red = R, grey = unknown.
	   Much faster than scan_certificates; used to find the island.	*/
	vector<unsigned char> img(3*(long)N*(long)N,0);
	vector<wpair> hits;
	string u,v;
	cpx z;
	int i,j,st,nR,nW,nu;
	long k;
	nR=0; nW=0; nu=0;
	for(j=0;j<N;j++){
		for(i=0;i<N;i++){
			z = cpx(c.real()-hw+2.0*hw*((double)i+0.5)/(double)N,
				c.imag()+hw-2.0*hw*((double)j+0.5)/(double)N);
			k = 3*((long)j*(long)N+(long)i);
			hits.clear();
			st = find_crossing(z,depth,u,v,cap,&hits);
			if(st==0){
				nW++;
				img[k]=255; img[k+1]=255; img[k+2]=255;
			} else if(st==1){
				nR++;
				img[k]=0; img[k+1]=18; img[k+2]=139;
			} else {
				nu++;
				img[k]=120; img[k+1]=120; img[k+2]=120;
			};
		};
		cout << "\r  " << j+1 << "/" << N << "  W " << nW << "  R " << nR
			<< "  unknown " << nu;
		cout.flush();
	};
	cout << "\n";
	ofstream out(filename,ios::binary);
	out << "P6\n" << N << " " << N << "\n255\n";
	out.write((const char *)&img[0],img.size());
	out.close();
	cout << "  wrote " << filename << "\n";
};

/* ---------------- verifying gamma(z) inside the envelope ----------------

   Every certificate above assumes gamma(z) is contained in the envelope; this
   is NOT implied by Lemma 2.2, since the ball B of that lemma is much larger
   than the envelope (R is about 1.31 near the island).  We discharge the
   hypothesis directly.  Since gamma = union over |w|=k of w gamma and
   gamma subset B, we have

	gamma(z) subset union over |w|=k of wB,

   and wB is the disk of radius |alpha_w| R about w(1/2).  Evaluating w(1/2)
   and alpha_w in disk arithmetic over a disk of parameters, and bounding R
   above there, gives a finite check covering all z at once.		*/

bool verify_gamma_in_envelope(cpx z0, double rho, int k, double &worst){
	vector<csim> cur,next;
	ball z,c,half;
	double az,az1,den1,den2,R_hi,rad,cl;
	int j;
	long i;

	worst=1.0e300;
	az  = abs(z0)+rho;
	az1 = abs(z0-1.0)+rho;
	den1 = 1.0-az;			// lower bound for 1-|z|
	den2 = 1.0-az1;			// lower bound for 1-|1-z|
	if(den1<=0.0 || den2<=0.0){
		cout << "  the parameter disk is too big for Lemma 2.2\n";
		return(false);
	};
	R_hi = az1/(2.0*den1);
	if(az/(2.0*den2) > R_hi){
		R_hi = az/(2.0*den2);
	};
	cout << "  parameter disk |z-(" << z0.real() << "," << z0.imag()
		<< ")| <= " << rho << " ;  R <= " << R_hi << "\n";

	z = bmake(z0,rho);
	half = bpt(cpx(0.5,0.0));
	cur.clear();
	cur.push_back(csim_id());
	for(j=0;j<k;j++){
		next.clear();
		for(i=0;i<(long)cur.size();i++){
			next.push_back(csim_mul(cur[i],csim_f(z)));
			next.push_back(csim_mul(cur[i],csim_g(z)));
		};
		cur=next;
	};
	cout << "  checking " << cur.size() << " disks (words of length "
		<< k << ")\n";
	for(i=0;i<(long)cur.size();i++){
		c = badd(bmul(cur[i].alpha,half),cur[i].beta);	// w(1/2)
		rad = bsup(cur[i].alpha)*R_hi + c.r;		// radius of wB, inflated
		if(env_inside_polygon(c.c)==false){
			cout << "  FAIL: centre of wB is outside the envelope at word "
				<< i << " , " << c.c << "\n";
			worst=-1.0;
			return(false);
		};
		cl = env_boundary_distance(c.c) - rad;
		if(cl<worst){
			worst=cl;
		};
		if(cl<=0.0){
			cout << "  FAIL: disk sticks out at word " << i << " , centre "
				<< c.c << " radius " << rad << "\n";
			return(false);
		};
	};
	return(true);
};

/* ------------- map shaded by the depth needed to certify -------------

   Same classification as locate_map, but records for each parameter the total
   length |u|+|v| of the shortest certifying pair, and shades the blue by it.
   Multithreaded: the pixels are independent.				*/

#include <thread>
#include <atomic>

void map_worker(cpx c, double hw, int N, int depth, long cap, int t, int nt,
		vector<int> *stat, vector<int> *len, std::atomic<int> *done){
	vector<wpair> hits;
	string u,v;
	cpx z;
	int i,j,k,st,best;
	for(j=t;j<N;j=j+nt){
		for(i=0;i<N;i++){
			z = cpx(c.real()-hw+2.0*hw*((double)i+0.5)/(double)N,
				c.imag()+hw-2.0*hw*((double)j+0.5)/(double)N);
			hits.clear();
			st = find_crossing(z,depth,u,v,cap,&hits);
			(*stat)[(long)j*N+i]=st;
			best=0;
			if(st==1){
				best=1000000;
				for(k=0;k<(int)hits.size();k++){
					if((int)(hits[k].u.size()+hits[k].v.size())<best){
						best=(int)(hits[k].u.size()+hits[k].v.size());
					};
				};
			};
			(*len)[(long)j*N+i]=best;
		};
		(*done)++;
	};
};

void map_depth(cpx c, double hw, int N, int depth, long cap, int nt,
		const char *filename){
	vector<int> stat(( long)N*N,0), len((long)N*N,0);
	vector<std::thread> th;
	std::atomic<int> done(0);
	vector<unsigned char> img(3*(long)N*(long)N,0);
	int i,j,lo,hi,nW,nR,nU;
	long k;
	double t;

	cout << "mapping " << N << "x" << N << " at depth " << depth
		<< " on " << nt << " threads\n";
	cout.flush();
	for(i=0;i<nt;i++){
		th.push_back(std::thread(map_worker,c,hw,N,depth,cap,i,nt,&stat,&len,&done));
	};
	while(done < N){
		cout << "\r  rows " << done << "/" << N;
		cout.flush();
		usleep(500000);
	};
	for(i=0;i<nt;i++){
		th[i].join();
	};
	cout << "\r  rows " << N << "/" << N << "        \n";

	lo=1000000; hi=0; nW=0; nR=0; nU=0;
	for(k=0;k<(long)N*N;k++){
		if(stat[k]==1){
			nR++;
			if(len[k]<lo){ lo=len[k]; };
			if(len[k]>hi){ hi=len[k]; };
		} else if(stat[k]==0){
			nW++;
		} else {
			nU++;
		};
	};
	cout << "  W " << nW << "  R " << nR << "  undecided " << nU
		<< " ;  |u|+|v| ranges over " << lo << " .. " << hi << "\n";
	for(j=0;j<N;j++){
		for(i=0;i<N;i++){
			k=(long)j*N+i;
			if(stat[k]==0){
				img[3*k]=255; img[3*k+1]=255; img[3*k+2]=255;
			} else if(stat[k]<0){
				img[3*k]=255; img[3*k+1]=0; img[3*k+2]=255;	// flag
			} else {
				t = (hi>lo ? (double)(len[k]-lo)/(double)(hi-lo) : 0.0);
				t = 0.12+0.88*t;		// keep the easiest clearly not white
				img[3*k]  =(unsigned char)(255.0*(1.0-t)+0.0*t);
				img[3*k+1]=(unsigned char)(255.0*(1.0-t)+18.0*t);
				img[3*k+2]=(unsigned char)(255.0*(1.0-t)+139.0*t);
			};
		};
	};
	ofstream out(filename,ios::binary);
	out << "P6\n" << N << " " << N << "\n255\n";
	out.write((const char *)&img[0],img.size());
	out.close();
	// raw data, so the shading can be retuned without recomputing:
	// one short per pixel, -1 undecided, 0 in W, else |u|+|v|
	string dat = string(filename)+".dat";
	ofstream od(dat.c_str(),ios::binary);
	od << N << " " << N << "\n";
	for(k=0;k<(long)N*N;k++){
		short val = (stat[k]==0 ? 0 : (stat[k]<0 ? -1 : (short)len[k]));
		od.write((const char *)&val,sizeof(short));
	};
	od.close();
	cout << "  wrote " << filename << " and " << dat << "\n";
};

/* ------------- the W test, run on its own to high depth -------------

   find_crossing returns as soon as a crossing appears, so it can never report
   status 0 at a parameter where some pair also crosses.  To test membership of
   W independently we run the same enumeration with the crossing test switched
   off: if the list of viable pairs empties then f gamma and g gamma meet only
   at z, and z is in the interior of W.  Running both tests at the same
   parameter is a consistency check: they must never both succeed.	*/

int certify_W(cpx z, int depth, long cap, int &level, long &maxlist){
	vector<wpair> cur,next;
	wpair P;
	sim A;
	int j,i,k;
	string uu,vv;
	const char *L="fg";

	level=-1;
	maxlist=0;
	cur.clear();
	for(i=0;i<2;i++){
		for(j=0;j<2;j++){
			P = make_wpair(z,string("f")+L[i],string("g")+L[j]);
			A = make_sim(P.alpha,P.beta);
			if(overlap(A,z)){
				cur.push_back(P);
			};
		};
	};
	next.clear();
	for(i=0;i<(int)cur.size();i++){
		for(j=0;j<4;j++){
			uu = cur[i].u + L[j/2];
			vv = cur[i].v + L[j%2];
			if(uu=="ffg" && vv=="ggf"){
				continue;
			};
			P = make_wpair(z,uu,vv);
			A = make_sim(P.alpha,P.beta);
			if(overlap(A,z)){
				next.push_back(P);
			};
		};
	};
	cur=next;
	for(j=2;j<depth;j++){
		if((long)cur.size()>maxlist){
			maxlist=(long)cur.size();
		};
		if(cur.size()==0){
			level=j;
			return(1);			// certified in the interior of W
		};
		next.clear();
		for(i=0;i<(int)cur.size();i++){
			for(k=0;k<2;k++){
				if(abs(cur[i].alpha)>1.0){
					uu=cur[i].u;
					vv=cur[i].v+L[k];
				} else {
					uu=cur[i].u+L[k];
					vv=cur[i].v;
				};
				P = make_wpair(z,uu,vv);
				A = make_sim(P.alpha,P.beta);
				if(overlap(A,z)){
					next.push_back(P);
				};
			};
			if((long)next.size()>cap){
				return(-1);		// gave up
			};
		};
		cur=next;
	};
	return(0);					// still viable pairs at this depth
};

/* ================= the corrected stable crossing test =================

   The old cross() checked that the segment [0,1] meets the segment
   [A(0),A(1)] and that the four endpoints avoid the opposite envelope.  That
   is NOT sufficient: it does not ensure that four rays can be chosen whose
   ends link at infinity.  We instead exhibit the rays.

   A sufficient certificate: four straight rays from 0, A(0), 1, A(1) whose
   arguments are in alternating cyclic order, such that the rays from 0 and 1
   miss A Gamma, the rays from A(0) and A(1) miss Gamma, and each ray from one
   side is disjoint from each ray from the other.  Then the two proper paths
   meet only along gamma and A gamma, and Lemma 2.4 forces gamma to meet
   A gamma.								*/

double pt_seg_dist(cpx w, cpx a, cpx b){
	cpx d=b-a;
	double n2=norm(d),t;
	if(n2<=0.0){ return(abs(w-a)); };
	t=((w.real()-a.real())*d.real()+(w.imag()-a.imag())*d.imag())/n2;
	if(t<0.0){ t=0.0; };
	if(t>1.0){ t=1.0; };
	return(abs(w-(a+t*d)));
};
double pt_ray_dist(cpx w, cpx p, cpx d){		// d a unit vector
	double t=(w.real()-p.real())*d.real()+(w.imag()-p.imag())*d.imag();
	if(t<0.0){ t=0.0; };
	return(abs(w-(p+t*d)));
};
bool seg_ray_hit(cpx p, cpx d, cpx a, cpx b){
	cpx e=b-a;
	double den=d.real()*e.imag()-d.imag()*e.real();
	cpx w=a-p;
	double t,s;
	if(fabs(den)<1.0e-15){ return(false); };
	t=(w.real()*e.imag()-w.imag()*e.real())/den;		// along the ray
	s=(w.real()*d.imag()-w.imag()*d.real())/den;		// along the segment
	return(t>=0.0 && s>=0.0 && s<=1.0);
};
double ray_poly_dist(cpx p, cpx d, vector<cpx> &P){
	double best=1.0e300,t;
	int i,n=(int)P.size();
	for(i=0;i<n;i++){
		if(seg_ray_hit(p,d,P[i],P[(i+1)%n])){ return(0.0); };
	};
	for(i=0;i<n;i++){
		t=pt_ray_dist(P[i],p,d);
		if(t<best){ best=t; };
		t=pt_seg_dist(p,P[i],P[(i+1)%n]);
		if(t<best){ best=t; };
	};
	return(best);
};
double ray_ray_dist(cpx p, cpx d, cpx q, cpx e){
	double den=d.real()*e.imag()-d.imag()*e.real();
	cpx w=q-p;
	double t,s,b1,b2;
	if(fabs(den)>1.0e-15){
		t=(w.real()*e.imag()-w.imag()*e.real())/den;
		s=(w.real()*d.imag()-w.imag()*d.real())/den;
		if(t>=0.0 && s>=0.0){ return(0.0); };		// they meet
	};
	b1=pt_ray_dist(q,p,d);
	b2=pt_ray_dist(p,q,e);
	return(b1<b2 ? b1 : b2);
};

struct raycert{
	bool ok;
	cpx d[4];			// unit directions from 0, A(0), 1, A(1)
	double margin;			// the smallest of the eight clearances
};

vector<cpx> env_polygon(){
	vector<cpx> P;
	int i;
	for(i=0;i<ENV_N;i++){ P.push_back(env_vertex(i)); };
	return(P);
};

raycert find_rays(cpx alpha, cpx beta, int K){
	/* search directions on a grid of K for an alternating quadruple; the
	   radial choice from the crossing point of the two segments is tried
	   first, since it is automatically alternating and pairwise disjoint. */
	raycert R;
	static vector<cpx> G=env_polygon();
	vector<cpx> AG;
	cpx P[4],d[4],X;
	double t,m,best;
	int i,j,k,a,b,c,e;
	R.ok=false;
	R.margin=0.0;
	AG.reserve(ENV_N);
	for(i=0;i<ENV_N;i++){ AG.push_back(alpha*G[i]+beta); };
	P[0]=cpx(0.0,0.0); P[1]=beta; P[2]=cpx(1.0,0.0); P[3]=alpha+beta;

	// radial candidate
	if(beta.imag()*(alpha+beta).imag() < 0.0 && fabs(alpha.imag())>0.0){
		t = -beta.imag()/alpha.imag();
		X = cpx(beta.real()+t*alpha.real(),0.0);
		if(X.real()>0.0 && X.real()<1.0){
			best=1.0e300;
			for(i=0;i<4;i++){
				d[i]=P[i]-X;
				if(abs(d[i])<1.0e-14){ best=-1.0; break; };
				d[i]=d[i]/abs(d[i]);
			};
			if(best>0.0){
				for(i=0;i<4;i++){
					m = ray_poly_dist(P[i],d[i],(i%2==0 ? AG : G));
					if(m<best){ best=m; };
				};
				for(i=0;i<4;i+=2){
					for(j=1;j<4;j+=2){
						m = ray_ray_dist(P[i],d[i],P[j],d[j]);
						if(m<best){ best=m; };
					};
				};
				if(best>0.0){
					R.ok=true;
					R.margin=best;
					for(i=0;i<4;i++){ R.d[i]=d[i]; };
					return(R);
				};
			};
		};
	};
	if(K<=0){ return(R); };
	// coarse grid search, keeping the best margin
	vector<vector<int> > adm(4);
	vector<vector<double> > adm_m(4);
	for(i=0;i<4;i++){
		for(k=0;k<K;k++){
			cpx dd(cos(TWOPI*k/K),sin(TWOPI*k/K));
			m = ray_poly_dist(P[i],dd,(i%2==0 ? AG : G));
			if(m>0.0){ adm[i].push_back(k); adm_m[i].push_back(m); };
		};
		if(adm[i].size()==0){ return(R); };
	};
	best=0.0;
	for(a=0;a<(int)adm[0].size();a++){
	for(b=0;b<(int)adm[1].size();b++){
	for(c=0;c<(int)adm[2].size();c++){
	for(e=0;e<(int)adm[3].size();e++){
		int k0=adm[0][a],k1=adm[1][b],k2=adm[2][c],k3=adm[3][e];
		int b1=((k1-k0)%K+K)%K, b2=((k2-k0)%K+K)%K, b3=((k3-k0)%K+K)%K;
		if(!(0<b1 && b1<b2 && b2<b3)){ continue; };
		m = adm_m[0][a];
		if(adm_m[1][b]<m){ m=adm_m[1][b]; };
		if(adm_m[2][c]<m){ m=adm_m[2][c]; };
		if(adm_m[3][e]<m){ m=adm_m[3][e]; };
		if(m<=best){ continue; };
		cpx dd[4];
		dd[0]=cpx(cos(TWOPI*k0/K),sin(TWOPI*k0/K));
		dd[1]=cpx(cos(TWOPI*k1/K),sin(TWOPI*k1/K));
		dd[2]=cpx(cos(TWOPI*k2/K),sin(TWOPI*k2/K));
		dd[3]=cpx(cos(TWOPI*k3/K),sin(TWOPI*k3/K));
		double mm=m;
		for(i=0;i<4;i+=2){
			for(j=1;j<4;j+=2){
				double q=ray_ray_dist(P[i],dd[i],P[j],dd[j]);
				if(q<mm){ mm=q; };
			};
		};
		if(mm>best){
			best=mm; R.ok=true; R.margin=mm;
			for(i=0;i<4;i++){ R.d[i]=dd[i]; };
		};
	};};};};
	return(R);
};

bool cross_rays(sim A, cpx z, int K){
	/* The conditions of cross_env are all NECESSARY for the radial ray
	   construction -- the segments must cross for the radial centre to exist,
	   and each foot must be outside the opposite envelope for its ray to
	   escape -- so we use it as a cheap pre-filter.		*/
	raycert R;
	if(abs(A.alpha)<=0.0){ return(false); };
	if(cross_env(A,z)==false){ return(false); };
	R = find_rays(A.alpha,A.beta,K);
	return(R.ok);
};

certificate certified_radius_rays(string u, string v, cpx z, double eps0){
	/* Fix the four directions at their values at z.  Then the cyclic order
	   cannot change, and the certificate persists as long as the eight
	   clearances stay positive.  The rays from 0 and 1 are fixed sets and
	   only A Gamma moves; the rays from A(0), A(1) move with their bases. */
	certificate cert;
	csim Ac,Ab;
	raycert R;
	vector<cpx> G=env_polygon(),AG;
	cpx alpha,beta,P[4];
	double m[8],L[8],Lp,Lb0,Lb1,maxv,e;
	int i,j,n;

	cert.ok=false; cert.u=u; cert.v=v; cert.eps=0.0;
	cert.worst_margin=0.0; cert.worst_lipschitz=0.0; cert.binding=-1;

	Ac = csim_pair(u,v,bpt(z));
	if(binf(Ac.alpha)<=0.0){ return(cert); };
	alpha=Ac.alpha.c; beta=Ac.beta.c;
	R = find_rays(alpha,beta,72);
	if(R.ok==false){ return(cert); };

	AG.reserve(ENV_N);
	for(i=0;i<ENV_N;i++){ AG.push_back(alpha*G[i]+beta); };
	P[0]=cpx(0.0,0.0); P[1]=beta; P[2]=cpx(1.0,0.0); P[3]=alpha+beta;

	Ab = csim_pair(u,v,bmake(z,eps0));
	maxv=0.0;
	for(i=0;i<ENV_N;i++){ if(abs(G[i])>maxv){ maxv=abs(G[i]); }; };
	Lb0 = bsup(Ab.dbeta);
	Lb1 = bsup(Ab.dalpha)+bsup(Ab.dbeta);
	Lp  = maxv*bsup(Ab.dalpha)+bsup(Ab.dbeta);

	m[0]=ray_poly_dist(P[0],R.d[0],AG);	L[0]=Lp;
	m[1]=ray_poly_dist(P[2],R.d[2],AG);	L[1]=Lp;
	m[2]=ray_poly_dist(P[1],R.d[1],G);	L[2]=Lb0;
	m[3]=ray_poly_dist(P[3],R.d[3],G);	L[3]=Lb1;
	n=4;
	for(i=0;i<4;i+=2){
		for(j=1;j<4;j+=2){
			m[n]=ray_ray_dist(P[i],R.d[i],P[j],R.d[j]);
			L[n]=(j==1 ? Lb0 : Lb1);
			n++;
		};
	};
	for(i=0;i<8;i++){
		if(m[i]<=0.0){ return(cert); };
	};
	e=eps0;
	for(i=0;i<8;i++){
		if(L[i]<=0.0){ continue; };
		if(m[i]/L[i] < e){
			e=m[i]/L[i]; cert.binding=i;
			cert.worst_margin=m[i]; cert.worst_lipschitz=L[i];
		};
	};
	cert.eps=e; cert.ok=true;
	return(cert);
};
