/*	wiggle.cc

	May 11 2022 version 0.01

	Copyright Danny Calegari

	released under the terms of the GNU GPL
	
*/

// standard libraries to include

#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <vector>
#include <complex>
#include <sstream>
#include <limits>
#include <ctime>
#include <assert.h>
#include <unistd.h>
#include <algorithm>


using namespace std;

// preprocessor definitions

#define debug 	1							//	debug is true by default
#define verbose 1							//	verbose is true by default
#define PI 		3.1415926535897932384626433
#define TWOPI	6.283185307179586476925286766559
#define ACC		0.00000000001				// accuracy 1.0e-11
#define cpx		complex<double>				// complex number
#define cvec 	vector<cpx>					// vector of complex numbers
#define polynomial	vector<cpx>

cpx I (0.0,1.0);

#include "graphics.cc"
#include "polynomial.cc"
#include "roots.cc"
#include "figures.cc"
#include "prune.cc"
#include "certify.cc"
#include "interface.cc"

int main(int argc, char *argv[]){

	char c;
	std::cout << std::setprecision(12) << std::fixed;

	/* non-interactive figure generation; needs no X display.
	   wiggle -fig [maxlen] [rmax] [N] [outfile.ppm]			*/
	if(argc>1 && string(argv[1])=="-fig"){
		int maxlen = (argc>2 ? atoi(argv[2]) : 10);
		double rmax = (argc>3 ? atof(argv[3]) : 0.55);
		int N = (argc>4 ? atoi(argv[4]) : 1000);
		const char *out = (argc>5 ? argv[5] : "roots.ppm");
		write_roots_image(maxlen,rmax,N,out);
		return(0);
	};

	/* certified radius for a stable crossing at a single parameter.
	   wiggle -cert x y [depth]							*/
	if(argc>1 && string(argv[1])=="-cert"){
		double x = (argc>2 ? atof(argv[2]) : 0.3409);
		double y = (argc>3 ? atof(argv[3]) : 0.43486);
		int depth = (argc>4 ? atoi(argv[4]) : 200);
		long cap = (argc>5 ? atol(argv[5]) : 5000);
		report_certificate(cpx(x,y),depth,cap);
		return(0);
	};

	/* map of certificates.  wiggle -scan x y halfwidth N [depth] [out.ppm] */
	if(argc>1 && string(argv[1])=="-scan"){
		double x = (argc>2 ? atof(argv[2]) : 0.3409);
		double y = (argc>3 ? atof(argv[3]) : 0.43486);
		double hw = (argc>4 ? atof(argv[4]) : 0.0001);
		int N = (argc>5 ? atoi(argv[5]) : 200);
		int depth = (argc>6 ? atoi(argv[6]) : 60);
		const char *out = (argc>7 ? argv[7] : "scan.ppm");
		scan_certificates(cpx(x,y),hw,N,depth,out);
		return(0);
	};

	/* verify gamma(z) inside the envelope.  wiggle -env x y rho k */
	if(argc>1 && string(argv[1])=="-env"){
		double x = (argc>2 ? atof(argv[2]) : 0.3409);
		double y = (argc>3 ? atof(argv[3]) : 0.43486);
		double rho = (argc>4 ? atof(argv[4]) : 0.0001);
		int k = (argc>5 ? atoi(argv[5]) : 14);
		double worst;
		bool okk = verify_gamma_in_envelope(cpx(x,y),rho,k,worst);
		cout << (okk ? "VERIFIED" : "NOT VERIFIED")
			<< " : gamma(z) is contained in the envelope for all z in the disk"
			<< "\n  smallest clearance " << worst << "\n";
		return(0);
	};

	/* status map.  wiggle -locate x y halfwidth N [depth] [cap] [out.ppm] */
	if(argc>1 && string(argv[1])=="-locate"){
		double x = (argc>2 ? atof(argv[2]) : 0.3409);
		double y = (argc>3 ? atof(argv[3]) : 0.43486);
		double hw = (argc>4 ? atof(argv[4]) : 0.00005);
		int N = (argc>5 ? atoi(argv[5]) : 200);
		int depth = (argc>6 ? atoi(argv[6]) : 300);
		long cap = (argc>7 ? atol(argv[7]) : 20000);
		const char *out = (argc>8 ? argv[8] : "locate.ppm");
		locate_map(cpx(x,y),hw,N,depth,cap,out);
		return(0);
	};

	/* the W test on its own.  wiggle -w x y [depth] [cap] */
	if(argc>1 && string(argv[1])=="-w"){
		double x = (argc>2 ? atof(argv[2]) : 0.3409);
		double y = (argc>3 ? atof(argv[3]) : 0.43486);
		int depth = (argc>4 ? atoi(argv[4]) : 2000);
		long cap = (argc>5 ? atol(argv[5]) : 2000000);
		int level; long maxlist;
		int r = certify_W(cpx(x,y),depth,cap,level,maxlist);
		cout << "z = " << cpx(x,y) << "\n";
		cout << "  W test: " << (r==1 ? "LIST EMPTIED -> certified in int(W)"
			: (r==0 ? "still viable pairs at max depth -> not in int(W) by this test"
			: "list cap exceeded -> inconclusive"));
		if(r==1){ cout << " at level " << level; };
		cout << " ; max list " << maxlist << "\n";
		return(0);
	};

	/* shaded map.  wiggle -map x y halfwidth N [depth] [cap] [threads] [out.ppm] */
	if(argc>1 && string(argv[1])=="-map"){
		double x = (argc>2 ? atof(argv[2]) : 0.34090875);
		double y = (argc>3 ? atof(argv[3]) : 0.43484625);
		double hw = (argc>4 ? atof(argv[4]) : 0.00005);
		int N = (argc>5 ? atoi(argv[5]) : 400);
		int depth = (argc>6 ? atoi(argv[6]) : 300);
		long cap = (argc>7 ? atol(argv[7]) : 20000);
		int nt = (argc>8 ? atoi(argv[8]) : 8);
		const char *out = (argc>9 ? argv[9] : "map.ppm");
		map_depth(cpx(x,y),hw,N,depth,cap,nt,out);
		return(0);
	};

	/* certified loop.  wiggle -loop x y rad [depth] [coarse] [out.txt] */
	if(argc>1 && string(argv[1])=="-loop"){
		double x = (argc>2 ? atof(argv[2]) : 0.3409);
		double y = (argc>3 ? atof(argv[3]) : 0.43486);
		double rad = (argc>4 ? atof(argv[4]) : 0.00006);
		int depth = (argc>5 ? atoi(argv[5]) : 60);
		int coarse = (argc>6 ? atoi(argv[6]) : 2000);
		const char *out = (argc>7 ? argv[7] : "loop.txt");
		certified_loop(cpx(x,y),rad,depth,coarse,out);
		return(0);
	};

	setup_graphics();
	cout << "plot [r]oots or draw [m]ap? ";
	cin >> c;
	if(c=='r'){
		plot_roots(23,false);
	} else {
		user_interface();
	};
	cout << "quit? ";
	cin >> c;
	
	return(0);
};
