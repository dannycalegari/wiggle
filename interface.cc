/* interface.cc */
/*
   interactive mouse and keyboard driven environments
*/

void user_interface(){
	bool finished;
	finished=false;
	point p;
	double px,py;
	double xl,xr,yl,yr,xc,yc,w;
	int step,depth;
	
	step=16;
	depth=200;
	xc=0.3409;
	yc=0.43486;
	w=0.0001;
//	xc=0.5;
//	yc=0.0;
//	w=1.0;

	xl = xc - w/2.0;
	xr = xc + w/2.0;
	yl = yc - w/2.0;
	yr = yc + w/2.0;

	draw_prune_map(depth,step,xl,xr,yl,yr);

	
	while(finished==false){			
		XNextEvent(display, &report);
		switch(report.type) {
			case ButtonPress:	
				p=mouse_location();
				px = (double) p.x;
				py = (double) p.y;
				xc = xl + (px*w/800.0);
				yc = yr - (py*w/800.0);
				cout << "new center " << xc << " + " << yc << " i\n";
				w=w/2.0;
				xl = xc - w/2.0;
				xr = xc + w/2.0;
				yl = yc - w/2.0;
				yr = yc + w/2.0;
				erase_field();
				draw_prune_map(depth,step,xl,xr,yl,yr);
				break;
			case KeyPress:
			 //   erase_field();
			 	if(XLookupKeysym(&report.xkey, 0) == XK_q){ // quit  
                    finished=true;
                    XCloseDisplay(display);
                    cout << "Quitting program. \n";
                    exit(0);                
                    break;
                } else if(XLookupKeysym(&report.xkey, 0) == XK_z){ 
                	// refine
                	if(step>1){
                		step=step/2;
                	};
           			erase_field();
					draw_prune_map(depth,step,xl,xr,yl,yr);
					break;     	
                } else if(XLookupKeysym(&report.xkey, 0) == XK_x){ 
                	// coarsen
                	if(step<64){
                		step=step*2;
					};
           			erase_field();
					draw_prune_map(depth,step,xl,xr,yl,yr);
					break;     	
               } else if(XLookupKeysym(&report.xkey, 0) == XK_a){ 
                	// increase depth
                	depth++;
                	cout << "depth " << depth << "\n";
           			erase_field();
					draw_prune_map(depth,step,xl,xr,yl,yr);
					break;     	
                } else if(XLookupKeysym(&report.xkey, 0) == XK_s){ 
                	// decrease depth
                	if(depth > 10){
	                	depth--;
	                };
           			erase_field();
					draw_prune_map(depth,step,xl,xr,yl,yr);
					break;     	
                } 
                
            default:
            	break;
        };
    };
};
