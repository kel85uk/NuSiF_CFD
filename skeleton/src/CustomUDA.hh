/** User defined function for accessing User Defined Array UDA inside the timeloop: 
	 Current Example: Scalar transport with convection-diffusion...
*/

/* Run to solve the inner points */
for(int i=1; i<= grid_.Nx(); ++i)
 for(int j=1; j<= grid_.Nx(); ++j)
    if( mesh_.geom_isFluid(i,j) )
 {
  real LAPLUDA = (grid_.UDA()(i+1,j)-2.0*grid_.UDA()(i,j)+grid_.UDA()(i-1,j))/(grid_.dx()*grid_.dx()) +
	  (grid_.UDA()(i,j+1)-2.0*grid_.UDA()(i,j)+grid_.UDA()(i,j-1))/(grid_.dy()*grid_.dy());
  real DUUDADX = ( (grid_.U()(i,j)*0.5*(grid_.UDA()(i,j)+grid_.UDA()(i+1,j)) -
			grid_.U()(i-1,j)*0.5*(grid_.UDA()(i-1,j)+grid_.UDA()(i,j))) +
	     gam_*(std::abs(grid_.U()(i,j))*0.5*(grid_.UDA()(i,j)-grid_.UDA()(i+1,j)) -
			std::abs(grid_.U()(i-1,j))*0.5*(grid_.UDA()(i-1,j)-grid_.UDA()(i,j)))
	  )/grid_.dx();
  real DVUDADY = ( (grid_.V()(i,j)*0.5*(grid_.UDA()(i,j)+grid_.UDA()(i,j+1)) -
			grid_.V()(i,j-1)*0.5*(grid_.UDA()(i,j-1)+grid_.UDA()(i,j))) +
	     gam_*(std::abs(grid_.V()(i,j))*0.5*(grid_.UDA()(i,j)-grid_.UDA()(i,j+1)) -
			std::abs(grid_.V()(i,j-1))*0.5*(grid_.UDA()(i,j-1)-grid_.UDA()(i,j)))
	  )/grid_.dy();
  grid_.UDA()(i,j) = grid_.UDA()(i,j)+dt_*(LAPLUDA/RE_/D_UDA - DUUDADX - DVUDADY);
 }
  
	  
/* Run to apply boundary conditions */
// Apply BC at the external boundaries
for(int j=0;j<=grid_.Ny()+1;++j){
	switch(WbctypeUDA_){
		case 1: {grid_.UDA()(0,j) = grid_.UDA()(1,j) - UDA_WBC_*grid_.dx(); break;}
		default: {grid_.UDA()(0,j) = 2*UDA_WBC_ - grid_.UDA()(1,j); break;}
	}
	switch(EbctypeUDA_){
		case 1: {grid_.UDA()(grid_.Nx()+1,j) = grid_.UDA()(grid_.Nx(),j) - UDA_EBC_*grid_.dx(); break;}
		default: {grid_.UDA()(grid_.Nx()+1,j) = 2*UDA_EBC_ - grid_.UDA()(grid_.Nx(),j); break;}	
	}	
}
for(int i=0;i<=grid_.Nx()+1;++i){
	switch(NbctypeUDA_){
		case 1: {grid_.UDA()(i,grid_.Ny()+1) = grid_.UDA()(i,grid_.Ny()) - UDA_NBC_*grid_.dy(); break;}
		default: {grid_.UDA()(i,grid_.Ny()+1) = 2*UDA_EBC_ - grid_.UDA()(i,grid_.Ny()); break;}	
	}
	switch(SbctypeUDA_){
		case 1: {grid_.UDA()(i,0) = grid_.UDA()(i,1) - UDA_NBC_*grid_.dy(); break;}
		default: {grid_.UDA()(i,0) = 2*UDA_EBC_ - grid_.UDA()(i,1); break;}	
	}
}

for(int i=1;i<=grid_.Nx();i++)
		for(int j=1;j<=grid_.Ny();j++)
			if(mesh_(i,j) >= B_N && mesh_(i,j) <= B_SE)
				switch (mesh_(i,j))
				{
					case B_N:  { 
									switch(ObctypeUDA_){
										case 1: {grid_.UDA()(i,j) = grid_.UDA()(i,j+1) - UDA_OBC_*grid_.dy();break;}
										default: {grid_.UDA()(i,j) = 2*UDA_OBC_ - grid_.UDA()(i,j+1);break;}
									}
									break;
						 }
					case B_E:  { 
										switch(ObctypeUDA_){
										case 1: {grid_.UDA()(i,j) = grid_.UDA()(i+1,j) - UDA_OBC_*grid_.dx();break;}
										default: {grid_.UDA()(i,j) = 2*UDA_OBC_ - grid_.UDA()(i+1,j);break;}
									}
									break;
						 }
					case B_S:  { 
									switch(ObctypeUDA_){
										case 1: {grid_.UDA()(i,j) = grid_.UDA()(i,j-1) - UDA_OBC_*grid_.dy();break;}
										default: {grid_.UDA()(i,j) = 2*UDA_OBC_ - grid_.UDA()(i,j-1);break;}
									}
									break;
						 }
					case B_W:  { 
									switch(ObctypeUDA_){
										case 1: {grid_.UDA()(i,j) = grid_.UDA()(i-1,j) - UDA_OBC_*grid_.dx();break;}
										default: {grid_.UDA()(i,j) = 2*UDA_OBC_ - grid_.UDA()(i-1,j);break;}
									}
									break;
						 }
					case B_NE: { 
									switch(ObctypeUDA_){
										case 1: {grid_.UDA()(i,j) = grid_.UDA()(i,j+1) - UDA_OBC_*grid_.dy();break;}
										default: {grid_.UDA()(i,j) = 2*UDA_OBC_ - grid_.UDA()(i,j+1);break;}
									}
									break;
						 }
					case B_SE: { 
									switch(ObctypeUDA_){
										case 1: {grid_.UDA()(i,j) = grid_.UDA()(i+1,j) - UDA_OBC_*grid_.dx();break;}
										default: {grid_.UDA()(i,j) = 2*UDA_OBC_ - grid_.UDA()(i+1,j);break;}
									}
									break;
								}
					case B_SW: { 
									switch(ObctypeUDA_){
										case 1: {grid_.UDA()(i,j) = grid_.UDA()(i,j-1) - UDA_OBC_*grid_.dy();break;}
										default: {grid_.UDA()(i,j) = 2*UDA_OBC_ - grid_.UDA()(i,j-1);break;}
									}
									break;
						 }
					case B_NW: { 
									switch(ObctypeUDA_){
										case 1: {grid_.UDA()(i,j) = grid_.UDA()(i-1,j) - UDA_OBC_*grid_.dx();break;}
										default: {grid_.UDA()(i,j) = 2*UDA_OBC_ - grid_.UDA()(i-1,j);break;}
									}
									break;
						 }
					default : break;
				}
