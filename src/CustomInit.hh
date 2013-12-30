/** User defined function for initial condition modifications */

switch(specialProblems){
	case 1: // Just modify the cells that were initialized by default
		for (int j = 1; j <= grid_.Ny(); ++j)
			for (int i = 1; i <= grid_.Nx(); ++i){
				if(mesh_.geom_isFluid(i,j) && grid_.YC()(i,j) <= conf_.getRealParameter("RectangleY2")){
					grid_.U()(i,j) = 0.;
					grid_.V()(i,j) = 0.;
					grid_.p()(i,j) = 1.;
				}
			}
		break;
	default: break;
}
