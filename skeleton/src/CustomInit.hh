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

if(customUDA){
	grid_.UDA().fill(conf_.getRealParameter("UDA_init"));
	D_UDA = conf_.getRealParameter("D_UDA");
	NbctypeUDA_ = conf_.getIntParameter("boundary_conditionUDA_N");
	EbctypeUDA_ = conf_.getIntParameter("boundary_conditionUDA_E");
	SbctypeUDA_ = conf_.getIntParameter("boundary_conditionUDA_S");
	WbctypeUDA_ = conf_.getIntParameter("boundary_conditionUDA_W");
	ObctypeUDA_ = conf_.getIntParameter("boundary_conditionUDA_O");
	UDA_NBC_ = conf_.getRealParameter("boundary_UDA_N");
	UDA_EBC_ = conf_.getRealParameter("boundary_UDA_E");
	UDA_SBC_ = conf_.getRealParameter("boundary_UDA_S");
	UDA_WBC_ = conf_.getRealParameter("boundary_UDA_W");
	UDA_OBC_ = conf_.getRealParameter("boundary_UDA_O");
}
