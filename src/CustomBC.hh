/** User defined function for boundary condition modifications */
/** Example parabolic profile for Poiseuille flow
for (int j = 1; j <= ny; ++j)
	U_WBCv(j) = 2*U_WBC*(1-(grid_.YC()(1,j)-0.5)*(grid_.YC()(1,j)-0.5)/(0.5*0.5)); */

switch(specialProblems){
	case 1:
		for (int j = 1; j <= ny; ++j)
			if(grid_.YC()(1,j) >= conf_.getRealParameter("RectangleY2"))
				U_WBCv(j) = 1;
			else
				U_WBCv(j) = 0.;
		break;
	default: break;
}
