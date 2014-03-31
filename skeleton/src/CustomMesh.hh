/** User defined function for mesh modifications */

/** Example: 
		Read from PNG file */

conf_.registerIntParameter("readFromPNG");
bool configread = conf_.readFile("EnhancedFluidSimTest.par");
CHECK_MSG(configread, "Could not open file 'EnhancedFluidSimTest.par' which has to be in the current directory.");
if(specialProblems != 0){ // Added routine to save the resulting geometry created by the special geometry functions
	int Nx = grid_.Nx();
	int Ny = grid_.Ny();
	GrayScaleImage tempMesh(Nx,Ny);
	for (int j = 0; j < Ny; ++j)
		for (int i = 0; i < Nx; ++i)
			if (mesh_.geom_isBoundary(i+1,j+1)){
				const int yFlip = tempMesh.size(1) - j - 1;
				tempMesh.getElement(i,yFlip) = 0;
			}
	tempMesh.save("geom_fin.png");
}
if(conf_.getIntParameter("readFromPNG") == 1){
	CHECK_MSG( specialProblems == 0, "Error! You can either choose to readfromPNG or do a SpecProb, not both!");
	GrayScaleImage tempMesh("customMesh.png");
	GrayScaleImage newMesh;
	mesh_.geom_init(); // Cleans the geometry to just the boundary edges
	int Nx = grid_.Nx(); // Get grid width (internal)
	int Ny = grid_.Ny(); // Get grid height (internal)
	newMesh =	tempMesh.getResizedImage(Nx,Ny); // Resize the image to fit the mesh
	newMesh.save("geom_fin.png");
	for (int i = 1; i <= Nx; ++i)
		for (int j = 1; j <= Ny; ++j)
			if (newMesh(i-1,j-1) == 0)
				mesh_(i,j) = C_B;
			else
				mesh_(i,j) = C_F;
}

// Create an error
/*mesh_(20,20) = C_B;
mesh_(21,20) = C_B;
mesh_(20,21) = C_B;
mesh_(21,21) = C_B;
meshvis_ = mesh_;
meshvis_.geom_print();
mesh_.geom_finalize(); */

