*
* Eric Lewis
*  February 2018
*
* Program that reads in stored matrix results as a .txt file, 
*  loads into Mata as a matrix called distmat
*  
*  Argument is the name of the distance matrix 
*	 Should not include the .txt extension
*	 Can include the file path

program define read_distance_matrix
	version 14
	syntax anything

	mata {
		fh = _fopen("`anything'.txt","rw")
		distmat = fgetmatrix(fh)
		fclose(fh)
	}

end


	
	
