*  Eric Lewis
*  February 2018
*
* A program that creates the distance matrix text file
*
*  The argument `anything' is the name of the file that the matrix will be saved to
*     It can include a file path.
*     It should not include an extension (.txt is automatically added by the code)
* 
*  Latitude is the latitude in decimal degrees
*  Longitude is the longitude in decimal degress
*  Indexvar is a variable that indexes the rows of the distance matrix, allowing
*    later Stata code to match observations with rows/columns in the matrix using
*    the index variable
*


program create_distance_matrix
	version 14
	syntax anything [if] [in] , LATitude(varname numeric) LONgitude(varname numeric) ///
		INDEXvar(varname numeric)

	* drops certain observations
	preserve
	if "`if'"!="" qui keep `if'
	if "`in'"!="" qui keep `in'

	* sorts by index variable
	sort `indexvar'

	* makes sure that indexvar will serve as an index to matrix:
	qui count if missing(`indexvar')
	if r(N)>0 {
		di as error "indexvar `indexvar' contains missing observations"
		exit
	}
	qui count if `indexvar'[1]!=1
	if r(N)>0 {
		list `indexvar' in 1
		di as error "minimum value of indexvar `indexvar' must be 1"
		exit
	}
	qui count if `indexvar'[_N]!=_N
	if r(N)>0 {
		list `indexvar' in _N
		di as error "maximum value of indexvar `indexvar' must be equal to total number of observations"
		exit
	}
	qui count if `indexvar' - `indexvar'[_n-1] != 1 & _n>1
	if r(N)>0 {
		sum `indexvar'
		di as error "indexvar `indexvar' must serve as an index."
		di as reror "`indexvar' must uniquely identify observations, "
		di as error "take integer values, and range from 1 to the number of observations"
		exit
	}
	
	
	mata {
		// reads data into mata
		lat = st_data(.,"`latitude'")
		lon = st_data(.,"`longitude'")
		
		// counts total observations
		n = rows(lat)
		
		// creates blank distance matrix that will get filled in:
		distmat = J(n,n,0) 

		// sets up parameters for how quickly to loop through observations
		incrementsize = 1000
		loop_ceil = ceil(n/incrementsize)
		
		// loops through blocks of rows of distmat to fill them in:
		for (i = 1; i <= loop_ceil; i++){
			lower_i = incrementsize*i-incrementsize+1
			upper_i = min((incrementsize*i, n))

			n_i = upper_i-lower_i+1

			// loops through blocks of columns of distmat to fill them in:
			for (j = i; j<= loop_ceil; j++){
				"on interation (i="+strofreal(i)+", j="+strofreal(j)+")"

				lower_j = incrementsize*j-incrementsize+1
				upper_j = min((incrementsize*j, n))
	
				n_j = upper_j-lower_j+1
			
				latmat_i = J(1, n_j, lat[lower_i..upper_i])*pi()/180
				latmat_j = (J(1, n_i, lat[lower_j..upper_j]))'*pi()/180
				lonmat_i = J(1, n_j, lon[lower_i..upper_i])*pi()/180
				lonmat_j = (J(1, n_i, lon[lower_j..upper_j]))'*pi()/180
	
				latdif_ij = latmat_i - latmat_j
				londif_ij = lonmat_i - lonmat_j
				
				a_ij = sin(latdif_ij/2):^2+cos(latmat_i):*cos(latmat_j):*(sin(londif_ij/2):^2)
				distmat_ij = 2*6371*asin(a_ij:^.5)	// 6371 turns the distance into km at earth's surface

				// fills in part of the distmat matrix
				distmat[lower_i..upper_i,lower_j..upper_j] = distmat_ij
				if (j>i) {
					distmat[lower_j..upper_j,lower_i..upper_i] = distmat_ij'
				}
			}
		}

		// saves the matrix as a txt file
		fh = fopen("`anything'.txt","rw")
		fputmatrix(fh, distmat)
		fclose(fh)
		
	}
	// end of Mata commands

	mata: mata describe distmat
	
	* restores data to unsorted, and also restores any observations dropped by if or in
	restore
	
end	
	
	
	
