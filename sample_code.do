*
* Sample code that:
*    (1) Creates a fake data set (create_sample_data.do)
*    (2) Creates and saves a distance matrix using create_distance_matrix
*    (3) Reads in the distance matrix using read_distance_matrix
*    (4) Demonstrates the use of the spatial regression with and without the
*         use of the pre-saved distance matrix



// change directory if needed
// cd <your directory>

// reads in programs. 
// Alternatively, store these programs where sysdir keeps your personal .ado files:
do create_distance_matrix.ado
do read_distance_matrix.ado
do spatial_reg.ado

// creates sample data:
clear
do create_sample_data


// creates the distance matrix (called newdistmat.txt) using Mata
// saves it as newdistmat.txt
set trace on
create_distance_matrix newdistmat, latitude(latitude) longitude(longitude) indexvar(id)
set trace off

// Creates a variable that is a contant -- unlike Stata commands from Stata Corp, my code 
// does not automatically create a _cons variable
gen constant = 1

// check whether the saved distance matrix called is indeed in Mata memory (called distmat if in Mata)
capture mata: mata describe distmat
if _rc!=0 di "distmat is not in Mata memory"
if _rc==0 {
	mata: mata describe distmat
	// it's good to check whether distmat has the right dimension, especially
	// if your work uses multiple different distance matrices
}


// Reads the distance matrix (newdistmat.txt) into Mata preparatory to running regression
// Saves it in Mata as distmat
// This is not necessary if the distmat matrix has not been cleared from Mata memory
read_distance_matrix newdistmat


// OLS:
	* without spatial standard errors
reg y x

	* calling program where distmat is used and id indexes observations:
spatial_reg y x constant, indexvar(id) distcutoff(50) small diprog

	* calling program where distances are calculated within regression rather than using saved distmat
spatial_reg y x constant, lat(latitude) lon(longitude) distcutoff(50) small diprog

	* bartlett weights:
spatial_reg y x constant, lat(latitude) lon(longitude) distcutoff(50) small diprog bartlett

	* it's okay if the dimensionality of distmat is larger than number of observations
	* the important thing is that indexvar indexes the observations within the Mata matrix distmat
spatial_reg y x constant if y>0, indexvar(id) distcutoff(50) small diprog


// 2SLS:
	* without spatial standard errors
ivregress 2sls y (x = z)

	* with spatial standard errors:
spatial_reg y constant (x = z), indexvar(id) distcutoff(50) small diprog



// Estimates can be stored using estimates store and then displayed in a table:
reg y x constant, noconstant
estimates store ols

spatial_reg y x constant, indexvar(id) distcutoff(50) small diprog
estimates store spat_ols

esttab ols spat_ols

// many standard post-estimation parameters are stored
estimates restore spat_ols
ereturn list









