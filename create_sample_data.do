* create sample data



set seed 1042
set obs 100

gen latitude = runiform()*2 - 1
gen longitude = runiform()*2 - 1

gen x = rnormal()
gen z = (rnormal() + x)/2
gen y = (rnormal() + x)/2

gen id = _n

save "sample_data.dta", replace
