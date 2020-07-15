#!/bin/awk -f

BEGIN {
	divergent_trans = 0
	low_bfmi = 0
	time = "NA"
	max_rhat = "NA"
	min_ess_bulk = "NA"
	min_ess_tail = "NA"

	num_countries = 0
 	num_total_subregions = 0

	country = "NA"
	country_code = "NA"

	iter[1] = 0
	iter[2] = 0
	iter[3] = 0
	iter[4]	= 0

	run_status = "GOOD"

	# print "job_id\tcountry_index\tdivergent_trans\tmax_rhat\tmin_ess_bulk\tmin_ess_tail" 
}

tolower($0) ~ /error/ {
	run_status = "FAIL"
}

match($0, /Running model with ([[:digit:]]+) countries and ([[:digit:]])+ total subnational entities/, run_matches) {
	num_countries = strtonum(run_matches[1])
	num_total_subregions = strtonum(run_matches[2])

	if (num_countries > 1) {
		country = "multi"
	}
}

match($0, /(\w[^\[]+) \[(\w+)\]: [[:digit:]]+ sub regions./, country_match) {
	if (num_countries == 1) {
		country = country_match[1] 
		country_code = country_match[2]
	}
}

match($0, /Chain ([[:digit:]]): Iteration:[^\[]+\[\s*([[:digit:]]+)/, iter_matches) {
	iter[strtonum(iter_matches[1])] = iter_matches[2]	
}

/Sampling: [[:digit:]]+(\.[[:digit:]]+)? sec elapsed/ {
	time = $2
	run_status = "DONE"
}

/There were [[:digit:]]+ divergent/ {
	divergent_trans = $4
}

/There were [[:digit:]]+ chains where the estimated Bayesian Fraction of Missing Information was low/ {
	low_bfmi = $4
}

/Maximum Rhat/ { 
	max_rhat = $4
}
/Minimum ESS Bulk/ {
	min_ess_bulk = $5
}

/Minimum ESS Tail/ {
	min_ess_tail = $5
}

END {
	match(FILENAME, /([[:digit:]]+)(_([[:digit:]]+))?\.log$/, m)

	print m[1], "\t", m[3], "\t", country_code, "\t", divergent_trans, "\t", low_bfmi, "\t", max_rhat, "\t", min_ess_bulk, "\t", min_ess_tail, "\t", iter[1] , "\t", iter[2], "\t", iter[3], "\t", iter[4] "\t", run_status, "\t", time , "\t", country 
}
