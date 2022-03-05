capture program drop classify
program classify, rclass
	version 14
	syntax varlist(min=1) [, var(string asis) PREDict]
	
	if "`predict'" == "predict" {	
		mata: y = st_data(., "`varlist'")
		mata: yunique = uniqrows(y)
		mata: ycount = rows(yunique)

		mata: st_local("ycount", strofreal(ycount))
		mata: stata(`"predict p1-p`ycount'"')
		
	}
	
	mata: classify("`varlist'","`var'")
end