capture program drop classify
program classify, rclass
	version 14
	syntax [, var(string asis) PREdict]
	
	if "`predict'" == "predict" {
		mata: yname = st_global("e(depvar)")	
		mata: y = st_data(., yname)
		mata: yunique = uniqrows(y)
		mata: ycount = rows(yunique)

		mata: st_local("ycount", strofreal(ycount))
		mata: stata(`"predict p1-p`ycount'"')

		var = y p1-p`ycount
		
	}
	
	mata: classify("`var'", "`predict'" == "predict")
end