capture program drop classify
program classify, rclass
	version 14
	capture syntax varlist(min=1 max=2) [, PROBS(varlist) PREDict POWER(real 1.5) PSEUDOspherical(real 1.5) Fbeta(real 1.5)]
	if _rc { 
		capture syntax[, MAT(string) POWER(real 1.5) PSEUDOspherical(real 1.5) Fbeta(real 1.5)]
		if _rc { 
			noisily display as err "Wrong input given. Please consult the help file for the correct way of input."

		}
		else{
			mata: excel = xl()
			capture mata: excel.create_book("Classify Metrics", "Probability scores")
			mata: excel.load_book("Classify Metrics")
			capture mata: excel.add_sheet("General metrics")
			capture mata: excel.add_sheet("Class specific metrics")	
			mata: excel.close_book()
		
			mata: classify("Contingency", "Contingency", "`power'", "`pseudospherical'", "`fbeta'", "`mat'")
		}
	}
	else{
		
		mata: excel = xl()
		capture mata: excel.create_book("Classify Metrics", "Probability scores")
		mata: excel.load_book("Classify Metrics")
		capture mata: excel.add_sheet("General metrics")
		capture mata: excel.add_sheet("Class specific metrics")	
		mata: excel.close_book()

		if "`predict'" == "predict" {	
			mata: y = st_data(., "`varlist'")
			mata: yunique = uniqrows(y)
			mata: ycount = rows(yunique)
	
			mata: st_local("ycount", strofreal(ycount))
			mata: stata(`"predict p1-p`ycount'"')
			capture mata: stata(`"ds, has(varlabel "*Pr*") skip(1)"')
			mata: probs = st_global("r(varlist)")
			mata: classify("`varlist'", probs, "`power'", "`pseudospherical'", "`fbeta'")
		} 
		
		if "`predict'" != "predict" {
			mata: classify("`varlist'","`probs'", "`power'", "`pseudospherical'", "`fbeta'")
		}
	}


end