capture program drop classify
program classify, rclass
	version 14
	syntax varlist(min=1) [, var(string asis) PREDict]

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
		mata: var = st_global("r(varlist)")
		mata: classify("`varlist'", var)
	} 
	
	if "`predict'" != "predict" {
		mata: classify("`varlist'","`var'")
	}
	

end