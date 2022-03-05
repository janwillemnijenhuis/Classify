capture program drop classify
program classify, rclass
	version 14
	syntax varlist(min=1) [, var(string asis)]
	mata: classify("`varlist'","`var'")
end
	
