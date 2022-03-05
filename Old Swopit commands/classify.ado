capture program drop classify
program classify, rclass
	version 14
	syntax [, var(string asis)]
	mata: classify("`var'")
end
	
