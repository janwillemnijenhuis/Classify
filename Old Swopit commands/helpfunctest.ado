version 14 

mata:

function matrix_mse(residual) {
	result = mean(rowsum(residual :* residual))
	return(result)
}
function running_rowsum(input_matrix) {
	result = input_matrix
	for(i=2; i<=cols(input_matrix); i++) {
		result[,i] = result[,i-1] + result[,i]
	}
	return(result)
}

void print_vector(text, elements) {
	printf(text)
	for(i=1;i<=cols(elements);i++) {
		printf("{bf:%9.4f}", elements[i])
	}
	printf("\n")
}	

void print_matrix(contents, rownames, colnames, | uline, lline, mline, digits, rowtitle, coltitle) {
	// because Stata cannot display matrices with dots in colnames, we need our own printing function!
	n = rows(contents)
	m = cols(contents)
	if (rownames == . | rows(rownames) == 0) {
		rowname_width = 0
		rowname_flag = 0
	} else {
		rowname_width = max(strlen(rownames) \ 10)
		rowname_flag = 1
	}
	if (uline == . | rows(uline) == 0) {
		uline = 0
	} 
	if (lline == . | rows(lline) == 0) {
		lline = 0
	}
	if (mline == . | rows(mline) == 0) {
		mline = (n > 1)
	}
	if (digits == . | rows(digits) == 0) {
		digits = 4
	}
	_colnames = colnames
	if (cols(_colnames) > 1){
		_colnames = _colnames'
	}
	
	if (rowtitle == . | rows(rowtitle) == 0) {
		rowtitle_rows = 0
	} else {
		// todo: ensure that rowname_flag is true
		rowtitle_rows = rows(rowtitle)
		rowname_width = max((strlen(rowtitle) \ rowname_width))
	}
	if (coltitle == . | rows(coltitle) == 0) {
		coltitle_rows = 0
	} else {
		// todo: ensure that rowname_flag is true
		coltitle_rows = rows(coltitle)
	}
	
	colwidths = rowmax((strlen(_colnames) :+ 3 , J(rows(_colnames), 1, 6)))
	// todo: support word wrap for long colnames and maybe row and col titles
	// todo: make colwidths depend on the contents
	// todo: support lines before totals
	numberf = strofreal(digits) + "f"
	if (rowname_flag) {
		hline = "{hline " + strofreal(rowname_width+1)+ "}{c +}{hline " + strofreal(sum(colwidths :+ 1) + 2)+ "}\n"
	} else {
		hline = "{hline " + strofreal(rowname_width+1+1+sum(colwidths :+ 1) + 2) + "}\n"
	}
	// print header
	if (uline) {
		printf(hline)
	}
	
	if (rowtitle_rows > 1) {
		for(i=1; i <= rowtitle_rows; i++) {
			// todo: take into accoutn possible difference in vlines
			printf("%" + strofreal(rowname_width) + "s {c |}", rowtitle[i])
			// todo: make coltitle centered
			if (coltitle_rows > 0) {
				coltitle_current =  i + coltitle_rows - rowtitle_rows + 1
				if ((coltitle_current > 0) & (coltitle_current <= coltitle_rows)) {
					printf(coltitle[coltitle_current])
				}
			}
			if (i < rowtitle_rows) {
				printf("\n")
			}
		}
	} else if (rowtitle_rows==1){
	    printf("%" + strofreal(rowname_width) + "s {c |} ", rowtitle)	    
	} else if (rowname_flag) {
		printf("%" + strofreal(rowname_width) + "s {c |} ", "")
	}
	for(j=1; j<=m; j++){
		displayas("txt")
		printf("%" + strofreal(colwidths[j]) + "s ", colnames[j])
	}
	printf("\n")
	if (mline) {
		printf(hline)
	}
	// print the rest of the table
	if (coltitle_rows==1){
		displayas("txt")
	    printf("%"+strofreal(rowname_width)+ "s {c |}\n",coltitle)
	} else if (coltitle_rows>1){
	    "A higher (>1) number of words for column title is not yet supported"
	}
	for(i=1; i<=n; i++) {
		if (rowname_flag) {
			displayas("txt")
			printf("%" + strofreal(rowname_width)+ "s {c |} ", rownames[i])
		}
		for(j=1; j<=m; j++){
			displayas("res")
			printf("%" + strofreal(colwidths[j]) + "." + numberf + " ", contents[i, j])
		}
		printf("\n")
	}
	if (lline) {
		printf(hline)
	}
}

end
