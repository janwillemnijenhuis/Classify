version 14
mata
function classify(string varlist, string AtClasses) {
	list = tokens(AtClasses)
	y = st_data(., varlist)
	X = st_data(., list)
	allcat = uniqrows(y)
	rowMax = rowmax(X)
	nrow = rows(X)
	ncol = cols(X)
	
	q = J(nrow, ncol, 0)
	for(i=1; i<=ncol; i++) {
			q[.,i] = (y :== allcat[i])
	}

	conf_mat = J(ncol, ncol, 0)
	prediction = J(nrow, 1, 0)
	for (i=1;i<=nrow;i++) {
		for (j=1;j<=ncol;j++) {
			if (X[i,j] == rowMax[i]) {
				prediction[i] = allcat[j]
			}
		}
	}
	
	for (i=1;i<=nrow;i++) {
		j = prediction[i]
		k = y[i]
		for (l=1;l<=ncol;l++) {
			cat = allcat[l]
			if (cat == j) {
				s = l
			} 
			if (cat == k) {
				t = l
			}
		}
		conf_mat[s,t]=conf_mat[s,t]+1
	}

	colstripes = "y=" :+ strofreal(allcat) 
	rowstripes = "y=" :+ strofreal(allcat) 
	rowtitle = "True"
	coltitle = "Predicted"
	
	accuracy = sum((y :== prediction))/nrow
	brier_vec = J(nrow, 1, 0)
	for (i=1;i<=nrow;i++) {
		for (j=1;j<=ncol;j++) {
			if (y[i] == j) {
				brier_vec[i] = brier_vec[i] + (1 - X[i,j])^2
			} else {
				brier_vec[i] = brier_vec[i] + (0 - X[i,j])^2
			}
		}
	}
	brier_score = sum(brier_vec)/nrow

	ranked_probability_score = matrix_mse(running_rowsum(X) - running_rowsum(q))

	tot = sum(conf_mat)
 	all_pos = colsum(conf_mat)'
 	all_neg = tot :- all_pos
 	pre_pos = rowsum(conf_mat)
 	pre_neg = tot :- pre_pos
	
 	true_pos = diagonal(conf_mat)
 	fals_pos = pre_pos - true_pos
 	fals_neg = all_pos - true_pos
 	true_neg = pre_neg - fals_neg
	
 	noise = fals_pos :/ all_neg
 	recall = true_pos :/ all_pos
 	precision = true_pos :/ pre_pos
 	n2s = noise :/ recall
	
 	result = precision, recall, n2s
 	colname = "Precision" \  "Recall" \  "Adj. noise-to-signal"
 	rowname = "y=" :+ strofreal(allcat) 
 	print_matrix(result, rowname, colname,., ., ., 4, ., .)
 	displayas("txt")
 	printf("\n")
 	printf("Accuracy                 = {bf:%9.4f} \n", accuracy)
 	printf("Brier score              = {bf:%9.4f} \n", brier_score)
 	printf("Ranked probability score = {bf:%9.4f} \n", ranked_probability_score)
 	printf("\nConfusion Matrix\n")
 	print_matrix(conf_mat, rowstripes, colstripes,., ., ., 0, rowtitle, coltitle)
}

end