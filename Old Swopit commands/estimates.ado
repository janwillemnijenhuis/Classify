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
	
	// computation of multiclass measures
	accuracy = sum((y :== prediction))/nrow
	misclassification_rate = 1 - accuracy
	balanced_acc_row = J(1, ncol, 0)
	colsum_conf_mat = colsum(conf_mat) // this is a_1...a_i
	rowsum_conf_mat = rowsum(conf_mat) // this is p_1...p_i
	a2 = 0
	p2 = 0
	pa = 0
	p_plus_a_2 = 0
	for (i=1;i<=ncol;i++) {
	    p_plus_a_2 = p_plus_a_2 + (colsum_conf_mat[i]*rowsum_conf_mat[i])^2
	    a2 = a2 + colsum_conf_mat[i]^2
		p2 = p2 + rowsum_conf_mat[i]^2
		pa = pa + colsum_conf_mat[i]*rowsum_conf_mat[i]
	    balanced_acc_row[i] = conf_mat[i,i]/colsum_conf_mat[i]
	}
	balanced_acc = sum(balanced_acc_row)/ncol
	corr_matthew = (sum((y :== prediction))*nrow - pa)/(sqrt((nrow^2-p2)*(nrow^2-a2)))
	cohen_kappa = (sum((y :== prediction))*nrow - pa)/(nrow^2-pa)
	scott_pi = (sum((y :== prediction))*nrow - 0.25*p_plus_a_2)/(nrow^2-0.25*p_plus_a_2)
	peirce_skill_score = (sum((y :== prediction))*nrow - pa)/(nrow^2-a2)
	clayton_skill_score = (sum((y :== prediction))*nrow - pa)/(nrow^2-p2)
	
	// computation of class specific measures
	true_pos_k = J(1,ncol,0)
	true_neg_k = J(1,ncol,0)
	false_pos_k = J(1,ncol,0)
	false_neg_k = J(1,ncol,0)
	for (i=1;i<=ncol;i++) {
	    true_pos_k[i] = conf_mat[i,i]
		false_pos_k[i] = rowsum_conf_mat[i] - conf_mat[i,i]
		false_neg_k[i] = rowsum_conf_mat[i] - conf_mat[i,i]
		true_neg_k[i] = nrow - true_pos_k[i] - false_pos_k[i] - false_neg_k[i]
	}
	precision_k = true_pos_k:/(true_pos_k :+ false_pos_k)
	nega_pred_val_k = true_neg_k:/(true_neg_k:+false_neg_k)
	recall_k = true_pos_k:/(true_pos_k:+false_neg_k)
	specificity_k = true_neg_k:/(true_neg_k:+false_pos_k)
	balanced_acc_k = true_pos_k:/(2*(true_pos_k:+false_neg_k)) + true_neg_k:/(2*(true_neg_k:+false_pos_k))
	prevalence_k = (true_pos_k:+false_neg_k):/nrow
	false_pos_rate_k = false_pos_k:/(true_neg_k:+false_pos_k)
	false_alarm_rate_k = false_pos_k:/(true_pos_k:+false_pos_k)
	false_neg_rate_k = false_neg_k:/(true_pos_k:+false_neg_k)
	false_omm_rate_k = false_neg_k:/(true_neg_k:+false_neg_k)
	bias_score_k = (true_pos_k:+false_pos_k):/(true_pos_k:+false_neg_k)
	pos_lik_ratio_k = (true_pos_k:*(false_pos_k:+true_neg_k)):/(false_pos_k:*(true_pos_k:+false_neg_k))
	neg_lik_ratio_k = (false_neg_k:*(false_pos_k:+true_neg_k)):/(true_neg_k:*(true_pos_k:+false_neg_k))
	youden_j_k = recall_k :+ specificity_k :- 1
	markedness_k = precision_k :+ nega_pred_val_k :- 1
	informedness_k = recall_k :+ specificity_k :- 1
	diagnostic_odds_k = (true_pos_k:*true_neg_k):/(false_pos_k:*false_neg_k)
	yule_q_coeff_k = (true_pos_k:*true_neg_k :- false_neg_k:*false_pos_k):/(true_pos_k:*true_neg_k :+ false_neg_k:*false_pos_k)
	yule_y_coeff_k = (sqrt(true_pos_k:*true_neg_k) :- sqrt(false_neg_k:*false_pos_k)):/(sqrt(true_pos_k:*true_neg_k) :+ sqrt(false_neg_k:*false_pos_k))
	fowlkes_mallows_k = true_pos_k :/ sqrt(((true_pos_k :+ false_neg_k):*(true_pos_k :+ false_pos_k)))
	f1_score_k = 2*true_pos_k:/(2*true_pos_k :+ false_pos_k :+ false_neg_k)
	beta = 0.5 // we should actually ask the user for this value
	fb_score_k = ((1+beta^2):*precision_k:*recall_k):/(beta^2:*precision_k :+ recall_k)
	matthew_corr_k = (true_pos_k:*true_neg_k :- false_pos_k:*false_neg_k):/(sqrt((true_pos_k:+false_pos_k):*(true_pos_k:+false_neg_k):*(true_neg_k:+false_pos_k):*(true_neg_k:+false_neg_k)))
	threat_k = true_pos_k :/ (true_pos_k :+ false_neg_k :+ false_pos_k)
	a = ((true_pos_k:+false_pos_k):*(true_pos_k:+false_neg_k)):/nrow
	gilbert_skill_score_k = (true_pos_k:-a):/(true_pos_k:+false_neg_k:+false_pos_k:-a)
	scott_pi_k = ((nrow:*true_pos_k):-0.25:*(2:*true_pos_k:+false_pos_k:+false_neg_k):^2):/(nrow^2 :- (0.25:*(2:*true_pos_k:+false_pos_k:+false_neg_k):^2))
	peirce_skill_score_k = ((true_pos_k:*true_neg_k):-(false_neg_k:*false_pos_k)):/((true_pos_k:+false_neg_k):*(false_pos_k:+true_neg_k))
	cohen_kappa_k = ((true_pos_k :- (true_pos_k :+ false_pos_k)) :* ((true_pos_k:+false_neg_k):/nrow) :/ (nrow :- (true_pos_k :+ false_pos_k) :* (true_pos_k :+ false_neg_k)/nrow))
	clayton_skill_score_k = ((true_pos_k :* true_neg_k) :- (false_neg_k :* false_pos_k)) :/ ((true_pos_k :+false_pos_k) :* (false_neg_k :+ true_neg_k))
	extr_dep_score_k = (2 :* log((true_pos_k :+ false_neg_k) :/ nrow)) :/ log(true_pos_k :/ nrow) :- 1
	symmetric_extr_dep_score_k = (log((true_pos_k :+ false_neg_k) :/ nrow) :+ log((true_pos_k :+ true_neg_k) :/ nrow)) :/ log(true_pos_k :/ nrow) :- 1
	prev_threshold_k = sqrt(false_pos_rate_k):/(sqrt(false_pos_rate_k):+sqrt(recall_k))
	adj_noise_to_signal_k = false_pos_rate_k:/recall_k
	
	// macro measures
	prec_macro = sum(precision_k)/ncol
	recall_macro = sum(recall_k)/ncol
	f1_macro = 2*(prec_macro*recall_macro)/(prec_macro+recall_macro)
	fb_macro = (1+beta^2)*(prec_macro*recall_macro)/(beta^2*(prec_macro+recall_macro))
	fowlkes_mallows_index = sqrt(prec_macro+recall_macro)
	
	// weighted averages (do we want to do this?)
	// also do we want each confusion matrix class specific?
		
	// brier score
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
	displayas("txt")
	printf("\nConfusion Matrix\n")
 	print_matrix(conf_mat, rowstripes, colstripes,., ., ., 0, rowtitle, coltitle)
	printf("\n")
 	colname = "Precision" \  "Recall" \  "Adj. noise-to-signal"
 	rowname = "y=" :+ strofreal(allcat) 
 	print_matrix(result, rowname, colname,., ., ., 4, ., .)
 	displayas("txt")
 	printf("\nMulticlass metrics\n")
 	printf("Accuracy                 = {bf:%9.4f} \n", accuracy)
	printf("Misclassification rate	 = {bf:%9.4f} \n", misclassification_rate)
	printf("Balanced accuracy   	 = {bf:%9.4f} \n", balanced_acc)
	printf("Correlation (Matthew)    = {bf:%9.4f} \n", corr_matthew)
	printf("Cohen's kappa coefficient= {bf:%9.4f} \n", cohen_kappa)
	printf("Scott's pi coefficient   = {bf:%9.4f} \n", scott_pi)
	printf("Peirce's skill score     = {bf:%9.4f} \n", peirce_skill_score)
	printf("Clayton skill score      = {bf:%9.4f} \n", clayton_skill_score)
	       
	printf("\nClass specific metrics\n")
	print_vector("Precision                = ", precision_k)
	print_vector("Negative predicted value = ", nega_pred_val_k)
	print_vector("Recall                   = ", recall_k)
	print_vector("Specificity              = ", specificity_k)
	print_vector("Balanced accuracy        = ", balanced_acc_k)
	print_vector("Prevalence               = ", prevalence_k)
	print_vector("False positive rate      = ", false_pos_k)
	print_vector("False alarm rate         = ", false_alarm_rate_k)
	print_vector("False negative rate      = ", false_neg_rate_k)
	print_vector("False omission rate      = ", false_omm_rate_k)
	print_vector("Bias score               = ", bias_score_k)
	print_vector("Positive likelihood ratio= ", pos_lik_ratio_k)
	print_vector("Negative likelihood ratio= ", neg_lik_ratio_k)
	print_vector("Youden's J statistic     = ", youden_j_k)  
	print_vector("Markedness               = ", markedness_k)
	print_vector("Informedness             = ", informedness_k)
	print_vector("Diagnostic odds ratio    = ", diagnostic_odds_k)
	print_vector("Yule's Q coefficient     = ", yule_q_coeff_k)
	print_vector("Yule's Y coefficient     = ", yule_y_coeff_k)
	print_vector("Fowlkes-Mallows index    = ", fowlkes_mallows_k)
	print_vector("F1-score                 = ", f1_score_k)
	print_vector("FB-score	               = ", fb_score_k)
	print_vector("Correlation (Matthews)   = ", matthew_corr_k)
	print_vector("Threat score             = ", threat_k)
	print_vector("Gilbert Skill Score      = ", gilbert_skill_score_k)
	print_vector("Scott's pi coefficient   = ", scott_pi_k)
	print_vector("Peirce's skill score     = ", peirce_skill_score_k)
	print_vector("Cohen's kappa coefficient= ", cohen_kappa_k)
	print_vector("Clayton Skill Score      = ", clayton_skill_score_k)
	print_vector("Extreme dependency score = ", extr_dep_score_k)
	print_vector("Symm. extr. dep. score   = ", symmetric_extr_dep_score_k)
	print_vector("Prevalence threshold     = ", prev_threshold_k)
	print_vector("Adj. N2S ratio           = ", adj_noise_to_signal_k)
	
	printf("\nMacro averages of class-specific metrics\n")
	printf("F1-score                 = {bf:%9.4f}\n", f1_macro)
	printf("FB-score                 = {bf:%9.4f}\n", fb_macro)
	printf("Fowlkes-Mallows index    = {bf:%9.4f}\n", fowlkes_mallows_index)
	
	
	printf("\nOther\n")
 	printf("Brier score              = {bf:%9.4f} \n", brier_score)
 	printf("Ranked probability score = {bf:%9.4f} \n", ranked_probability_score)
	
	
}

end