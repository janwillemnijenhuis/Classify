version 14
mata
function classify(string varlist, string AtClasses) {
	list = tokens(AtClasses)
	y = st_data(., varlist)
	X = st_data(., list)
	allcat = uniqrows(y)
	ncat = rows(allcat)
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
	rogot_goldberg = 0
	gktm = 0 // goodman kruskal tau multi
	gkg_1= 0
	gkg_2=0
	theil_u_1=0
	theil_u_2=0
	stuart_1=0
	stuart_2=0
	D_j = J(1,ncol,0)
	colsum_conf_mat = colsum(conf_mat) // this is a_1...a_i
	rowsum_conf_mat = rowsum(conf_mat) // this is p_1...p_i
	a2 = 0
	p2 = 0
	pa = 0
	p_plus_a_2 = 0
	n_ij_abs = 0
	n_ij_k =0
	n_ij_k_abs = 0
	pcm = 0 // pearsion chi multi
	lik_g = 0 // likelihood g
	adj_r_1=0
	adj_r_2=0
	adj_r_3=0
	kendall_1=0
	kendall_2=0
	for (i=1;i<=ncol;i++) {
	    p_plus_a_2 = p_plus_a_2 + (colsum_conf_mat[i]*rowsum_conf_mat[i])^2
	    a2 = a2 + colsum_conf_mat[i]^2
		p2 = p2 + rowsum_conf_mat[i]^2
		pa = pa + colsum_conf_mat[i]*rowsum_conf_mat[i]
	    balanced_acc_row[i] = conf_mat[i,i]/colsum_conf_mat[i]
		rogot_goldberg = rogot_goldberg + balanced_acc_row[i] + conf_mat[i,i]/rowsum_conf_mat[i] 
		gktm_2 = 0 
		for (j=1;j<=ncol;j++) {
		    kendall_1 = kendall_1 + rowsum_conf_mat[i]*(rowsum_conf_mat[i]-1)/2
			kendall_2 = kendall_2 + colsum_conf_mat[i]*(colsum_conf_mat[i]-1)/2
		    // this part is for gerity skill score
			n_r = 0
			for (r=1;r<=j;r++) {
			    n_r = n_r+colsum_conf_mat[r]
			}
			D_j[j] = (nrow-n_r)/n_r
			// this is for tonnies measure
			if (abs(i-j)==1) {
			    n_ij_abs = n_ij_abs+conf_mat[i,j]
			}
			if (i+j-1==ncol) {
			    n_ij_k = n_ij_k + conf_mat[i,j]
			}
			if (i+j-1==ncol-1 || i+j-1==ncol+1) {
			    n_ij_k_abs = n_ij_k_abs + conf_mat[i,j]
			}
			//
		    n_ij = conf_mat[i,j]
			// this is all for the adj_rand_index
			if (n_ij>1){
				adj_r_1 = adj_r_1 + comb(n_ij,2)
			}
			if (rowsum_conf_mat[i]>1) {
				adj_r_2 = adj_r_2 + comb(rowsum_conf_mat[i],2)
			}
			if (colsum_conf_mat[j]>1) {
				adj_r_3 = adj_r_3 + comb(colsum_conf_mat[j],2)
			}
		
			gktm_2 = gktm_2 + conf_mat[i,j]^2/rowsum_conf_mat[i]
			n_hmax_kmax = 0
			n_hmax_kmin = 0
			for (h=i+1;h<=ncol;h++){
			    for (k=j+1;k<=ncol;k++){
				    n_hmax_kmax = n_hmax_kmax + conf_mat[h,k]
				}
				for (k=1;k<j;k++) {
				    n_hmax_kmin = n_hmax_kmin + conf_mat[h,k]
				}
			}
			gkg_1 = gkg_1 + n_ij*n_hmax_kmax
			gkg_2 = gkg_2 + n_ij*n_hmax_kmin
			stuart_1 = stuart_1 + n_ij*n_hmax_kmax
			stuart_2 = stuart_2 + n_ij*n_hmax_kmin
			
			if ((n_ij/(rowsum_conf_mat[i]+colsum_conf_mat[j]))>0) {
				theil_u_1 = theil_u_1 + n_ij*log((n_ij/(rowsum_conf_mat[i]+colsum_conf_mat[j])))
				lik_g = lik_g + log((n_ij*nrow)/(rowsum_conf_mat[i]*colsum_conf_mat[j]))
			} // is this right?
			theil_u_2 = theil_u_2 + log(colsum_conf_mat[j])*colsum_conf_mat[j]
			
			pcm = pcm + (n_ij-rowsum_conf_mat[i]*colsum_conf_mat[j]/nrow)^2/(rowsum_conf_mat[i]*colsum_conf_mat[j]/nrow)
			
		}
		gktm = gktm + (nrow*gktm_2-colsum_conf_mat[i]^2)/(nrow^2-colsum_conf_mat[i]^2)
	}
	theil_u = -1*theil_u_1/theil_u_2
	
	gsk = 0 // gerity skill score
	for (i=1;i<=ncol;i++) {
	    for (j=1;j<=ncol;j++) {
		    w=0
		    if (i==j) {
			    w_jj = 0
				w_jj_1 = 0
				w_jj_2 = 0
				for (r=1;r<=j-1;r++) {
				    w_jj_1 = w_jj_1 + 1/D_j[r]
				}
				for (r=j;r<=ncol;r++) {
				    w_jj_2 = w_jj_2 + D_j[r]
				}
				w_jj = (w_jj_1+w_jj_2)/(ncol-1)
				w = w_jj
			} else {
			    w_ij = 0
				w_ij_1 = 0
				w_ij_2 = 0
				for (r=1;r<=j-1;r++) {
				    w_ij_1 = w_ij_1 + 1/D_j[r]
				}
				for (r=j;r<=ncol;r++) {
				    w_ij_2 = w_ij_2 + D_j[r]
				}
				w_ij = (w_ij_1+w_ij_2-(j-1))/(ncol-1)
				w = w_ij
			}
			gsk = gsk+w*conf_mat[i,j]
		}
	}
	gsk = gsk/nrow^2
	
	balanced_acc = sum(balanced_acc_row)/ncol
	corr_matthew = (sum((y :== prediction))*nrow - pa)/(sqrt((nrow^2-p2)*(nrow^2-a2)))
	cohen_kappa = (sum((y :== prediction))*nrow - pa)/(nrow^2-pa)
	scott_pi = (sum((y :== prediction))*nrow - 0.25*p_plus_a_2)/(nrow^2-0.25*p_plus_a_2)
	peirce_skill_score = (sum((y :== prediction))*nrow - pa)/(nrow^2-a2)
	clayton_skill_score = (sum((y :== prediction))*nrow - pa)/(nrow^2-p2)
	// 2nd draft, from page 29 onwards
	n_kk = sum((y:==prediction))
	rogot_goldberg = rogot_goldberg/(2*ncol)
	goodman_kruskal_lambda = (n_kk-rowmax(colsum_conf_mat))/(nrow-rowmax(colsum_conf_mat))
	// gktm above
	gkg = (gkg_1-gkg_2)/(gkg_1+gkg_2)
	stuart = (stuart_1-stuart_2)/((ncol-1)/ncol)
	tonnies = 2*n_kk + n_ij_abs - 2* n_ij_k-n_ij_k_abs
	koppen_ms = accuracy + (1/(2*nrow))*n_ij_abs
	msc = pcm / nrow
	lik_g = 2*lik_g
	c_contingency = sqrt((pcm/nrow)/(1+(pcm/nrow)))
	cramer_v_multi = sqrt((pcm/nrow)) // i took the easier approach
	cramer_v_bias_multi = sqrt(pcm/(nrow*(2-(1/(nrow-1)))))
	adj_rand_index = (adj_r_1-adj_r_2*adj_r_3/comb(nrow,2))/(0.5*(adj_r_2+adj_r_3)-adj_r_2*adj_r_3/comb(nrow,2))
	kendall_tau = (stuart_1-stuart_2)/sqrt((comb(nrow,2)-kendall_1)*(comb(nrow,2)-kendall_2))
	
	// computation of class specific measures
	true_pos_k = J(1,ncol,0)
	true_neg_k = J(1,ncol,0)
	false_pos_k = J(1,ncol,0)
	false_neg_k = J(1,ncol,0)
	for (i=1;i<=ncol;i++) {
	    true_pos_k[i] = conf_mat[i,i]
		false_pos_k[i] = rowsum_conf_mat[i] - conf_mat[i,i]
		false_neg_k[i] = colsum_conf_mat[i] - conf_mat[i,i]
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
	cohen_kappa_k = (true_pos_k :- (true_pos_k :+ false_pos_k) :* (true_pos_k:+false_neg_k):/nrow) :/ (nrow :- (true_pos_k :+ false_pos_k) :* (true_pos_k :+ false_neg_k)/nrow)
	clayton_skill_score_k = ((true_pos_k :* true_neg_k) :- (false_neg_k :* false_pos_k)) :/ ((true_pos_k :+false_pos_k) :* (false_neg_k :+ true_neg_k))
	extr_dep_score_k = (2 :* log((true_pos_k :+ false_neg_k) :/ nrow)) :/ log(true_pos_k :/ nrow) :- 1
	symmetric_extr_dep_score_k = (log((true_pos_k :+ false_neg_k) :/ nrow) :+ log((true_pos_k :+ false_pos_k) :/ nrow)) :/ log(true_pos_k :/ nrow) :- 1
	prev_threshold_k = sqrt(false_pos_rate_k):/(sqrt(false_pos_rate_k):+sqrt(recall_k))
	adj_noise_to_signal_k = false_pos_rate_k:/recall_k


	//From here after next draft
	consonni_todeschini = log(true_pos_k :+ 1) :/ log(true_pos_k :+ false_pos_k :+ false_neg_k :+ 1)
	van_der_maarel = (2 :* true_pos_k :- false_pos_k :- false_neg_k) :/ (2 :* true_pos_k :+ false_pos_k :+ false_neg_k)
	anderberg_one =  8 :* true_pos_k :/ (8 :* true_pos_k :+ false_pos_k :+ false_neg_k)
	anderberg_two = true_pos_k :/ (true_pos_k :+ 2:* (false_pos_k :+ false_neg_k))
	benini = (true_pos_k :- (true_pos_k :+ false_neg_k) :* (true_pos_k :+ false_pos_k)):/(true_pos_k :+ colmin((false_neg_k\ false_pos_k)) :- (true_pos_k :+ false_neg_k) :* (true_pos_k :+ false_pos_k))
	russel_rao = true_pos_k:/nrow
	kulczynski_one = true_pos_k :/ (false_pos_k :+ false_neg_k)
	kulczynski_two = 0.5 :* (true_pos_k :/ (true_pos_k :+ false_pos_k) + true_pos_k :/ (true_pos_k :+ false_neg_k))
	johnson = true_pos_k :/ (true_pos_k :+ false_pos_k) + true_pos_k :/ (true_pos_k :+ false_neg_k)
	mcconnaughey = (true_pos_k:^2 :- false_neg_k :* false_pos_k):/( (true_pos_k :+ false_pos_k):* (true_pos_k :+ false_neg_k) )
	driver_kroeber = true_pos_k :/ sqrt((true_pos_k :+ false_pos_k) :* (true_pos_k :+ false_neg_k))
	sorensen = 4 :* true_pos_k :/ (4 :* true_pos_k :+ false_pos_k :+ false_neg_k)
	jaccard = 3 :* true_pos_k :/ (3 :* true_pos_k :+ false_pos_k :+ false_neg_k)

	alpha_rb = 1
	rijsbergen = (false_neg_k :+ alpha_rb :* (false_pos_k :- false_neg_k)):/(true_pos_k :+ false_neg_k :+ alpha_rb :* (false_pos_k :- false_neg_k))
	
	upholt_s = (0.5:*(-f1_score_k :+ sqrt(f1_score_k:^2 :+ 8:* f1_score_k))):^(1:/nrow)


	//From here im just gonna redefine true_pos etc, we can change earlier later

	fn = false_neg_k 
	fp = false_pos_k
	tn = true_neg_k
	tp = true_pos_k

	gini_index = (tp :- (tp :+fp):*(tp :+fn)):/sqrt(1:-(tp:+fp):^2 :* (1 :- (tp :+ fn):^2))
	modified_gini = (tp :- (tp :+fp):*(tp :+fn)):/ (1 :- abs(fp :- fn):/2 :- (tp :+fp):*(tp :+fn))
	norm_coallacation = tp :/ (fn :+ fp :- tp)
	savage = 1 :- tp :/ (tp :+ colmax((fp\fn)))
	sokal_sneath = tp :/ (tp :+ 2:* (fp :+ fn))
	hellinger = 2:* sqrt(1:- tp:/sqrt( (tp :+ fp):*(tp :+ fn) ))
	lance_williams = (fn :+ fp):/ ( 2:* tp :+ fn :+ fp)
	mountford = 2:* tp :/ (2:*fp:*fn :+ tp :*fp :+ tp:*fn  )
	michelet = tp:^2 :/ (fp :* fn)
	sorgenfrei = tp:^2 :/ ((tp:+fp):*(tp:+fn))
	fager = tp:/sqrt(((tp:+fp):*(tp:+fn))) :- colmax((fp\fn))
	fager_mcgowan = tp:/sqrt(((tp:+fp):*(tp:+fn))) :- 1:/(2:*sqrt(colmax(((tp :+ fp)\(tp :+ fn)))))
	u_cost = log(1 :+ (colmin((fp\fn)) :+ tp):/(colmax((fp\fn)) :+ tp))
	s_cost = log(1:+ (colmin((fp\fn)) :+ tp):/ (1 :+ tp) ):^ (-0.5)
	r_cost = log(1:+ tp:/(tp:+fp)):*log(1:+tp:/(tp:+fn))
	t_cost = sqrt(u_cost :* r_cost :* s_cost)
	simpson = tp:/(tp:+ colmin((fp\fn)))
	kent_foster_one = (-fp:*fn):/(fp:*(tp:+fp):+fn:*(tp:+fn):+fp:*fn)
	braum_blanquet = tp:/(tp:+ colmax((fp\fn)))
	chord_dissimilarity = sqrt(2:*(1:-tp:/sqrt((tp:+fp):*(tp:+fn))))
	rousseau_skill = (4:*tp :- (2:*tp :+fn:+fp):^2):/( 2:* (2:*tp :+ fn :+ fp) :- (2:* tp :+ fn :+ fp):^2)
	anonymous_one = (2:*tn):/(2:*tn  :+ fp :+ fn)
	anonymous_two = (2:*tn  :- fp :- fn):/(2:*tn  :+ fp :+ fn)
	soergel_distance = (fn :+ fp) :/ (fn :+ fp:+tn)
	kent_foster_two = (-fp:*fn):/(fp:*(tn:+fn):+fn:*(tn:+fp):+fp:*fn)
	consonni_todeschini_one = log(1:+ tp :+ tn):/log(1:+nrow)
	consonni_todeschini_two = (log(1:+nrow) :- log(1:+fp:+fn)):/log(1:+nrow)
	consonni_todeschini_three = log(1:+tp):/log(1:+nrow)
	consonni_todeschini_five = (log(1:+tp:*tn):-log(1:+fp:*fn)):/log(1:+nrow^2/4)
	matching = (tp :+ tn):/ (2:*tp :+ fn :+ fp)
	balance_error_rate = 1:- tp:/(2:*(tp:+fn)) - tn:/(2:*(tn:+fp))
	austin_colwell = 2/pi():*asin(sqrt((tp:+tn):/nrow))
	dominance = tp:/(tp:+fn)  :- tn:/(tn:+fp)
	geometric_mean = sqrt(tp:/(tp:+fn)  :* tn:/(tn:+fp))
	adjusted_g_mean = colmin(J(1,cols(tp),1)\tp) :* (sqrt(tp:/(tp:+fn)  :* tn:/(tn:+fp)) :+ tn):/(1:+fp:+tn)
	optimization = (tp :+ tn):/nrow :- abs(tp:/(tp:+fn):-tp:/(tn:+fp)):/(tp:/(tp:+fn):+tp:/(tn:+fp))
	positive_lik = tp:*(fp:+tn):/(fp:*(tp:+fn))
	negative_lik = fn:*(fp:+tn):/(tn:*(tp:+fn))
	

	anderberg_a = colmax((tp\fp)) :+ colmax((fn\tn)) :+ colmax((tp\fn)) :+ colmax((fp\tn))
	anderberg_b = colmax(((tp:+fn)\(fp:+tn))) :+ colmax(((tp:+fp)\(fn:+tn)))
	anderberg_d = (anderberg_a :- anderberg_b):/(2:*nrow)

	ample_sim = abs(tp:/(tp:+fp) :- fn:/(fn:+tn))
	relative_dec = (colmax((tp\fp)) :+ colmax((fn\tn)) :- colmax(((tp:+fn)\(fp:+tn)))):/ (1 :- colmax(((tp:+fn)\(fp:+tn))))
	baulieu_one = (nrow:^2 :- (fp :- fn):^2):/(nrow:^2)
	baulieu_two = (nrow:^2 :- 4:*fp:*fn):/(nrow:^2)
	baulieu_three = (tp:*tn :- fp:*fn):/(nrow:^4)
	baulieu_four = (nrow:^2 :- nrow:*(fp:+fn):+(fp:-fn):^2):/(nrow:^2)
	maron_kuhns = (tp:*tn :- fp:*fn):/nrow
	benini_repulsion = (tp:*tn :- fp:*fn):/( colmin(((tp:+fn)\(tp:+fp))) :- (tp:+fn):*(tp:+fp))
	cole_one = (tp:*tn:-fp:*fn):/( (tn:+fp):*(tp:+fp))
	cole_two = (tp:*tn:-fp:*fn):/( (fn:+fn):*(tp:+fn))
	cole_three = (tp:*tn:-fp:*fn):/ colmin((( (tp:+fn):*(tp:+fp) )\( (tn:+fn):*(tn:+fp) )))
	cole_four = sqrt(2):*(tp:*tn:-fp:*fn):/sqrt((tp:*tn :-fp:*fn):^2 :- (tp:+fp):*(tp:+fn):*(tn:+fp):*(tn:+fn))
	digby = ( (tp:*tn):^(3/4) :- (fp:*fn):^(3/4) ):/( (tp:*tn):^(3/4) :+ (fp:*fn):^(3/4) )
	hawkin_dotson = 1/2:*(tp:/(tp:+fp:+fn) :+ tn:/(tn:+fp:+fn))
	pattern_diff = 4:*fp:*fn:/(nrow:^2)
	size_diff = (fp:+fn):^2:/(nrow:^2)
	shape_diff = (nrow:*(fp:+fn) :- (fp:-fn):^2):/(nrow:^2)
	sneath_patt_diff = 2:*sqrt(fp:*fn):/nrow
	koppen = ( (tp:+fp):*(nrow:-tp:-fp):-fn):/( (tp:+fp):*(nrow:-tp:-fp))
	binary_shape = (nrow:*(fp:+fn):-(fp:-fn):^2):/(nrow:^2)
	rogot_goldberg_one = tp:/(2:*tp:+fp:+fn) :+ tn:/(2:*tn:+fp:+fn)
	lacour = tp:*(fn:+tn):/( (tp:+fp):*fn )
	peirce_one = (tp:*tn :- fp:*fn):/( (tp:+fn):*(fp:*tn) )
	peirce_two = ( tp:*fp :+ fp:* fn ):/( tp:*fp :+ 2:*fp:*fn + fn:*tn)
	
	// from here is page 17 until the end (jw)
	log_odds = log((tp:/fp):*(tn:/fn))
	goodman_ass = log_odds :* sqrt(((tp:+fp):*(tp:+fn):*(tn:+fp):*(tn:+fn)))
	goodman_un_ass = 0.25:*log_odds
	// yule is computed before
	goodman_lambda = (tp:+tn - colmax((tp:+fn\fp:+tn))):/(nrow:-colmax((tp:+fn\fp:+tn)))
	goodman_kruskal_1 = (2:*colmin((tp\tn)):-fp:-fn):/(2:*colmin((tp\tn)):+fp:+fn)
	goodman_kruskal_2 = (colmax((tp\fn)):+colmax((fp\tn)):-colmax(((tp:+fp)\(tn:+fn)))):/(1:-colmax(((tp:+fp)\(tn:+fn))))
	goodman_kruskal_3 = (tp:+tn:-colmax((tp\tn)):-((fp:+fn):/2)):/(1:-colmax((tp\tn)):-((fp:+fn):/2))
	// goodman_kruskal_4 isnt formatted well in the draft, ask about this
	// yule y is computed before
	adj_f_score = sqrt(((5:*tp:^2):/((5:*tp:^2):+4:*tp:*fn+tp:*fp)):*((2:*tn):/(2:*tn:+fn:+fp)))
	pearson_mscc = ((tp:*tn:-fp:*fn):^2):/((tp:+fp):*(tp:+fn):+(tn:+fp):*(tn:+fn))
	// matthew corr k is computed before
	// what is phi in pearson #??? coefficient?
	gower = (tp:+tn):/sqrt((tp:+fp):*(tp:+fn):*(fn:+fp):*(tn:+fn))
	stiles = log((nrow:*(abs(tp:*tn:-fp:*fn):-nrow/2):^2):/((tp:+fp):*(tp:+fn):*(tn:+fp):*(tn:+fn))):/log(10)
	dennis = (tp:*tn:-fp:*fn):/sqrt((tp:+fp):*(tp:+fn):*nrow)
	loevinger = (tp:*tn:-fp:*fn):/(colmin((((tp:+fp):*(fp:+tn))\((tn:+fp):*(tn:+fn)))))
	// ??? this formula is the same as loevinger
	pearson_heron_2 = cos((pi():*sqrt(fp:*fn)):/(sqrt(tp:*tn):+sqrt(fp:*fn)))
	forbes_d = (tp:*nrow):/((tp:+fp):*(tp:+fn))
	alroy_forbes = (tp:*(nrow+sqrt(nrow))):/(tp:*(nrow+sqrt(nrow)):+1.5:*fp:*fn)
	// from here page 20
	forbes_2 = (tp:*nrow:-(tp:+fp):*(tp:+fn)):/(nrow:*colmin((tp:+fp)\(tp:+fn)):-(tp:+fp):*(tp:+fn))
	gini_ass = ((tp:*tn):-(fp:*fn)):/(nrow:-abs(fp:-fn):-(tp:+fp):*(tp:+fn):-(tn:+fp):*(tn:+fn))
	fleiss = ((tp:*tn:-fp:*fn):*(tp:+fp):*(fp:+tn):+(tp:+fn):*(fn:+tn)):/(2:*(tp:+fp):*(tp:+fn):*(tn:+fp):*(tn:+fn))
	kuhns_1 = (2:*(tp:*tn:-fp:*fn)):/(nrow:*(2:*tp:+fp:+fn))
	kuhns_2 = ((tp:*tn:-fp:*fn):*(tp:+fp):*(tp:+fn)):/(nrow:*tp:*(2:*tp:+fp:+fn:-(tp:+fp):*(tp:+fn):/nrow))
	goodman_conc = 2:*(sqrt(tp:*(fn:+tn):*(fp:+tn)):+sqrt(tn:*(tp:+fn):*(tp:+fp)):-sqrt(fp:*(fn:+tn):*(fn:+tp)):-sqrt(fn:*(tp:+fp):*(fp:+tn)))
	rogers_tanimoto = (tp:+tn):/(nrow:+fn:+fp)
	sokal_sneath_2 = 2:*(tp:+tn):/(nrow:+tp:+tn)
	sokal_sneath_3 = (tp:+tn):/(fp:+fn)
	sokal_sneath_4 = 0.25:*((tp:/(tp:+fp)):+(tp:/(tp:+fn)):+(tn:/(fp:+tn)):+(tn:/(fn:+tn)))
	sokal_sneath_5 = (tp:*tn):/(sqrt((tp:+fp):*(tp:+fn):*(tn:+fp):*(tn:+fn)))
	kocher_wang = (tp:*nrow):/((tp:+fp):*(fn:+tn))
	faith = (tp:+tn:/2):/nrow
	hamann = ((tp:+tn):-(fp:+fn)):/(nrow)
	hawkins_dotson = 0.5:*(tp:/(fn:+fp:+tp))+0.5:*(tn:/(fn:+fp:+tn))
	roux_1 = (tp:+tn):/(colmin(fp\fn):+colmin((nrow:-fp)\(nrow:-fn)))
	roux_2 = (nrow:-tp:*tn):/(sqrt((tp:+fp):*(tp:+fn):*(tn:+fp):*(tn:+fn)))
	maxwell_pilliner = 2:*(tp:*tn:-fp:*fn):/((tp:+fp):*(tn:+fn):+(tp:+fn):*(tn:+fp))
	unigram_subtuples = log((tp:*tn):/(fp:*fn)):-3.29:*sqrt((1:/tp):+(1:/fp):+(1:/fn):+(1:/tn))
	norm_google_dist = (colmax(log(fp)\log(fn)):-log(tp)):/(log(nrow):-colmin(log(fp)\log(fn)))
	eyraud_sim_ind = (tp:-(tp:+fp):*(tp:+fn)):/((tp:+fp):*(tp:+fn):*(tn:+fp):*(tn:+fn))
	fossum = (nrow:*(tp:-0.5):^2):/((tp:+fp):*(tp:+fn))
	baroni_urbani_buser_1 = (sqrt(tp:*tn):+tp):/(sqrt(tp:*tn):+tp:+fn:+fp)
	baroni_urbani_buser_2 = (sqrt(tp:*tn):+tp-(fp:+fn)):/(sqrt(tp:*tn):+tp:+fn:+fp)
	michael = 4:*(tp:*tn:-fp:*fn):/((tp:+tn):^2:+(fp:+fn):^2)
	clement = (tp:*(fn:+tn):/(tp:+fp)):+(tn:*(tp:+fp):/(fn:+tn))
	harris_lahey = tp:*(2:*tn:+fn:+fp):/(2:*(tn:+fp:+fn)):+(tn:*(2:*tp:+fn:+fp)):/(2:*(tn:+fp:+fn))
	kuder_richardson = (4:*(tp:*tn:-fp:*fn)):/((tp:+fp):*(fn:+tn):+(tp:+fn):*(fn:+tn):+2:*(tp:*tn:-fp:*fn))
	bc_diss = (4:*fn:*fp):/(nrow:^2)
	sokal_dist = sqrt((fp:+fn):/nrow)
	yule_dist = (2:*fp:*fn):/(tp:*tn:+fp:*fn)
	// gilber skill is already done before
	// so is peirce skill score
	quetelet = (tp:*tn:-fp:*fn):/((tp:+fn):*(fp:+tp))
	koppen = (tp:*tn:-fp:*fn):/((tp:+fp):*(fp:+tn))
	d = ((tp:+fp):*(tp:+fn):+(tn:+fp):*(tn:+fn)):/nrow
	heidke = (tp:+tn:-d):/(nrow:-d)
	shrank = (tp:+tn:-((fp:+fn):/2):-d):/nrow
	tarwid = (nrow:*tp:-(tp:+fp):*(tp:+fn)):/(nrow:*tp:+(tp:+fp):*(tp:+fn))
	weighted_rel_acc = ((4:*fn):/(1:+fn):^2):*((tp:/(tp:+fn)):-(fp:/(fp:+tn)))
	// scott pi has been done
	krippendorf_a = (1:-2:*(nrow:-1):*(fp:+fn)):/(2:*tp:+fp:+fn):*(2:*tn:+fn:+fp)
	// ???
	scott_coeff = ((4:*tp:*tn):-(fn:+fp):^2):/((2:*tp:+fp:+fn):*(2:*tn:+fp:+fn))
	// clayton has been done
	gilber_wells = log((tp:*nrow):/((tp:+fp):*(tp:+fn)))
	weighted_mut_inf = log(((tp:^3):*nrow):/((tp:+fp):*(tp:+fn)))
	//extreme dep done
	// so is symm extr 
	extremal_dep_ind = (log(fp:/(fp:+tn)):-log(tp:/(tp:+fn))):/(log(fp:/(fp:+tn)):+log(tp:/(tp:+fn)))
	symmetric_extremal = (log(fp:/(fp:+tn)):-log(tp:/(tp:+fn)):+log(fn:/(tp:+fn)):-log(tn:/(fp:+tn))):/(log(fp:/(fp:+tn)):+log(tp:/(tp:+fn)):+log(fn:/(tp:+fn)):+log(tn:/(fp:+tn)))
	dunning = 2:*(tp:*log(tp):+fp:*log(fp):+fn:*log(fn):+tn:*log(tn)):-(tp:+fp):*log(tp:+fp):-(tp:+fp):*log(tp:+fn):-(tn:+fp):*log(tn:+fp):-(tn:+fn):*log(tn:+fn):+(tp:+fp:+fn:+tn):*log(tp:+fp:+fn:+tn)
	prevalence_2 = sqrt(fp:/(fp:+tn)):/(sqrt(tp:/(tp:+fn)):+sqrt(fp:/(fp:+tn)))
	tarantula = (tp:/(fp:+tp)):/((tp:/(tp:+fp)):+(fn:/(fn:+tn)))
	discriminant_power = (sqrt(3)/pi()):/(log((tp:*(tn:+fp)):/((tp:+fn):*fp)):+log((tn:*(tp:+fn)):/((tn:+fp):*fn)))
	discr_distance = invnormal(tp:/(tp:+fn)):-invnormal(fp:/(fp:+tn))
	pearson_chi = (nrow:*(tp:*tn:-fp:*fn)):/((tp:+fp):*(tp:+fn):*(tn:+fp):*(tn:+fn)) // didnt check equality with the other thingy
	pearson_chi_yates = (nrow:*(abs(tp:*tn:-fp:*fn):-(nrow:/2)):^2):/((tp:+fp):*(tp:+fn):*(tn:+fp):*(tn:+fn))
	mean_sq_contingency = pearson_chi:/nrow // check the equality there otherwise this doesnt hold
	goodman_kruskal_tau = (nrow:*((tp:^2:/(tp:+fp)):+(fp:^2:/(tp:+fp)):+(fn:^2:/(fn:+tp)):+(tn:^2:/(fn:+tp))):-(tp:+fn):^2:-(fn:+tn):^2):/(2:*(tp:+fn):*(fn:+tn)) // i didnt check it
	likelihood_chi = 2:*tp:*log((tp:*nrow):/((tp:+fp):*(tp:+fn))):+2:*fp:*log((fp:*nrow):/((tp:+fp):*(fp:*tn))):+2:*tn:*log((tn:*nrow):/((tn:+fp):*(tn:+fn))):+2:*fn:*log((fn:*nrow):/((tp:+fn):*(fn:+tn)))
	pearson_c = sqrt((pearson_chi:/nrow):/(1:+pearson_chi:/nrow))
	cramer_v = sqrt(pearson_chi:/nrow)
	cramer_v_bias = sqrt(pearson_chi:/(nrow:*(2:-(1:/(nrow:-1)))))
	
	
	
	



	
	// weighted averages of class-specific metrics
	class_spec_metrics_mat = (precision_k\nega_pred_val_k\recall_k\specificity_k\balanced_acc_k\prevalence_k\false_pos_rate_k\false_alarm_rate_k\false_neg_rate_k\false_omm_rate_k\bias_score_k\pos_lik_ratio_k\neg_lik_ratio_k\youden_j_k\markedness_k\informedness_k\diagnostic_odds_k\yule_q_coeff_k\yule_y_coeff_k\fowlkes_mallows_k\f1_score_k\fb_score_k\matthew_corr_k\threat_k\gilbert_skill_score_k\scott_pi_k\peirce_skill_score_k\cohen_kappa_k\clayton_skill_score_k\extr_dep_score_k\symmetric_extr_dep_score_k\prev_threshold_k\adj_noise_to_signal_k)
	
	weighted_avg = J(33,1,0)
	
	for (i=1; i<=33; i++) {
		for (j=1; j<=ncol; j++) {
			weighted_avg[i] = weighted_avg[i] + class_spec_metrics_mat[i,j]*colsum_conf_mat[j]/nrow
		}
	}


	// brier/logarithmic/spherical score/power score/Psuedospherical score !! Later combine these all in the single loop!!

	//Beta set to 1.5 for now
	power_beta  = 1.5
	pseudo_beta = 1.5

	brier_vec = J(nrow, 1, 0)
	logscore_vec = J(nrow, 1, 0)
	spher_vec = J(nrow, 1, 0)
	power_vec = J(nrow, 1, 0)
	pseudo_vec = J(nrow, 1, 0)


	//Power score
	for (i=1;i<=nrow;i++) {
			power_vec[i] = power_vec[i] + 1/power_beta
	for (j=1;j<=ncol;j++) {
			if (y[i] == j) {
				power_vec[i] = power_vec[i] + (power_beta-1)/power_beta * X[i,j]^power_beta - X[i,j]^(power_beta-1)
			} else {
				power_vec[i] = power_vec[i] + (power_beta-1)/power_beta * X[i,j]^power_beta
			}
		}
	}

	//Pseudospherical score
	for (i=1;i<=nrow;i++) {
		numerator_i = 0
		denominator_i = 0
	for (j=1;j<=ncol;j++) {
			if (y[i] == j) {
				numerator_i = numerator_i + X[i,j]^(pseudo_beta-1)
				denominator_i = denominator_i + X[i,j]^pseudo_beta 
			} else {
				denominator_i = denominator_i + X[i,j]^pseudo_beta 
			}
		}
		pseudo_vec[i] = numerator_i/(denominator_i^((pseudo_beta-1)/pseudo_beta))
	}


	//Brier, log and spherical score
	for (i=1;i<=nrow;i++) {
	for (j=1;j<=ncol;j++) {
			if (y[i] == j) {
				brier_vec[i] = brier_vec[i] + (1 - X[i,j])^2
				logscore_vec[i] = logscore_vec[i] - log(X[i,j])
				spher_vec[i] = spher_vec[i] + X[i,j]/(sqrt((X[i,j])^2 + (1-X[i,j])^2))
			} else {
				brier_vec[i] = brier_vec[i] + (0 - X[i,j])^2
				logscore_vec[i] = logscore_vec[i] - log(1 - X[i,j])
				spher_vec[i] = spher_vec[i] + (1 - X[i,j])/(sqrt( (X[i,j])^2 + (1-X[i,j])^2))
			}
		}
	}

	brier_score = sum(brier_vec)/nrow
	log_score = sum(logscore_vec)/nrow
	spherical_score = sum(spher_vec)/nrow
	power_score = sum(power_vec)/nrow
	pseudo_score = 1 - sum(power_vec)/nrow

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
	printf("\nProbabilistic Forecasts metrics\n")
	printf("Brier score              = {bf:%9.4f} \n", brier_score)
	printf("Power score              = {bf:%9.4f} \n", power_score)
	printf("Logarithmic score        = {bf:%9.4f} \n", log_score)
	printf("Spherical score          = {bf:%9.4f} \n", spherical_score)
	printf("Pseudo spherical score   = {bf:%9.4f} \n", pseudo_score)
 	printf("Ranked probability score = {bf:%9.4f} \n", ranked_probability_score)


	if(ncat > 2){
		// macro measures
		prec_macro = sum(precision_k)/ncol
		recall_macro = sum(recall_k)/ncol
		f1_macro = 2*(prec_macro*recall_macro)/(prec_macro+recall_macro)
		fb_macro = (1+beta^2)*(prec_macro*recall_macro)/(beta^2*(prec_macro+recall_macro))
		fowlkes_mallows_index = sqrt(prec_macro*recall_macro)
	
		// weighted averages (do we want to do this?)
		// also do we want each confusion matrix class specific?

 		displayas("txt")
 		printf("\nMulticlass metrics\n")
 		printf("Accuracy                 = {bf:%9.4f} \n", accuracy)
		printf("Misclassification rate	 = {bf:%9.4f} \n", misclassification_rate)
		printf("Balanced accuracy   	 = {bf:%9.4f} \n", balanced_acc)
		printf("Rogot-Goldberg       	 = {bf:%9.4f} \n", rogot_goldberg)
		printf("Goodman-Kruskal lambda   = {bf:%9.4f} \n", goodman_kruskal_lambda)
		printf("Goodman-Kruskal tau      = {bf:%9.4f} \n", gktm)
		printf("Goodman-Kruskal gamma    = {bf:%9.4f} \n", gkg)
		printf("Correlation (Matthew)    = {bf:%9.4f} \n", corr_matthew)
		printf("Cohen's kappa coefficient= {bf:%9.4f} \n", cohen_kappa)
		printf("Scott's pi coefficient   = {bf:%9.4f} \n", scott_pi)
		printf("Peirce's skill score     = {bf:%9.4f} \n", peirce_skill_score)
		printf("Clayton skill score      = {bf:%9.4f} \n", clayton_skill_score)
		printf("Gerity skill score       = {bf:%9.4f} \n", gsk)
		printf("Theil U                  = {bf:%9.4f} \n", theil_u)
		printf("Tönnies                  = {bf:%9.4f} \n", tonnies)
		printf("Köppen                   = {bf:%9.4f} \n", koppen_ms)
		printf("Adj. Rand index          = {bf:%9.4f} \n", adj_rand_index)
		printf("Kendall tau              = {bf:%9.4f} \n", kendall_tau)
		printf("Pearson Chi              = {bf:%9.4f} \n", pcm)
		printf("Mean square contingency  = {bf:%9.4f} \n", msc)
		printf("Likelihood G             = {bf:%9.4f} \n", lik_g)
		printf("Mean square cont. C      = {bf:%9.4f} \n", c_contingency)
		printf("Cramer V                 = {bf:%9.4f} \n", cramer_v_multi)
    	printf("Cramer V (bias corr.)    = {bf:%9.4f} \n", cramer_v_bias_multi)

		
		printf("\nClass specific metrics\n")
		printrows = 1
	}else{
		printf("\nBinary metrics\n")
		printf("Accuracy                 = {bf:%9.4f} \n", accuracy)
		printf("Misclassification rate	 = {bf:%9.4f} \n", misclassification_rate)
		printrows = 2
	}

	print_vector("Precision                = ", precision_k[printrows::ncat])
	print_vector("Negative predicted value = ", nega_pred_val_k[printrows::ncat])
	print_vector("Recall                   = ", recall_k[printrows::ncat])
	print_vector("Specificity              = ", specificity_k[printrows::ncat])
	print_vector("Balanced accuracy        = ", balanced_acc_k[printrows::ncat])
	print_vector("Prevalence               = ", prevalence_k[printrows::ncat])
	print_vector("False positive rate      = ", false_pos_k[printrows::ncat])
	print_vector("False alarm rate         = ", false_alarm_rate_k[printrows::ncat])
	print_vector("False negative rate      = ", false_neg_rate_k[printrows::ncat])
	print_vector("False omission rate      = ", false_omm_rate_k[printrows::ncat])
	print_vector("Bias score               = ", bias_score_k[printrows::ncat])
	print_vector("Positive likelihood ratio= ", pos_lik_ratio_k[printrows::ncat])
	print_vector("Negative likelihood ratio= ", neg_lik_ratio_k[printrows::ncat])
	print_vector("Youden's J statistic     = ", youden_j_k[printrows::ncat])  
	print_vector("Markedness               = ", markedness_k[printrows::ncat])
	print_vector("Informedness             = ", informedness_k[printrows::ncat])
	print_vector("Diagnostic odds ratio    = ", diagnostic_odds_k[printrows::ncat])
	print_vector("Yule's Q coefficient     = ", yule_q_coeff_k[printrows::ncat])
	print_vector("Yule's Y coefficient     = ", yule_y_coeff_k[printrows::ncat])
	print_vector("Fowlkes-Mallows index    = ", fowlkes_mallows_k[printrows::ncat])
	print_vector("F1-score                 = ", f1_score_k[printrows::ncat])
	print_vector("FB-score	               = ", fb_score_k[printrows::ncat])
	print_vector("Correlation (Matthews)   = ", matthew_corr_k[printrows::ncat])
	print_vector("Threat score             = ", threat_k[printrows::ncat])
	print_vector("Gilbert Skill Score      = ", gilbert_skill_score_k[printrows::ncat])
	print_vector("Scott's pi coefficient   = ", scott_pi_k[printrows::ncat])
	print_vector("Peirce's skill score     = ", peirce_skill_score_k[printrows::ncat])
	print_vector("Cohen's kappa coefficient= ", cohen_kappa_k[printrows::ncat])
	print_vector("Clayton Skill Score      = ", clayton_skill_score_k[printrows::ncat])
	print_vector("Extreme dependency score = ", extr_dep_score_k[printrows::ncat])
	print_vector("Symm. extr. dep. score   = ", symmetric_extr_dep_score_k[printrows::ncat])
	print_vector("Prevalence threshold     = ", prev_threshold_k[printrows::ncat])
	print_vector("Adj. N2S ratio           = ", adj_noise_to_signal_k[printrows::ncat])

	//From here after draft
	print_vector("Consonni-Todeschini #4   = ", consonni_todeschini[printrows::ncat])
	print_vector("Van der Maarel	       = ", van_der_maarel[printrows::ncat])
	print_vector("Anderberg #1             = ", anderberg_one[printrows::ncat])
	print_vector("Anderberg #2             = ", anderberg_two[printrows::ncat])
	print_vector("Benini                   = ", benini[printrows::ncat])
	print_vector("Russel-Rao               = ", russel_rao[printrows::ncat])
	print_vector("Kulczynski #1            = ", kulczynski_one[printrows::ncat])
	print_vector("Kulczynski #2            = ", kulczynski_two[printrows::ncat])
	print_vector("Johnson	               = ", johnson[printrows::ncat])
	print_vector("McConnaughey	       = ", mcconnaughey[printrows::ncat])
	print_vector("Driver-Kroeber           = ", driver_kroeber[printrows::ncat])
	print_vector("Sorensen	               = ", sorensen[printrows::ncat])
	print_vector("Jaccard 3w               = ", jaccard[printrows::ncat])
	print_vector("Rijsbergen effectiveness = ", rijsbergen[printrows::ncat])
	print_vector("Upholt S		       = ", upholt_s[printrows::ncat])
	print_vector("Gini index               = ", gini_index[printrows::ncat])
	print_vector("Modified Gini index      = ", modified_gini[printrows::ncat])
	print_vector("Normalized collacation   = ", norm_coallacation[printrows::ncat])
	print_vector("savage	               = ", savage[printrows::ncat])
	print_vector("Sokal-Sneath #1          = ", sokal_sneath[printrows::ncat])
	print_vector("Hellinger                = ", hellinger[printrows::ncat])
	print_vector("Lance-Williams           = ", lance_williams[printrows::ncat])
	print_vector("Mountford	               = ", mountford[printrows::ncat])
	print_vector("Michelet	               = ", michelet[printrows::ncat])
	print_vector("Sorgenfrei	       = ", sorgenfrei[printrows::ncat])
	print_vector("Fager	               = ", fager[printrows::ncat])
	print_vector("Fager-McGowan            = ", fager_mcgowan[printrows::ncat])
	print_vector("U cost		       = ", u_cost[printrows::ncat])
	print_vector("S cost		       = ", s_cost[printrows::ncat])
	print_vector("R cost		       = ", r_cost[printrows::ncat])
	print_vector("T combined cost	       = ", t_cost[printrows::ncat])
	print_vector("Simpson		       = ", simpson[printrows::ncat])
	print_vector("Kent-Foster #1	       = ", kent_foster_one[printrows::ncat])
	print_vector("Braum-Blanquet	       = ", braum_blanquet[printrows::ncat])
	print_vector("Chord dissimilarity      = ", chord_dissimilarity[printrows::ncat])
	print_vector("Rousseau skill	       = ", rousseau_skill[printrows::ncat])
	print_vector("Anonymous #1	       = ", anonymous_one[printrows::ncat])
	print_vector("Anonymous #2	       = ", anonymous_two[printrows::ncat])
	print_vector("Soergel distance         = ", soergel_distance[printrows::ncat])
	print_vector("Kent-Foster #2	       = ", kent_foster_two[printrows::ncat])
	print_vector("Consonni-Todeschini #1   = ", consonni_todeschini_one[printrows::ncat])
	print_vector("Consonni-Todeschini #2   = ", consonni_todeschini_two[printrows::ncat])
	print_vector("Consonni-Todeschini #3   = ", consonni_todeschini_three[printrows::ncat])
	print_vector("Consonni-Todeschini #5   = ", consonni_todeschini_five[printrows::ncat])
	print_vector("Matching  	       = ", matching[printrows::ncat])
	print_vector("Balance error rate       = ", balance_error_rate[printrows::ncat])
	print_vector("Austin-Colwell	       = ", austin_colwell[printrows::ncat])
	print_vector("Dominance		       = ", dominance[printrows::ncat])
	print_vector("Geometric mean	       = ", geometric_mean[printrows::ncat])
	print_vector("Adjusted geometric mean  = ", adjusted_g_mean[printrows::ncat])
	print_vector("Optimization precision   = ", optimization[printrows::ncat])
	print_vector("Positive LR	       = ", positive_lik[printrows::ncat])
	print_vector("Negative LR	       = ", negative_lik[printrows::ncat])
	print_vector("Anderberg D coefficient  = ", anderberg_d[printrows::ncat])
	print_vector("AMPLE similarity         = ", ample_sim[printrows::ncat])
	print_vector("Baulieu #1  	       = ", baulieu_one[printrows::ncat])
	print_vector("Baulieu #2  	       = ", baulieu_two[printrows::ncat])
	print_vector("Baulieu #3  	       = ", baulieu_three[printrows::ncat])
	print_vector("Baulieu #4  	       = ", baulieu_four[printrows::ncat])
	print_vector("Maron-Kuhns  	       = ", maron_kuhns[printrows::ncat])
	print_vector("Benini (repulsion)       = ", benini_repulsion[printrows::ncat])
	print_vector("Cole #1  	       	       = ", cole_one[printrows::ncat])
	print_vector("Cole #2  	       	       = ", cole_two[printrows::ncat])
	print_vector("Cole #3  	       	       = ", cole_three[printrows::ncat])
	print_vector("Cole #4  	       	       = ", cole_four[printrows::ncat])
	print_vector("Digby	       	       = ", digby[printrows::ncat])
	print_vector("Hawkin-Dotson	       = ", hawkin_dotson[printrows::ncat])
	print_vector("Pattern difference       = ", pattern_diff[printrows::ncat])
	print_vector("Size difference          = ", size_diff[printrows::ncat])
	print_vector("Shape difference 	       = ", shape_diff[printrows::ncat])
	print_vector("Sneath pattern diff.     = ", sneath_patt_diff[printrows::ncat])
	print_vector("Koppen 1870	       = ", koppen[printrows::ncat])
	print_vector("Binary shape diss.       = ", binary_shape[printrows::ncat])
	print_vector("Rogot-Goldberg #1        = ", rogot_goldberg_one[printrows::ncat])
	print_vector("Lacour		       = ", lacour[printrows::ncat])
	print_vector("Peirce #1		       = ", peirce_one[printrows::ncat])
	print_vector("Peirce #2		       = ", peirce_two[printrows::ncat])
	
	// from here page 17 and onwards
	print_vector("Log odds		       = ", log_odds[printrows::ncat])
	print_vector("Goodman association w.   = ", goodman_ass[printrows::ncat])
	print_vector("Goodman association uw.  = ", goodman_un_ass[printrows::ncat])
	print_vector("Goodman-Kruskal lambda   = ", goodman_lambda[printrows::ncat]) // another version? see paper
	print_vector("Goodman-Kruskal #1       = ", goodman_kruskal_1[printrows::ncat])
	print_vector("Goodman-Kruskal #2       = ", goodman_kruskal_2[printrows::ncat])
	print_vector("Goodman-Kruskal #3       = ", goodman_kruskal_3[printrows::ncat])
	print_vector("Adj. F-score             = ", adj_f_score[printrows::ncat])
	print_vector("Pearson MSCC             = ", pearson_mscc[printrows::ncat])
	print_vector("Gower                    = ", gower[printrows::ncat])
	print_vector("Stiles                   = ", stiles[printrows::ncat])
	print_vector("Dennis                   = ", dennis[printrows::ncat])
	print_vector("Loevinger                = ", loevinger[printrows::ncat])
	print_vector("Pearson Heron #2         = ", pearson_heron_2[printrows::ncat])
	print_vector("Forbes D                 = ", forbes_d[printrows::ncat]) // will we also be in forbes after this?
	print_vector("Alroy corr. Forbes       = ", pearson_heron_2[printrows::ncat])
	// page 20 starts here
	print_vector("Forbes #2                = ", forbes_2[printrows::ncat])
	print_vector("Gini association         = ", gini_ass[printrows::ncat])
	print_vector("Fleiss                   = ", fleiss[printrows::ncat])
	print_vector("Kuhns #1                 = ", kuhns_1[printrows::ncat])
	print_vector("Kuhns #2                 = ", kuhns_2[printrows::ncat])
	print_vector("Goodman concomitance     = ", goodman_conc[printrows::ncat])
	print_vector("Rogers-Tanimoto          = ", rogers_tanimoto[printrows::ncat])
	print_vector("Sokal-Sneath #2          = ", sokal_sneath_2[printrows::ncat])
	print_vector("Sokal-Sneath #3          = ", sokal_sneath_3[printrows::ncat])
	print_vector("Sokal-Sneath #4          = ", sokal_sneath_4[printrows::ncat])
	print_vector("Sokal-Sneath #5          = ", sokal_sneath_5[printrows::ncat])
	print_vector("Kocher-Wong              = ", kocher_wang[printrows::ncat])
	print_vector("Faith                    = ", faith[printrows::ncat])
	print_vector("Hamann                   = ", hamann[printrows::ncat])
	print_vector("Hawkin-Dotson            = ", hawkins_dotson[printrows::ncat])
	print_vector("Roux #1                  = ", roux_1[printrows::ncat])
	print_vector("Roux #2                  = ", roux_2[printrows::ncat])
	print_vector("Unigram subtuples        = ", unigram_subtuples[printrows::ncat])
	print_vector("Norm. google distance    = ", norm_google_dist[printrows::ncat])
	print_vector("Eyraud similarity        = ", eyraud_sim_ind[printrows::ncat])
	print_vector("Fossum                   = ", fossum[printrows::ncat])
	print_vector("Baroni-Urbani-Buser #1   = ", baroni_urbani_buser_1[printrows::ncat])
	print_vector("Baroni-Urbani-Buser #2   = ", baroni_urbani_buser_2[printrows::ncat])
	print_vector("Michael                  = ", michael[printrows::ncat])
	print_vector("Clement                  = ", clement[printrows::ncat])
	print_vector("Harris-Lahey             = ", harris_lahey[printrows::ncat])
	print_vector("Kuder-Richardson         = ", kuder_richardson[printrows::ncat])
	print_vector("BC dissimilarity         = ", bc_diss[printrows::ncat])
	print_vector("Sokal distance           = ", sokal_dist[printrows::ncat])
	print_vector("Yule distance            = ", yule_dist[printrows::ncat])
	print_vector("Quetelet                 = ", quetelet[printrows::ncat])
	print_vector("Köppen                   = ", koppen[printrows::ncat])
	print_vector("Heidke skill score       = ", heidke[printrows::ncat])
	print_vector("Shrank skill score       = ", shrank[printrows::ncat])
	print_vector("Tarwid                   = ", tarwid[printrows::ncat])
	print_vector("Relative accuracy (wgtd) = ", weighted_rel_acc[printrows::ncat])
	print_vector("Krippendorf alpha        = ", krippendorf_a[printrows::ncat])
	print_vector("Scott                    = ", scott_coeff[printrows::ncat])
	print_vector("Gilbert-Wells            = ", gilber_wells[printrows::ncat])
	print_vector("Weighted mut. Inf.       = ", weighted_mut_inf[printrows::ncat])
	print_vector("Extremal dep.            = ", extremal_dep_ind[printrows::ncat])
	print_vector("Symmetric extremal dep.  = ", symmetric_extremal[printrows::ncat])
	print_vector("Dunning                  = ", dunning[printrows::ncat])
	print_vector("Prevalence #2            = ", prevalence_2[printrows::ncat])
	print_vector("Tarantula                = ", tarantula[printrows::ncat])
	print_vector("Discriminant power       = ", discriminant_power[printrows::ncat])
	print_vector("Discrimination distance  = ", discr_distance[printrows::ncat])
	print_vector("Pearson Chi              = ", pearson_chi[printrows::ncat])
	print_vector("Pearson Chi (Yates)      = ", pearson_chi_yates[printrows::ncat])
	print_vector("Mean square contingency  = ", mean_sq_contingency[printrows::ncat])
	print_vector("Goodman Kruskal Tau      = ", goodman_kruskal_tau[printrows::ncat])
	print_vector("Likelihood Chi           = ", likelihood_chi[printrows::ncat])
	print_vector("Pearson C                = ", pearson_c[printrows::ncat])
	print_vector("Cramer V                 = ", cramer_v[printrows::ncat])
	print_vector("Cramer V (bias corr.)    = ", cramer_v_bias[printrows::ncat])
	print_vector("Roux #1                  = ", roux_1[printrows::ncat])
	print_vector("Roux #1                  = ", roux_1[printrows::ncat])
	print_vector("Roux #1                  = ", roux_1[printrows::ncat])

	// this is page 29 onwards

	//Untill here

	if(ncat > 2){
		printf("\nMacro averages of class-specific metrics\n")
		printf("F1-score                 = {bf:%9.4f}\n", f1_macro)
		printf("FB-score                 = {bf:%9.4f}\n", fb_macro)
		printf("Fowlkes-Mallows index    = {bf:%9.4f}\n", fowlkes_mallows_index)
		
		printf("\nWeighted averages of class specific metrics\n")
		printf("Precision                = {bf:%9.4f} \n", weighted_avg[1])
		printf("Negative predicted value = {bf:%9.4f} \n", weighted_avg[2])
		printf("Recall                   = {bf:%9.4f} \n", weighted_avg[3])
		printf("Specificity              = {bf:%9.4f} \n", weighted_avg[4])
		printf("Balanced accuracy        = {bf:%9.4f} \n", weighted_avg[5])
		printf("Prevalence               = {bf:%9.4f} \n", weighted_avg[6])
		printf("False positive rate      = {bf:%9.4f} \n", weighted_avg[7])
		printf("False alarm rate         = {bf:%9.4f} \n", weighted_avg[8])
		printf("False negative rate      = {bf:%9.4f} \n", weighted_avg[9])
		printf("False omission rate      = {bf:%9.4f} \n", weighted_avg[10])
		printf("Bias score               = {bf:%9.4f} \n", weighted_avg[11])
		printf("Positive likelihood ratio= {bf:%9.4f} \n", weighted_avg[12])
		printf("Negative likelihood ratio= {bf:%9.4f} \n", weighted_avg[13])
		printf("Youden's J statistic     = {bf:%9.4f} \n", weighted_avg[14])  
		printf("Markedness               = {bf:%9.4f} \n", weighted_avg[15])
		printf("Informedness             = {bf:%9.4f} \n", weighted_avg[16])
		printf("Diagnostic odds ratio    = {bf:%9.4f} \n", weighted_avg[17])
		printf("Yule's Q coefficient     = {bf:%9.4f} \n", weighted_avg[18])
		printf("Yule's Y coefficient     = {bf:%9.4f} \n", weighted_avg[19])
		printf("Fowlkes-Mallows index    = {bf:%9.4f} \n", weighted_avg[20])
		printf("F1-score                 = {bf:%9.4f} \n", weighted_avg[21])
		printf("FB-score                 = {bf:%9.4f} \n", weighted_avg[22])
		printf("Correlation (Matthews)   = {bf:%9.4f} \n", weighted_avg[23])
		printf("Threat score             = {bf:%9.4f} \n", weighted_avg[24])
		printf("Gilbert Skill Score      = {bf:%9.4f} \n", weighted_avg[25])
		printf("Scott's pi coefficient   = {bf:%9.4f} \n", weighted_avg[26])
		printf("Peirce's skill score     = {bf:%9.4f} \n", weighted_avg[27])
		printf("Cohen's kappa coefficient= {bf:%9.4f} \n", weighted_avg[28])
		printf("Clayton Skill Score      = {bf:%9.4f} \n", weighted_avg[29])
		printf("Extreme dependency score = {bf:%9.4f} \n", weighted_avg[30])
		printf("Symm. extr. dep. score   = {bf:%9.4f} \n", weighted_avg[31])
		printf("Prevalence threshold     = {bf:%9.4f} \n", weighted_avg[32])
		printf("Adj. N2S ratio           = {bf:%9.4f} \n", weighted_avg[33])
	}
	
	
}

end