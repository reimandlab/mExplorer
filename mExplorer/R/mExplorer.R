#' Selection of process-specific regulators from high-throughput data using multinomial regression models.
#' 
#' @param dframe Data frame of predictors. Row and column names are required for identifying samples (genes) and predictors (gene regulators), respectively.
#' @param response Vector of factors. Names of vector need to correspond to rownames in dframe.
#' @param interactions If enabled, pairs of predictors as interactions will be evaluated (much slower).
#' @param significance Significance cutoff for p-values from log likelihood ratio tests.
#' @param n_cores Number of processor cores to engage in computation. Use all available cores by default (n_cores=0).
#' @param multitest Method to perform multiple testing correction for p-values from predictor evaluation. See p.adjust() for details.

#' @return Vector of scores, with names corresponding to predictors.
#' @references  TBA
#' @author  Juri Reimand <juri.reimand@@utoronto.ca>
#' @examples
#'  data(yeastCCgenes)
#' 	data(yeastTFdata)
#'  mExplorer(yeastTFdata, yeastCCgenes)
#' @export

mExplorer = function(
			dframe, 
			response, 
			interactions = F, 
			significance = 0.05,
			n_cores = 1, # 0 use all available cores, int-- use requested cores
			multitest = "BY"
		) {
	
	if(n_cores==0) {
		n_cores = NULL
	}
	
	# response vector prepared
	response = prepare_response(response, rownames(dframe))
	all_resp_levels = setdiff(levels(response), ".")

	# dframe of predictors prepared
	dframe = prepare_dframe(dframe)
	all_pred_levels = setdiff(unique(c(as.matrix(dframe))), ".")
	
	# collection of submodels prepared
	submodels = colnames(dframe)
	if (interactions) {
		submodels = c(submodels, combn(submodels, 2, paste, collapse=":"))
	}
	
	# fit null model
	fit0 = nnet::multinom(response~1, data=dframe, maxit=5000, trace=0)
	ll_fit0 = logLik(fit0)
	df_fit0 = fit0$edf
	
	# fit alternative models and test against null model
	regressed = parallel::mclapply(1:length(submodels), 
		mExplorer:::run_regression, 
			submodels, dframe, response, ll_fit0, df_fit0, all_pred_levels, all_resp_levels,
		mc.cores=n_cores)
	pval_list = sapply(regressed, '[[', 'chi_p')
	resp_coef_list = do.call("rbind", lapply(regressed, '[[', 'resp_coefs'))
	pred_coef_list = do.call("rbind", lapply(regressed, '[[', 'pred_coefs'))
	names(pval_list) = rownames(resp_coef_list) = rownames(pred_coef_list) = submodels
	
	# collect pvalues from regression into scores
	scores = p.adjust(unlist(pval_list), method=multitest)
	scores = sort(scores[scores<=significance])
	resp_coef_list = resp_coef_list[names(scores),,drop=F]
	pred_coef_list = pred_coef_list[names(scores),,drop=F]
	colnames(resp_coef_list) = all_resp_levels
	colnames(pred_coef_list) = all_pred_levels
	return(list(scores=scores, resp_coefs=resp_coef_list, pred_coefs=pred_coef_list))
}

# pkey=229 
# pkey = 81 # one level in predictor
run_regression = function(pkey, submodels, dframe, response, ll_fit0, df_fit0, all_pred_levels, all_resp_levels) {
	my_formula = as.formula(paste("response ~", submodels[pkey]))
	fit1 = nnet::multinom(my_formula, data=dframe, maxit=5000, trace=0)
	
	ci = tryCatch(confint(fit1), error=function(x) NULL)
	co = coef(fit1)
	ci_sign = resp_specificity = pred_specificity = NULL
	if (length(all_resp_levels)==1) {
		resp_specificity = structure(1, names=all_resp_levels)
		pred_specificity = rep(0, length(co))
		if (!is.null(ci)) {
			ci_sign = sign(ci[,1])==sign(ci[,2])
			co[is.na(ci_sign) | !ci_sign | co<0] = 0
			co["(Intercept)"] = 0
			pred_specificity = co/sum(co)
		}
	} else if (length(levels(dframe[,submodels[pkey]]))<=2) {
		resp_specificity = rep(0, nrow(co))
		if (!is.null(ci)) {
			ci_sign = t(apply(ci, 3, function(x) sign(x[,1])==sign(x[,2])))
			co[is.na(ci_sign) | !ci_sign | co<0] = 0
			co[,"(Intercept)"] = 0
			resp_specificity = (x<-apply(co, 1, sum))/sum(x)
		}
		pred_specificity = structure(1, names=setdiff(levels(dframe[,submodels[pkey]]), "."))
	} else {
		resp_specificity = rep(0, nrow(co))
		pred_specificity = rep(0, ncol(co))
		if (!is.null(ci)) {
			ci_sign = t(apply(ci, 3, function(x) sign(x[,1])==sign(x[,2])))
			co[is.na(ci_sign) | !ci_sign | co<0] = 0
			co[,"(Intercept)"] = 0
			resp_specificity = (x<-apply(co, 1, sum))/sum(x) # wow functional kung-fu 
			pred_specificity = ((x<-apply(co, 2, sum))/sum(x))[-1] # wow functional kung-fu 
		}
	}
	names(pred_specificity) = gsub(submodels[pkey], "", names(pred_specificity))
	pred_specificity = pred_specificity[all_pred_levels]
	# NB! 2x log likelihood equivalent to anova
	# NB! log.p=T option returns natural logarithm of P, converting to log10
	pp = (pchisq(2*(logLik(fit1)-ll_fit0), fit1$edf-df_fit0, lower.tail=F))[1]
	pp = ifelse(any(sign(co)>0), pp, 1)
	cat(pkey, " ")
		
	return(list(chi_p=pp, resp_coefs=resp_specificity, pred_coefs=pred_specificity))
}



prepare_response = function(gene_levels, dframe_genes) {

	if (!any(names(gene_levels) %in% dframe_genes)) {
		stop("MEXPLORER_ERR: No response names recognized in predictors matrix")
	}

	names_matched = length(which(names(gene_levels) %in% dframe_genes))
	cat(paste("MEXPLORER_MSG: Response names recognized: ", names_matched, "/", length(gene_levels), sep=""), fill=T)

	if (length(grep("\\.", gene_levels))) {
		cat("MEXPLORER_MSG: Background level `.' defined in response", fill=T)
	}

	gene_levels = gene_levels[names(gene_levels) %in% dframe_genes]
	
	response = structure(rep(".", length(dframe_genes)), names=dframe_genes)
	response[names(gene_levels)] = gene_levels
	response = factor(response, levels=c(".", setdiff(unique(response), ".")), ordered=T)

	if (length(levels(response))<2) {
		stop("MEXPLORER_ERR: Uninformative response with <2 levels")
	}

	return(response)
}



prepare_dframe = function(dframe) {
	# remove predictors with one distinct level
	suitable_predictors = apply(dframe, 2, function(x) length(unique(x))> 1)
	dframe = dframe[,suitable_predictors, drop=F]
	if (length(which(!suitable_predictors))>0) {
		cat(paste("MEXPLORER_MSG: Ignoring single-level predictors:", names(which(!suitable_predictors))), fill=T)
	}
	if (ncol(dframe)<1) {
		stop("MEXPLORER_ERR: Too few suitable predictors")
	}
	dframe
}


#' Example vector of yeast transcription factors for m:Explorer
#'
#' @docType data
#' @keywords datasets
#' @name yeastCCgenes
#' @usage data(yeastCCgenes)
#' @format A named character vector with 186 elements
NULL

#' Example predictor data for m:Explorer
#'
#' @docType data
#' @keywords datasets
#' @name yeastTFdata
#' @usage data(yeastTFdata)
#' @format A data frame with 6253 observations of 18 variables
NULL
