
******************
* Apply catalogue to source data to generate stats on predictive performance (sensitivity, specificity)
******************

****************************************************************************************************
****************************************************************************************************
* For Catalogue Version 2
****************************************************************************************************
****************************************************************************************************

set more off, permanently

* import the v2 catalogue produced by stata (or by R but then you'll need to make sure the heading are the same etc)
clear
import excel "List_of_graded_variants_BDQ-CFZ-INH-DLM-relaxed_25Apr2023.xlsx", sheet("Sheet1") firstrow allstring

rename INITIALCONFIDENCEGRADING Initial_Confidence_Grading
rename FINALCONFIDENCEGRADING Final_Confidence_Grading
rename Additionalgradingcriteriaappl Additionalgradingcriteria
keep drug variant Initial_Confidence_Grading Additionalgradingcriteria Final_Confidence_Grading 

bysort drug variant: gen n=_N
drop if n==2 & Final_Confidence_Grading=="3) Uncertain significance"
drop n

save "Graded_mut_STATA_version2_relaxed.dta", replace

* take all the data we used to derive the catalogue. This is an interim file from my main algorithm (the file of cleaned, merged data that we start with, before the samples with katG315 and rpoB450 S mutations are dropped). 
* RIF genotypes without phenotypes are not yet included and will have to be merged in here. 

* merge in orphan_genotypes:
use "master_data_file_pheno_complete.dta", clear
merge m:1 sample_id drug variant using "orphan_genotypes_maxaf.dta"
drop _m
replace phenotype="D" if phenotype==""
save "master_data_file_pheno_complete_with_orphans.dta", replace


foreach a in 75 25 {
use "master_data_file_pheno_complete_with_orphans.dta", clear

** merge in file with RIF genotypes here. Call them pheno=="D" (for dummy observation) **

keep sample_id drug  variant pheno maxaf effect	tier			
gen het=1 if maxaf<.`a'
replace het=0 if het!=1

* merge in v2 catalogue:
merge m:1 drug variant using "Graded_mut_STATA_version2_relaxed.dta"
drop if _m==2
drop _m

* generate new variable called 'gene'
split variant, p("_")
rename variant1 gene
drop variant2 variant3

*** EXPERT RULES START ***

replace Initial_Confidence_Grading="3) Uncertain significance" if Initial_Confidence_Grading==""

* identify the RRDR:
* then anything else in the RRDR that meets criteria:
gen aa = regexs(0) if(regexm(variant, "[0-9][0-9][0-9]"))
destring aa, replace
gen rrdr=1 if strpos(variant,"rpoB_p.") & aa>=426 & aa<=452 & inlist(effect,"synonymous_variant","missing","initiator_codon_variant","stop_retained_variant")==0
drop aa
* and identify those indels where the first aa in the nomenclature falls outside of the rrdr but the second one falls within it (e.g. rpoB_p.Phe425_Gly426del)
gen aa = regexs(1) if(regexm(variant,"([0-9][0-9][0-9])[a-zA-Z]*$"))
destring aa, replace
replace rrdr=1 if strpos(variant,"rpoB_p.") & aa>=426 & aa<=452 & inlist(effect,"synonymous_variant","missing","initiator_codon_variant","stop_retained_variant")==0
replace rrdr=0 if rrdr!=1
replace Final_Confidence_Grading="2) Assoc w R - Interim" if Initial_Confidence_Grading=="3) Uncertain significance" & rrdr==1 & Additionalgradingcriteria!="Borderline" & Final_Confidence_Grading!="1) Assoc w R"
drop aa

* LoF for BDQ/CFZ-Rv0678+pepQ; CAP-tlyA; DLM-ddn+fbiA+fbiB+fbiC+fgd1+Rv2983; ETH:ethA; INH-katG; PZA-pncA; STM-gid
replace Final_Confidence_Grading="2) Assoc w R - Interim" if Initial_Confidence_Grading=="3) Uncertain significance" & inlist(effect,"feature_ablation","frameshift","start_lost","stop_gained")  & inlist(drug,"Bedaquiline","Clofazimine","Capreomycin","Delamanid","Ethionamide","Isoniazid","Pyrazinamide","Streptomycin") & (inlist(gene,"Rv0678","pepQ","tlyA","ddn","fbiA") | inlist(gene,"fbiB","fbiC","fgd1","Rv2983","ethA","katG","pncA","gid")) & tier==1 &  het!=1 & Final_Confidence_Grading!="1) Assoc w R"


* implement epistasis rules:
*KAN: Ignore eis_c.-14C>T and eis_c.-10G>A if found in combination with LoF mutation in eis at a frequency >975%.
*AMK: Ignore eis_c.-14C>T if found in combination with LoF mutation in eis at a frequency >75%.
*BDQ/CFZ: Ignore Rv0678 mutations if found in combination with LoF mutation in mmpL5 at a frequency >75%.

* for Kanamycin
gen z=1 if inlist(variant,"eis_c.-14C>T","eis_c.-10G>A") & drug=="Kanamycin" &  het!=1
replace z=0 if z!=1
bysort sample_id drug: egen mut=max(z)
drop z 

gen z=1 if inlist(effect,"feature_ablation","frameshift","start_lost","stop_gained") & gene=="eis" & maxa>.75
replace z=0 if z!=1
bysort sample_id drug: egen fs=max(z)

gen epistasis=1 if mut+fs==2
drop z mut fs

* for Amikacin
gen z=1 if inlist(variant,"eis_c.-14C>T") & drug=="Amikacin" &  het!=1
replace z=0 if z!=1
bysort sample_id drug: egen mut=max(z)
drop z 

gen z=1 if inlist(effect,"feature_ablation","frameshift","start_lost","stop_gained") & gene=="eis" & maxa>.75
replace z=0 if z!=1
bysort sample_id drug: egen fs=max(z)

replace epistasis=1 if mut+fs==2
drop z mut fs

* for bdq/clf
gen z=1 if strpos(variant,"Rv0678") & inlist(drug,"Bedaquiline","Clofazimine") & inlist(Final_Confidence_Grading,"1) Assoc w R","2) Assoc w R - Interim") &  het!=1
replace z=0 if z!=1
bysort sample_id drug: egen mut=max(z)
drop z 

gen z=1 if inlist(effect,"feature_ablation","frameshift","start_lost","stop_gained") & gene=="mmpL5" & maxa>.75
replace z=0 if z!=1
bysort sample_id drug: egen fs=max(z)

replace epistasis=1 if mut+fs==2
drop z mut fs

replace epistasis=0 if epistasis!=1

*** EXPERT RULES END ***


** Generate prediction based on AwR. This will be 'group 1'.
gen g1=1 if Final_Confidence_Grading=="1) Assoc w R" & het==0 & epistasis==0 & pheno!="D"
replace g1=0 if g1!=1
bysort sample_id drug: egen group1=max(g1)	

* quantify the number of different variants contributing to this group:
bysort drug variant g1: gen n=_n 
replace n=0 if n!=1
replace n=0 if g1!=1
bysort drug: egen tot_g1=sum(n)
drop n

** Generate prediction based on AwRi only (ie strains don't already have an AwR mutation). This will be 'group 2'.
gen g2=1 if Final_Confidence_Grading=="2) Assoc w R - Interim" & het==0 & epistasis==0 & pheno!="D"
replace g2=0 if g2!=1
bysort sample_id drug: egen group2=max(g2)
replace group2=0 if group1==1

* quantify the number of different variants contributing to this group:
bysort drug variant g2: gen n=_n 
replace n=0 if n!=1
replace n=0 if g2!=1
bysort drug: egen tot_g2=sum(n)

drop g1 g2 n

** Any rif genoR strain that does not have group1 or group2 mutations, but does have a coding mutation for BDQ in Rv0678; for INH in katG; or for PZA in pncA. This will be 'group 3'.
gen grif=1 if drug=="Rifampicin" & inlist(Final_Confidence_Grading,"1) Assoc w R","2) Assoc w R - Interim")
replace grif=1 if drug=="Rifampicin" & rrdr==1 
drop rrdr
replace grif=0 if grif!=1
bysort sample_id: egen Grif=max(grif)
drop grif

gen g3=1 if Grif==1 & inlist(drug,"Bedaquiline","Clofazimine","Isoniazid","Pyrazinamide","Delamanid") & inlist(gene,"Rv0678","katG","pncA","ddn","fbiA","fbiB","fbiC","fgd1","Rv2983") & group1==0 & group2==0 & inlist(effect,"feature_ablation","frameshift","inframe_deletion","inframe_insertion","missense_variant","start_lost","stop_gained") & epistasis==0 & pheno!="D" & inlist(Final_Confidence_Grading,"4) Not assoc w R - Interim","5) Not assoc w R")==0 & tier==1 &  het!=1
drop tier
replace g3=0 if g3!=1
bysort sample_id drug: egen group3=max(g3)
replace group3=0 if group1==1 | group2==1

* quantify the number of different variants contributing to this group:
bysort drug variant g3: gen n=_n 
replace n=0 if n!=1
replace n=0 if g3!=1
bysort drug: egen tot_g3=sum(n)

drop g3 gene n

* reduce data set to one row per sample/drug
keep drug sample_id phenotype variant group1 group2 group3 Grif tot*
bysort sample_id drug: gen n=_n
keep if n==1
drop n sample_id variant

* generate new variables to allow us to loop over rif-s, rif-r, and all samples:
gen rifs=1 if Grif==0
replace rifs=0 if rifs!=1
gen rifr=1 if Grif==1
replace rifr=0 if rifr!=1
gen rifall=1 
drop Grif

* compute the phenotypic R rate for all:
bysort drug phenotype: gen Number_rifall_R=_N
replace Number_rifall_R=0 if phenotype!="R"
bysort drug: egen number_rifall_r=max(Number_rifall_R)
drop Number_rifall_R
gen denom_all=1 if pheno!="D"
replace denom_all=0 if denom_all!=1
bysort drug: egen denominator_rifall=sum(denom_all)
drop denom_all

gen percent_pheno_R=number_rifall_r/denominator_rifall*100

* compute the phenotypic R rate for rif-S strains:
bysort drug phenotype rifs: gen Number_rifs_R=_N
replace Number_rifs_R=0 if phenotype!="R"
replace Number_rifs_R=0 if rifs==0
bysort drug: egen number_rifs_r=max(Number_rifs_R)
drop Number_rifs_R
gen denom_rifs=1 if pheno!="D" & rifs==1
replace denom_rifs=0 if denom_rifs!=1
bysort drug: egen denominator_rifs=sum(denom_rifs)
drop denom_rifs

gen percent_rifs_pheno_R=number_rifs_r/denominator_rifs*100

* compute the phenotypic R rate for rif-R strains:
bysort drug phenotype rifr: gen Number_rifr_R=_N
replace Number_rifr_R=0 if phenotype!="R"
replace Number_rifr_R=0 if rifr==0
bysort drug: egen number_rifr_r=max(Number_rifr_R)
drop Number_rifr_R
gen denom_rifr=1 if pheno!="D" & rifr==1
replace denom_rifr=0 if denom_rifr!=1
bysort drug: egen denominator_rifr=sum(denom_rifr)
drop denom_rifr

gen percent_rifr_pheno_R=number_rifr_r/denominator_rifr*100

preserve

bysort drug: gen n=_n
keep if n==1
keep drug denominator_rifr number_rifr_r denominator_rifs number_rifs_r denominator_rifall number_rifall_r


gen rifall_lb=.
gen rifall_ub=.

gen rifr_lb=.
gen rifr_ub=.

gen rifs_lb=.
gen rifs_ub=.

foreach q in r s all {
	forval i=1/`=_N' {
cii proportions denominator_rif`q'[`i'] number_rif`q'_r[`i']  
	 replace rif`q'_lb=r(lb)*100 in `i'
	 replace rif`q'_ub=r(ub)*100 in `i'
	compress
save "drug_ci_`q'.dta", replace
}
}

foreach q in r s all {
use  "drug_ci_`q'.dta", clear
keep drug rif`q'_lb rif`q'_ub
save "drug_ci_`q'.dta", replace
}

use "drug_ci_r.dta", clear
merge 1:1 drug using "drug_ci_s.dta"
drop _m
merge 1:1 drug using "drug_ci_all.dta"
keep drug rifall_lb rifall_ub rifr_lb rifr_ub rifs_lb rifs_ub
save "drug_ci.dta", replace
restore

merge m:1 drug using "drug_ci.dta"
drop _m

save "prepped_data_for_sens_spec_v2_`a'_relaxed.dta", replace
}

** Three analyses follow:

** 1. Generate the outputs for sens / spec / npv / ppv for maxaf .75 and .25; for rif-s, rif-r, and rif-all; and for groups 1, 2 and 3 in turn :

cap erase "sens_spec_`c(current_date)'_v2_relaxed.xlsx"

foreach a in 75 25 {

foreach y in rifs rifr rifall {
	
foreach q in 1 2 3 {

use "prepped_data_for_sens_spec_v2_`a'_relaxed.dta", clear

* generate tp tn fp fn for each group:
gen TP=1 if pheno=="R" & group`q'==1 & `y'==1
gen TN=1 if pheno=="S" & group`q'==0 & `y'==1
gen FP=1 if pheno=="S" & group`q'==1 & `y'==1
gen FN=1 if pheno=="R" & group`q'==0 & `y'==1


replace TP=0 if TP!=1
replace TN=0 if TN!=1
replace FP=0 if FP!=1
replace FN=0 if FN!=1

bysort drug: egen tp=sum(TP)
bysort drug: egen tn=sum(TN)
bysort drug: egen fp=sum(FP)
bysort drug: egen fn=sum(FN)


drop TP TN FP FN

* generate predictions for each set:

gen sensitivity_group`q'=tp/(tp+fn)
gen specificity_group`q'=tn/(tn+fp)

gen sensitivity_group`q'_lb=.
gen sensitivity_group`q'_ub=.

gen  totR=tp+fn
gen  totS=tn+fp

bysort drug: gen n=_n
keep if n==1
drop pheno 

 forval i=1/`=_N' {
capture	 cii proportions totR[`i'] tp[`i']  
	 replace sensitivity_group`q'_ub=r(ub)*100 in `i'
	 replace sensitivity_group`q'_lb=r(lb)*100 in `i'
	compress
}

gen specificity_group`q'_lb=.
gen specificity_group`q'_ub=.

quietly forval i=1/`=_N' {
 capture   cii proportions totS[`i']  tn[`i'] 
	replace specificity_group`q'_ub=r(ub)*100 in `i'
	replace specificity_group`q'_lb=r(lb)*100 in `i'
	compress
}

gen ppv=tp/(tp + fp)
gen npv=tn/(tn + fn)

gen denominatorS=tn+fn
gen denominatorR=tp+fp

gen npv_group`q'_lb=.
gen npv_group`q'_ub=.

gen ppv_group`q'_lb=.
gen ppv_group`q'_ub=.

quietly forval i=1/`=_N' {
capture	 cii proportions denominatorS[`i']  tn[`i'] 
	replace npv_group`q'_ub=r(ub)*100 in `i'
	replace npv_group`q'_lb=r(lb)*100 in `i'
	compress
}

quietly forval i=1/`=_N' {
capture	 cii proportions denominatorR[`i']  tp[`i'] 
	replace ppv_group`q'_ub=r(ub)*100 in `i'
	replace ppv_group`q'_lb=r(lb)*100 in `i'
	compress
}


export excel using "sens_spec_`c(current_date)'_v2_relaxed.xlsx", firstrow(variables) sheet("group`q'_`a'_`y'") sheetmodify 
}

}
}

** 2. Generate the outputs for sens / spec / npv / ppv for maxaf .75 and .25; for rif-s, rif-r, and rif-all, for groups 1 and 2 combined:


foreach a in 75 25 {

foreach y in rifs rifr rifall {

use "prepped_data_for_sens_spec_v2_`a'_relaxed.dta", clear

* drop the total number of different mutations contributing to each group as only relevant to the analysis immediately above:
drop tot*

* generate tp tn fp fn for each group:
gen TP=1 if pheno=="R" & (group1==1 | group2==1) & `y'==1
gen TN=1 if pheno=="S" & (group1==0 & group2==0) & `y'==1
gen FP=1 if pheno=="S" & (group1==1 | group2==1) & `y'==1
gen FN=1 if pheno=="R" & (group1==0 & group2==0) & `y'==1

replace TP=0 if TP!=1
replace TN=0 if TN!=1
replace FP=0 if FP!=1
replace FN=0 if FN!=1

bysort drug: egen tp=sum(TP)
bysort drug: egen tn=sum(TN)
bysort drug: egen fp=sum(FP)
bysort drug: egen fn=sum(FN)

drop TP TN FP FN

* generate predictions for each set:

gen sensitivity_group12=tp/(tp+fn)
gen specificity_group12=tn/(tn+fp)

gen sensitivity_group12_lb=.
gen sensitivity_group12_ub=.

gen  totR=tp+fn
gen  totS=tn+fp

bysort drug: gen n=_n
keep if n==1
drop pheno 

 forval i=1/`=_N' {
capture	 cii proportions totR[`i'] tp[`i']  
	 replace sensitivity_group12_ub=r(ub)*100 in `i'
	 replace sensitivity_group12_lb=r(lb)*100 in `i'
	compress
}

gen specificity_group12_lb=.
gen specificity_group12_ub=.

quietly forval i=1/`=_N' {
 capture   cii proportions totS[`i']  tn[`i'] 
	replace specificity_group12_ub=r(ub)*100 in `i'
	replace specificity_group12_lb=r(lb)*100 in `i'
	compress
}

gen ppv=tp/(tp + fp)
gen npv=tn/(tn + fn)

gen denominatorS=tn+fn
gen denominatorR=tp+fp

gen npv_group12_lb=.
gen npv_group12_ub=.

gen ppv_group12_lb=.
gen ppv_group12_ub=.

quietly forval i=1/`=_N' {
capture	 cii proportions denominatorS[`i']  tn[`i'] 
	replace npv_group12_ub=r(ub)*100 in `i'
	replace npv_group12_lb=r(lb)*100 in `i'
	compress
}

quietly forval i=1/`=_N' {
capture	 cii proportions denominatorR[`i']  tp[`i'] 
	replace ppv_group12_ub=r(ub)*100 in `i'
	replace ppv_group12_lb=r(lb)*100 in `i'
	compress
}


export excel using "sens_spec_`c(current_date)'_v2_relaxed.xlsx", firstrow(variables) sheet("group12_`a'_`y'") sheetmodify 

}
}

** Generate the outputs for sens / spec / npv / ppv for maxaf .75 and .25; for rif-s, rif-r, and rif-all, for groups 1, 2 and 3 combined:


foreach a in 75 25 {

foreach y in rifs rifr rifall {


use "prepped_data_for_sens_spec_v2_`a'_relaxed.dta", clear

* drop the total number of different mutations contributing to each group as only relevant to the first analysis above:
drop tot*

* generate tp tn fp fn for each group:
gen TP=1 if pheno=="R" & (group1==1 | group2==1 | group3==1) & `y'==1
gen TN=1 if pheno=="S" & (group1==0 & group2==0 & group3==0) & `y'==1
gen FP=1 if pheno=="S" & (group1==1 | group2==1 | group3==1) & `y'==1
gen FN=1 if pheno=="R" & (group1==0 & group2==0 & group3==0) & `y'==1

replace TP=0 if TP!=1
replace TN=0 if TN!=1
replace FP=0 if FP!=1
replace FN=0 if FN!=1

bysort drug: egen tp=sum(TP)
bysort drug: egen tn=sum(TN)
bysort drug: egen fp=sum(FP)
bysort drug: egen fn=sum(FN)

drop TP TN FP FN

* generate predictions for each set:

gen sensitivity_group123=tp/(tp+fn)
gen specificity_group123=tn/(tn+fp)

gen sensitivity_group123_lb=.
gen sensitivity_group123_ub=.

gen  totR=tp+fn
gen  totS=tn+fp

bysort drug: gen n=_n
keep if n==1
drop pheno 

 forval i=1/`=_N' {
capture	 cii proportions totR[`i'] tp[`i']  
	 replace sensitivity_group123_ub=r(ub)*100 in `i'
	 replace sensitivity_group123_lb=r(lb)*100 in `i'
	compress
}

gen specificity_group123_lb=.
gen specificity_group123_ub=.

quietly forval i=1/`=_N' {
capture    cii proportions totS[`i']  tn[`i'] 
	replace specificity_group123_ub=r(ub)*100 in `i'
	replace specificity_group123_lb=r(lb)*100 in `i'
	compress
}

gen ppv=tp/(tp + fp)
gen npv=tn/(tn + fn)

gen denominatorS=tn+fn
gen denominatorR=tp+fp

gen npv_group123_lb=.
gen npv_group123_ub=.

gen ppv_group123_lb=.
gen ppv_group123_ub=.

quietly forval i=1/`=_N' {
capture	 cii proportions denominatorS[`i']  tn[`i'] 
	replace npv_group123_ub=r(ub)*100 in `i'
	replace npv_group123_lb=r(lb)*100 in `i'
	compress
}

quietly forval i=1/`=_N' {
capture	 cii proportions denominatorR[`i']  tp[`i'] 
	replace ppv_group123_ub=r(ub)*100 in `i'
	replace ppv_group123_lb=r(lb)*100 in `i'
	compress
}


export excel using "sens_spec_`c(current_date)'_v2_relaxed.xlsx", firstrow(variables) sheet("group123_`a'_`y'") sheetmodify 

}
}





