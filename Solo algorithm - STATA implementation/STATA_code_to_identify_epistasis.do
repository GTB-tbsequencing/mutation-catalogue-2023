
**************
** EPISTASIS
**************

* import catalogue
clear
import excel "List_of_graded_variants.xlsx", sheet("Sheet1") firstrow allstring

save "List_of_graded_variants.dta", replace

* import the final V2 graded catalogue and prep it: 
use "Leonid_list_of_graded_variants_25Apr2023.dta", clear
keep if inlist(drug,"Amikacin","Kanamycin","Bedaquiline")

keep drug variant Initial_Confidence_Grading  Final_Confidence_Grading
save "catalogue_v2_summary.dta", replace

use "master_data_file_pheno_complete.dta", clear
merge m:1 drug variant using "catalogue_v2_summary.dta"

drop if _m==2
drop _m

merge m:1 drug variant using "v1_grades.dta"
drop if _m==2
drop _m
replace Additionalgradingcriteria_v1="" if inlist(Additionalgradingcriteria_v1, "Evidence from ALL dataset only", "Algorithm pass 2", "RRDR", "Indel or premature stop codon (LoF)", "Borderline")
gen Additionalgradingcriteria = Additionalgradingcriteria_v1

*implement the expert rules
* if pooled LoF = AwR, then any LoF "Uncertain" mutation in the same drug_gene should be upgraded to AwR-interim 
replace Final_Confidence_Grading="3) Uncertain significance" if Final_Confidence_Grading==""
gen z=1 if strpos(variant,"_lof") & Final_Confidence_Grading=="1) Assoc w R" 
replace z=0 if z!=1
bysort drug gene: egen zz=max(z)
replace Final_Confidence_Grading="2) Assoc w R - Interim" if Initial_Confidence_Grading=="3) Uncertain significance" & inlist(effect,"feature_ablation","frameshift","start_lost","stop_gained") & strpos(variant,"_lof")==0 & zz==1
drop z zz

* Upgrade to AwR-interim based on "recognized as DR marker" WHO-endorsed assay:
merge m:1 drug variant using "assay_mutations.dta"

replace Final_Confidence_Grading="2) Assoc w R - Interim" if Initial_Confidence_Grading=="3) Uncertain significance" & _m==3
replace Final_Confidence_Grading="2) Assoc w R - Interim" if _m==2

replace Additionalgradingcriteria="Recognized as DR marker through WHO-endorsed assay" if _m==3 & Initial_Confidence_Grading=="3) Uncertain significance"
replace Additionalgradingcriteria="Recognized as DR marker through WHO-endorsed assay" if _m==2
drop _m

* Upgrade to AwR-interim based on "allelic exchange" evidence
merge m:1 drug variant using "allelic_exchanges.dta"

replace Final_Confidence_Grading="2) Assoc w R - Interim" if Initial_Confidence_Grading=="3) Uncertain significance" & _m==3 
replace Final_Confidence_Grading="2) Assoc w R - Interim" if _m==2 

replace Additionalgradingcriteria="Evidence from allelic exchange experiments" if _m==3 & Initial_Confidence_Grading=="3) Uncertain significance" 
replace Additionalgradingcriteria="Evidence from allelic exchange experiments" if _m==2
drop _m

* Upgrade to AwR-interim based on previous guidance to be distinguished between "evidence in ver1" and "guidance before ver1" (e.g. Miotto ERJ2017 and WHO docs, already embedded in ver1)
replace Final_Confidence_Grading="2) Assoc w R - Interim" if Initial_Confidence_Grading=="3) Uncertain significance" & inlist(Additionalgradingcriteria_v1,"Previous WHO guidance","WHO-endorsed gDST assay")
replace Additionalgradingcriteria="WHO guidance before catalogue version 1" if  Initial_Confidence_Grading=="3) Uncertain significance" & inlist(Additionalgradingcriteria_v1,"Previous WHO guidance","WHO-endorsed gDST assay")

replace Final_Confidence_Grading="2) Assoc w R - Interim" if Initial_Confidence_Grading=="3) Uncertain significance" & inlist(final_grading_v1,"1) Assoc w R","2) Assoc w R - Interim") & Additionalgradingcriteria!="WHO guidance before catalogue version 1"
replace Additionalgradingcriteria="Evidence from catalogue version 1" if Initial_Confidence_Grading=="3) Uncertain significance" & inlist(final_grading_v1,"1) Assoc w R","2) Assoc w R - Interim") & Additionalgradingcriteria!="WHO guidance before catalogue version 1"

keep drug sample_id gene maxaf phenotype variant Final_Confidence_Grading effect

save "epistasis_input.dta", replace

************************************************
* Epistasis calculation 1: eis_c.-14C>T (Amikacin)
************************************************

use "epistasis_input.dta", clear

* Exclude any strain with at least one of the following mutations at ANY mutation frequency:
gen q=1 if inlist(variant,"rrs_n.1401A>G","rrs_n.1402C>T","rrs_n.1484G>T")
replace q=0 if q!=1
bysort sample_id drug: egen qq=max(q)
drop if qq==1
drop q qq

* Compute AMK PPV with binomial 95% CI for any strain with eis_c.-14C>T at a frequency >90% with
* a) LoF mutation in eis at a frequency >90% vs.
* b) No coding mutations in eis at ANY frequency, except synonymous_variant or group 4/5

keep if drug=="Amikacin"
gen q=1 if variant=="eis_c.-14C>T" & maxaf>.9
replace q=0 if q!=1
bysort sample_id: egen qq=max(q)
keep if qq==1
drop q

* identify the strains with relevant lof mutation:
gen LOF=1 if inlist(effect,"feature_ablation","frameshift","start_lost","stop_gained") & gene=="eis" & maxaf>.9
replace LOF=0 if LOF!=1
bysort sample_id: egen lof=max(LOF)
drop LOF

* identify the strains with no coding mutations (CM):
gen CM=0 if gene=="eis" & inlist(effect,"synonymous_variant","stop_retained_variant","initiator_codon_variant","upstream_gene_variant")
replace CM=0 if gene=="eis" & inlist(Final_Confidence_Grading,"4) Not assoc w R - Interim","5) Not assoc w R")
replace CM=1 if CM!=0 & gene=="eis"
replace CM=0 if CM!=1
bysort sample_id: egen cm=max(CM)
drop CM
gen ncm=1 if cm==0
replace ncm=0 if ncm!=1
drop cm

* make predictions through Final_Confidence_Grading and compute PPV
gen r=1 if inlist(Final_Confidence_Grading,"1) Assoc w R","2) Assoc w R - Interim")
replace r=0 if r!=1
bysort sample_id: egen R=max(r)
drop r

gen TP=1 if pheno=="R" & R==1
gen FP=1 if pheno=="S" & R==1

replace TP=0 if TP!=1
replace FP=0 if FP!=1

keep sample_id drug TP FP lof ncm
bysort sample_id drug: gen n=_n
keep if n==1
drop n
keep drug TP FP lof ncm 
save "epistasis_AMI.dta", replace

cap erase "epistasis_`c(current_date)'.xlsx"

foreach a in lof ncm {

use "epistasis_AMI.dta", clear

keep if `a'==1

bysort drug `a': egen tp_`a'=sum(TP)
bysort drug `a': egen fp_`a'=sum(FP)

bysort drug `a': gen n=_n
keep if n==1

* generate predictions for each set:

gen	ppv_`a'=tp/(tp+fp) if `a'==1

gen denominator=tp+fp if `a'==1

gen ppv_`a'_lb=. if `a'==1
gen ppv_`a'_ub=. if `a'==1
	
quietly forval i=1/`=_N' {
capture	 cii proportions denominator[`i']  tp[`i'] 
	replace ppv_`a'_ub=r(ub) in `i'
	replace ppv_`a'_lb=r(lb) in `i'
	compress
}

keep drug `a' tp fp ppv* 

export excel using "epistasis_`c(current_date)'.xlsx", firstrow(variables) sheet("AMI_eis14CT_`a'") sheetmodify 
	
}


************************************************
* Epistasis calculation eis variants for Kanamycin:
************************************************


foreach z in eis_c.-14C>T eis_c.-12C>T eis_c.-37G>T eis_c.-10G>A eis_c.-8delC {

use "epistasis_input.dta", clear
	
* Exclude any strain with at least one of the following mutations at ANY mutation frequency:
gen q=1 if inlist(variant,"rrs_n.1401A>G","rrs_n.1402C>T","rrs_n.1484G>T","eis_c.-14C>T","eis_c.-10G>A","eis_c.-12C>T","eis_c.-37G>T","eis_c.-8delC")
replace q=0 if q!=1
replace q=0 if variant=="`z'"
bysort sample_id drug: egen qq=max(q)
drop if qq==1
drop q qq


* Compute AMK PPV with binomial 95% CI for any strain with eis_c.-14C>T at a frequency >90% with
* a) LoF mutation in eis at a frequency >90% vs.
* b) No coding mutations in eis at ANY frequency, except synonymous_variant or group 4/5

keep if drug=="Kanamycin"
gen q=1 if variant=="`z'" & maxaf>.9
replace q=0 if q!=1
bysort sample_id: egen qq=max(q)
keep if qq==1
drop q

* identify the strains with relevant lof mutation:
gen LOF=1 if inlist(effect,"feature_ablation","frameshift","start_lost","stop_gained") & gene=="eis" & maxaf>.9
replace LOF=0 if LOF!=1
bysort sample_id: egen lof=max(LOF)
drop LOF

* identify the strains with no coding mutations (CM):
gen CM=0 if gene=="eis" & inlist(effect,"synonymous_variant","stop_retained_variant","initiator_codon_variant","upstream_gene_variant")
replace CM=0 if gene=="eis" & inlist(Final_Confidence_Grading,"4) Not assoc w R - Interim","5) Not assoc w R")
replace CM=1 if CM!=0 & gene=="eis"
replace CM=0 if CM!=1
bysort sample_id: egen cm=max(CM)
drop CM
gen ncm=1 if cm==0
replace ncm=0 if ncm!=1
drop cm

* make predictions through Final_Confidence_Grading and compute PPV
gen r=1 if inlist(Final_Confidence_Grading,"1) Assoc w R","2) Assoc w R - Interim")
replace r=0 if r!=1
bysort sample_id: egen R=max(r)
drop r

gen TP=1 if pheno=="R" & R==1
gen FP=1 if pheno=="S" & R==1

replace TP=0 if TP!=1
replace FP=0 if FP!=1

keep sample_id drug TP FP lof ncm
bysort sample_id drug: gen n=_n
keep if n==1
drop n
keep drug TP FP lof ncm 
save "epistasis_KAN_`z'.dta", replace

}

foreach z in eis_c.-14C>T eis_c.-10G>A {

foreach a in lof ncm {

use "epistasis_KAN_`z'.dta", clear

keep if `a'==1

bysort drug `a': egen tp_`a'=sum(TP)
bysort drug `a': egen fp_`a'=sum(FP)

bysort drug `a': gen n=_n
keep if n==1

* generate predictions for each set:

gen	ppv_`a'=tp/(tp+fp) if `a'==1

gen denominator=tp+fp if `a'==1

gen ppv_`a'_lb=. if `a'==1
gen ppv_`a'_ub=. if `a'==1
	
quietly forval i=1/`=_N' {
	cii proportions denominator[`i']  tp[`i'] 
	replace ppv_`a'_ub=r(ub) in `i'
	replace ppv_`a'_lb=r(lb) in `i'
	compress
}

keep drug `a' tp fp ppv* 

export excel using "epistasis_`c(current_date)'.xlsx", firstrow(variables) sheet("KAN_`z'_`a'") sheetmodify 
	
}
}

* for eis_c.-12C>T eis_c.-37G>T there are no associated lof mutations

foreach z in eis_c.-12C>T eis_c.-37G>T eis_c.-8delC {

use "epistasis_KAN_`z'.dta", clear

bysort drug ncm: egen tp_ncm=sum(TP)
bysort drug ncm: egen fp_ncm=sum(FP)

bysort drug ncm: gen n=_n
keep if n==1
keep if ncm==1

* generate predictions for each set:

gen	ppv_ncm=tp/(tp+fp) if ncm==1

gen denominator=tp+fp if ncm==1

gen ppv_ncm_lb=. if ncm==1
gen ppv_ncm_ub=. if ncm==1
	
quietly forval i=1/`=_N' {
	cii proportions denominator[`i']  tp[`i'] 
	replace ppv_ncm_ub=r(ub) in `i'
	replace ppv_ncm_lb=r(lb) in `i'
	compress
}

keep drug ncm tp fp ppv* 

export excel using "epistasis_`c(current_date)'.xlsx", firstrow(variables) sheet("KAN_`z'_ncm") sheetmodify 
	
}

*************************************************************
* Epistasis calculation mmpL5/mmpS5 variants for Bedaquiline:
*************************************************************

* there are no strains with mmpS5 lof mutation at maxaf >90 AND an Rv0678 group1/2 mutation, so only perform this analysis for mmpL5

foreach k in L5 {
	
use "epistasis_input.dta", clear

keep if drug=="Bedaquiline"

* apply lof rule to pepQ
* if pooled LoF = AwR, then any LoF "Uncertain" mutation in the same drug_gene should be upgraded to AwR-interim 
gen z=1 if strpos(variant,"_lof") & Final_Confidence_Grading=="1) Assoc w R" 
replace z=0 if z!=1
bysort drug gene: egen zz=max(z)
replace Final_Confidence_Grading="2) Assoc w R - Interim" if Final_Confidence_Grading=="3) Uncertain significance" & inlist(effect,"feature_ablation","frameshift","start_lost","stop_gained") & strpos(variant,"_lof")==0 & zz==1 & gene=="pepQ"
drop z zz

* Exclude any strain with at least one of the following mutations at ANY mutation frequency, including in promotors:
gen q=1 if inlist(Final_Confidence_Grading,"1) Assoc w R","2) Assoc w R - Interim") & inlist(gene,"atpE","pepQ")
replace q=0 if q!=1
bysort sample_id drug: egen qq=max(q)
drop if qq==1
drop q qq

* Compute BDQ PPV with binomial 95% CI for any strain with and group1/2 Rv0678 mutation at a frequency >90% with
* a) LoF mutation in mmpL5 at a frequency >90% vs.
* b) No coding mutations in mmpL5 at ANY frequency, except synonymous_variant or group 4/5

* identify strains with Rv0678 group1/2 mutations
gen z=1 if gene=="Rv0678" & inlist(Final_Confidence_Grading,"1) Assoc w R","2) Assoc w R - Interim") & maxaf>.9
replace z=0 if z!=1
bysort sample_id: egen zz=max(z)
keep if zz==1
drop z zz

* among these, identify the strains with lof mutations in mmpL5 or mmpS5 with maxaf>90%:
gen LOF=1 if inlist(effect,"feature_ablation","frameshift","start_lost","stop_gained") & gene=="mmp`k'" & maxaf>.9
replace LOF=0 if LOF!=1
bysort sample_id: egen lof=max(LOF)
drop LOF

* identify the strains with no coding mutations (CM) in mmpL5 or mmpS5 at any frequency except synon or grou4/5:
gen CM=0 if gene=="mmp`k'" & inlist(effect,"synonymous_variant","stop_retained_variant","initiator_codon_variant","upstream_gene_variant")
replace CM=0 if inlist(Final_Confidence_Grading,"4) Not assoc w R - Interim","5) Not assoc w R") | gene!="mmp`k'"
replace CM=1 if CM!=0 & gene=="mmp`k'"
bysort sample_id gene: egen cm=max(CM)
drop CM
bysort sample_id: egen cm1=max(cm)
gen ncm=1 if cm1==0
replace ncm=0 if ncm!=1
drop cm cm1

* All samples are geno R by definition already, as we are only looking at samples ith Rv0678 and a group1/2 mutation in them:

gen TP=1 if pheno=="R" 
gen FP=1 if pheno=="S" 

replace TP=0 if TP!=1
replace FP=0 if FP!=1

keep sample_id drug TP FP lof ncm
bysort sample_id drug: gen n=_n
keep if n==1
drop n
keep drug TP FP lof ncm 
save "epistasis_BDQ_mmp`k'.dta", replace

foreach a in lof ncm {

use "epistasis_BDQ_mmp`k'.dta", clear

keep if `a'==1

bysort drug `a': egen tp_`a'=sum(TP)
bysort drug `a': egen fp_`a'=sum(FP)

bysort drug `a': gen n=_n
keep if n==1

* generate predictions for each of lof and ncm:

gen	ppv_`a'=tp/(tp+fp) if `a'==1

gen denominator=tp+fp if `a'==1

gen ppv_`a'_lb=. if `a'==1
gen ppv_`a'_ub=. if `a'==1
	
quietly forval i=1/`=_N' {
capture	 cii proportions denominator[`i']  tp[`i'] 
	replace ppv_`a'_ub=r(ub) in `i'
	replace ppv_`a'_lb=r(lb) in `i'
	compress
}

keep drug `a' tp fp ppv* 

export excel using "epistasis_`c(current_date)'.xlsx", firstrow(variables) sheet("BDQ_`k'_`a'") sheetmodify 
	
}


}
