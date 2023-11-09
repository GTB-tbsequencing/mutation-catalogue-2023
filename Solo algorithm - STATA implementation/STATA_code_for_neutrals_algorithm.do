****************************
****************************
* Identify neutral mutations in a pre-algorithmic step:
****************************
****************************

clear

set more off, permanently


foreach y in complete noNICD {	

forvalues a = 1(1)2 {

use "master_data_file_WHOa_`y'.dta", clear

gen byte set=`a'

* exclude samples with poor qc
merge m:1 sample_id using "excluded_after_qc.dta"
keep if _m==1
drop _m

* merge in uncontroversial R mutations and mask the samples containing them
merge m:1 drug variant using "urm.dta"
drop if _m==2

split variant, p("_")
gen aaposition=subinstr(variant2,".",">",.)
moss aaposition, match("([0-9]+)") regex
rename _match1 pos1
replace pos1="-"+pos1 if strpos(variant,"-")
destring pos1, replace
drop aaposition _count _pos	variant1 variant2 variant3 _m
* and identify those indels where the first aa in the nomenclature falls outside of the rrdr but the second one falls within it (e.g. rpoB_p.Phe425_Gly426del)
gen pos2 = regexs(1) if(regexm(variant,"([0-9][0-9][0-9])[a-zA-Z]*$"))
destring pos2, replace

replace urm=0 if urm!=1
* RIF
replace urm=1 if strpos(variant,"rpoB") & pos1>=426 & pos1<=452 & effect!="synonymous_variant" & strpos(variant,"-")==0
replace urm=1 if strpos(variant,"rpoB") & pos2>=426 & pos2<=452 & effect!="synonymous_variant" & strpos(variant,"-")==0
drop pos1 pos2
* rpoB borderline mutations:
replace urm=1 if inlist(variant,"rpoB_p.Leu430Pro", "rpoB_p.Asp435Tyr", "rpoB_p.His445Leu", "rpoB_p.His445Asn", "rpoB_p.His445Ser", "rpoB_p.Leu452Pro")
* INH, PZA, SM and ETH
replace urm=1 if strpos(variant,"katG") & inlist(effect,"stop_gained","frameshift","start_lost","feature_ablation") & tier==1
replace urm=1 if strpos(variant,"pncA") & inlist(effect,"stop_gained","frameshift","start_lost","feature_ablation") & tier==1
replace urm=1 if strpos(variant,"gid") & inlist(effect,"stop_gained","frameshift","start_lost","feature_ablation") & tier==1
replace urm=1 if strpos(variant,"ethA") & inlist(effect,"stop_gained","frameshift","start_lost","feature_ablation") & tier==1
replace urm=1 if strpos(variant,"panD") & inlist(effect,"stop_gained","frameshift","start_lost","feature_ablation") & tier==1
replace urm=1 if strpos(variant,"Rv0678") & inlist(effect,"stop_gained","frameshift","start_lost","feature_ablation","inframe_deletion","inframe_insertion","missense_variant","stop_lost") 
replace urm=1 if strpos(variant,"pepQ") & inlist(effect,"stop_gained","frameshift","start_lost","feature_ablation") & tier==1
replace urm=1 if strpos(variant,"Rv1979c") & inlist(effect,"stop_gained","frameshift","start_lost","feature_ablation") & tier==1
replace urm=1 if strpos(variant,"fgd1") & inlist(effect,"stop_gained","frameshift","start_lost","feature_ablation") & tier==1
replace urm=1 if strpos(variant,"ddn") & inlist(effect,"stop_gained","frameshift","start_lost","feature_ablation") & tier==1
replace urm=1 if strpos(variant,"fbiC") & inlist(effect,"stop_gained","frameshift","start_lost","feature_ablation") & tier==1
replace urm=1 if strpos(variant,"Rv2983") & inlist(effect,"stop_gained","frameshift","start_lost","feature_ablation") & tier==1
replace urm=1 if strpos(variant,"fbiA") & inlist(effect,"stop_gained","frameshift","start_lost","feature_ablation") & tier==1
replace urm=1 if strpos(variant,"fbiB") & inlist(effect,"stop_gained","frameshift","start_lost","feature_ablation") & tier==1
replace urm=1 if strpos(variant,"mshA") & inlist(effect,"stop_gained","frameshift","start_lost","feature_ablation") & tier==1
replace urm=1 if strpos(variant,"tlyA") & inlist(effect,"stop_gained","frameshift","start_lost","feature_ablation") & tier==1
replace urm=1 if strpos(variant,"fbiC") & inlist(effect,"stop_gained","frameshift","start_lost","feature_ablation") & tier==1


bysort sample_id drug: egen URM=max(urm)
drop if URM==1 & set==2
drop URM urm 

drop if inlist(variant,"missing","")

* work out ppv:
gen MUTR=1 if pheno=="R"
replace MUTR=0 if pheno=="S"
gen MUTS=1 if pheno=="S"
replace MUTS=0 if pheno=="R"

bysort drug variant: egen mut_R=sum(MUTR)
bysort drug variant: egen mut_S=sum(MUTS)

sort drug variant
egen i=group(drug variant)

gen allmut=mut_S+mut_R

gen ppv=.
gen ppv_lb=.
gen ppv_ub=.

quietly forval i=1/`=_N' {
	cii proportions `=allmut[`i']' `=mut_R[`i']'
	replace ppv=r(proportion) in `i'
	replace ppv_ub=r(ub) in `i'
	replace ppv_lb=r(lb) in `i'
}	

keep if ppv_ub<0.1

save "neutral_mutations_WHO_`a'_`y'.dta", replace

}

use "neutral_mutations_WHO_1_`y'.dta", clear

append using "neutral_mutations_WHO_2_`y'.dta"
keep drug variant set ppv_ub 

bysort drug variant set: gen total=_N
bysort drug variant set: gen n=_n
keep if n==1
drop n

save "neutral_mutations_`y'.dta", replace

drop ppv total
bysort variant drug : gen n=_n
reshape wide set, i(variant drug) j(n)
gen setA=1 if set1==1 | set2==1 
gen setB=1 if set1==2 | set2==2 
drop set1 set2  
replace setA=0 if setA!=1
replace setB=0 if setB!=1

save "neutral_mutations_deduplicated_`y'.dta", replace




***********
* Now identify other neutral mutations as SOLOs, masking setC neutral (setC is A+B+merker):
***********
set more off, permanently

use "master_data_file_WHOa_`y'.dta", clear


* exclude samples with poor qc
merge m:1 sample_id using "excluded_after_qc.dta"
keep if _m==1
drop _m

* mask tier 2 mutations
gen byte mask=1 if tier==2 | inlist(effect,"synonymous_variant","initiator_codon_variant","stop_retained_variant")
replace mask=0 if mask!=1
replace mask=0 if effect=="missing" & tier==1

* label Merker neutral mutations:
merge m:1 drug variant using "merker_neutrals.dta"
replace merker=0 if merker!=1
drop _m

merge m:1 drug variant using "neutral_mutations_deduplicated_`y'.dta"
replace setA=0 if setA==.
replace setB=0 if setB==.
gen setC=1 if setA==1 | setB==1 | merker==1
replace setC=0 if setC!=1
drop if _m==2
drop _m

* Quantify the number of non-masked variants per sample and characterise them according to rules of the algorithm (1st pass)
gen byte MUT=1 if mask==0 & variant!="" & setC==0
replace MUT=0 if MUT!=1

bysort sample_id drug: egen solo=sum(MUT)

replace solo=0 if solo!=1
drop if het==1 | inlist(variant,"","missing") | inlist(effect,"synonymous_variant","initiator_codon_variant","stop_retained_variant")

* identify algorithmic neutral mutations where the ppvSOLOub<10% and they occur as SOLOs =>5 times. This will become set D:
* the formula here is PPV SOLO = (𝑀𝑈𝑇 𝑆𝑂𝐿𝑂 𝑝ℎ𝑒𝑛𝑜 𝑅)/(𝑀𝑈𝑇 𝑆𝑂𝐿𝑂 𝑝ℎ𝑒𝑛𝑜 𝑅 + 𝑀𝑈𝑇 𝑆𝑂𝐿𝑂 𝑝ℎ𝑒𝑛𝑜 𝑆)
gen MUTR=1 if pheno=="R" & solo==1 & MUT==1
replace MUTR=0 if pheno=="S"
gen MUTS=1 if pheno=="S" & solo==1 & MUT==1
replace MUTS=0 if pheno=="R"

bysort drug variant: egen mut_R=sum(MUTR)
bysort drug variant: egen mut_S=sum(MUTS)

sort drug variant
egen i=group(drug variant)
compress

gen allmut=mut_S+mut_R
drop if allmut==0


gen ppvSOLO=.
gen ppvSOLO_lb=.
gen ppvSOLO_ub=.

quietly forval i=1/`=_N' {
	cii proportions `=allmut[`i']' `=mut_R[`i']'
	replace ppvSOLO=r(proportion) in `i'
	replace ppvSOLO_ub=r(ub) in `i'
	replace ppvSOLO_lb=r(lb) in `i'
}

keep if ppvSOLO_ub<0.1

keep drug variant 
bysort drug variant: gen n=_n
keep if n==1
drop n
gen setD1=1

save "WHO_setD1_`y'.dta", replace


*****
*Now repeat this, masking setC and setD1 neutrals identifying neutral mutations by their ppvSOLO:
*****

set more off, permanently

use "master_data_file_WHOa_`y'.dta", clear


* exclude samples with poor qc
merge m:1 sample_id using "excluded_after_qc.dta"
keep if _m==1
drop _m

* mask tier 2 mutations
gen byte mask=1 if tier==2 | inlist(effect,"synonymous_variant","initiator_codon_variant","stop_retained_variant")
replace mask=0 if mask!=1
replace mask=0 if effect=="missing" & tier==1

* label neutral mutations:
merge m:1 drug variant using "merker_neutrals.dta"
replace merker=0 if merker!=1
drop _m

merge m:1 drug variant using "neutral_mutations_deduplicated_`y'.dta"
replace setA=0 if setA==.
replace setB=0 if setB==.
gen setC=1 if setA==1 | setB==1 | merker==1
replace setC=0 if setC!=1
drop if _m==2
drop _m

merge m:1 drug variant using "WHO_setD1_`y'.dta"
drop if _m==2
replace setD1=0 if setD1!=1

drop _m

* Quantify the number of non-masked variants per sample and characterise them according to rules of the algorithm (1st pass)
gen byte MUT=1 if mask==0 & variant!="" & setC==0 & setD1==0
replace MUT=0 if MUT!=1

bysort sample_id drug: egen solo=sum(MUT)

replace solo=0 if solo!=1
drop if het==1 | inlist(variant,"","missing") | inlist(effect,"synonymous_variant","initiator_codon_variant","stop_retained_variant")

* identify algorithmic neutral mutations where the ppvSOLOub<10% and they occur as SOLOs =>5 times. This will become set D:
* the formula here is PPV SOLO = (𝑀𝑈𝑇 𝑆𝑂𝐿𝑂 𝑝ℎ𝑒𝑛𝑜 𝑅)/(𝑀𝑈𝑇 𝑆𝑂𝐿𝑂 𝑝ℎ𝑒𝑛𝑜 𝑅 + 𝑀𝑈𝑇 𝑆𝑂𝐿𝑂 𝑝ℎ𝑒𝑛𝑜 𝑆)
gen MUTR=1 if pheno=="R" & solo==1 & MUT==1
replace MUTR=0 if pheno=="S"
gen MUTS=1 if pheno=="S" & solo==1 & MUT==1
replace MUTS=0 if pheno=="R"

bysort drug variant: egen mut_R=sum(MUTR)
bysort drug variant: egen mut_S=sum(MUTS)

sort drug variant
egen i=group(drug variant)
compress

gen allmut=mut_S+mut_R
drop if allmut==0

gen ppvSOLO=.
gen ppvSOLO_lb=.
gen ppvSOLO_ub=.

quietly forval i=1/`=_N' {
	cii proportions `=allmut[`i']' `=mut_R[`i']'
	replace ppvSOLO=r(proportion) in `i'
	replace ppvSOLO_ub=r(ub) in `i'
	replace ppvSOLO_lb=r(lb) in `i'
}


keep if ppvSOLO_ub<0.1

keep drug variant
bysort drug variant: gen n=_n
keep if n==1
drop n
gen setD2=1

save "WHO_setD2_`y'.dta", replace

use "WHO_setD2_`y'.dta", clear

append using "WHO_setD1_`y'.dta"

append using "neutral_mutations_deduplicated_`y'.dta"
merge 1:1 drug variant using "neutral_mutations_catalogue_v1.dta"
gen v1=1 if _m!=1
drop _m
merge 1:1 drug variant using "merker_neutrals.dta"
order set*, alpha
order variant, first
order drug, first
drop _m
rename merker literature
gen setC=1 if (setB==1 | setA==1 | literature==1)

bysort drug variant: gen n=_n
keep if n==1

foreach b in A B C D1 D2 {
replace set`b'=0 if set`b'==.
}
replace v1=0 if v1==.
replace literature=0 if literature==.

save "detailed_neutral_sets_`y'.dta", replace

keep drug variant

gen neutral="true"
save "algorithmic_neutrals_`y'.dta", replace
}

