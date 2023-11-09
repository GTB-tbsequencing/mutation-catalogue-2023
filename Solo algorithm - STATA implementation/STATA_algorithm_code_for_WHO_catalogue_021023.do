*0* Start by running the conversion script to prepare the relevant files (LC)
********************************************************************************
run STATA_code_to_convert_old_to_new_nomenclature.do
********************************************************************************

cd FILE_CONTAINING_INPUT_DATA/full_genotypes/

set more off, permanently

*1* For each drug, import genotypic and phenotypic data:
********************************************************************************
********************************************************************************

* Import genotypic data (some of the paths below may need to be changed depending on the structure of the files you're working with):

clear

foreach d in Amikacin Bedaquiline Capreomycin  Delamanid Ethambutol Ethionamide Isoniazid Kanamycin Levofloxacin Linezolid Moxifloxacin Pyrazinamide Rifampicin  Streptomycin Clofazimine Pretomanid {

cd "drug_name=`d'/"

forvalues i = 1/2 {
filelist, dir("tier=`i'")
local dir = dirname
local oldname = filename 
local newname = "`oldname'.csv"
capture _renamefile "`dir'/`oldname'" "`dir'/`newname'"		 
if _rc!=0 di
capture import delimited using "`dir'/`newname'", clear
if _rc!=0 di
capture save "../../../../data_for_analysis/`d'_t`i'.dta", replace
if _rc!=0 di

}

cd ..

}

* and import phenotypic data, categorised as 'pheno' (which means WHO and ALL datasets), and then as CC and CCATU subsets:

clear

cd ../phenotypes/

foreach d in Amikacin Bedaquiline Capreomycin  Delamanid Ethambutol Ethionamide Isoniazid Kanamycin Levofloxacin Linezolid Moxifloxacin Pyrazinamide Rifampicin  Streptomycin Clofazimine Pretomanid {

cd "drug_name=`d'/"

filelist, dir("")

local dir = dirname
local oldname = filename 
local newname = "`oldname'.csv"
_renamefile "`dir'/`oldname'" "`dir'/`newname'"		 
import delimited using "`dir'/`newname'", clear

preserve
keep if inlist(phenotypic_category,"WHO","ALL")
save "../../../../data_for_analysis/`d'_pheno.dta", replace
restore

preserve
keep if phenotypic_category=="CC"
save "../../../../data_for_analysis/`d'_CC.dta", replace
restore

keep if phenotypic_category=="CC-ATU"
save "../../../../data_for_analysis/`d'_CCATU.dta", replace

cd ..

}

*2* Combine all the genotypic and phenotypic data in one big file:
********************************************************************************
********************************************************************************

* first merge each geno t1 and t2 file with its corresponding pheno file, and then append everything into one master_data_file:
 
cd ../../../data_for_analysis/

set more off, permanently

cap erase "master_data_file_pheno.dta"
cap erase "master_data_file_CC.dta"
cap erase "master_data_file_CCATU.dta"

foreach d in Amikacin  Bedaquiline Capreomycin  Delamanid Ethambutol Ethionamide Isoniazid Kanamycin Levofloxacin  Linezolid Moxifloxacin Pyrazinamide Rifampicin  Streptomycin Clofazimine Pretomanid {

foreach a in pheno CC CCATU {

forvalues i = 1/2 {

capture use "`d'_t`i'.dta", clear

if _rc!=0 di

capture replace variant_category="missing" if variant_category==""
capture replace predicted_effect="missing" if predicted_effect==""

capture merge m:1 sample_id using "`d'_`a'.dta"
if _N {
	capture drop if _m==1
	capture drop _m
}

capture gen drug="`d'" 
capture gen byte tier=`i' if resolved_symbol!=""
capture save "`d'_t`i'_geno_`a'.dta", replace
if _rc!=0 di

cap confirm file "master_data_file_`a'.dta"
if _rc!=0 save "master_data_file_`a'.dta", replace
else if _rc==0 {

append using "master_data_file_`a'.dta", force
compress

* drop rows with missing genomic data for isolates that have genomic data for that drug in another row:
gsort sample_id drug - resolved_symbol 
bysort sample_id drug: gen n=_n
drop if resolved_symbol=="" & n>1 
drop n

save "master_data_file_`a'.dta", replace

}
}
}
}

foreach a in pheno CC CCATU { 

use "master_data_file_`a'.dta", clear
 
* rename some variables with a friendlier name:
compress
gen variant=resolved_symbol+"_"+variant_category if resolved_symbol!=""
replace variant="missing" if strpos(variant,"missing")
rename resolved_symbol gene
rename variant_category mutation
rename phenotypic_category category_phenotype
rename predicted_effect effect
drop neutral box

save "master_data_file_`a'_complete.dta", replace
}

* generate additional versions of the data without the BDQ samples from NICD South Africa:

foreach a in pheno CC CCATU { 

use "master_data_file_`a'.dta", clear

merge m:1 sample_id using "nicd_sampleids.dta"
drop if _m==3 & drug=="Bedaquiline"
drop if _m==2
drop _m

* rename some variables to make them friendlier:
compress
gen variant=resolved_symbol+"_"+variant_category if resolved_symbol!=""
replace variant="missing" if strpos(variant,"missing")
rename resolved_symbol gene
rename variant_category mutation
rename phenotypic_category category_phenotype
rename predicted_effect effect
drop neutral box

save "master_data_file_`a'_noNICD.dta", replace
}


 

*3*  Identify which whole isolates to exclude based on katG315 or rpoB450 mutations and an S phenotype as very likely mislabelled. 
********************************************************************************
********************************************************************************
use "master_data_file_pheno_complete.dta", clear

gen excl=1 if variant=="katG_p.Ser315Thr" & pheno=="S" & drug=="Isoniazid"
replace excl=1 if variant=="rpoB_p.Ser450Leu" & pheno=="S" & drug=="Rifampicin"
replace excl=0 if excl!=1
bysort sample_id: egen exclude=max(excl)
keep if exclude==1
keep sample_id
duplicates drop

save "excluded_after_qc.dta", replace


*4* Start algorithm
********************************************************************************
********************************************************************************

foreach y in complete noNICD {

use "master_data_file_pheno_`y'.dta", clear
keep if category_phenotype=="WHO"

* define a het call as anything with a maxaf<0.75:
gen byte het=1 if maxaf<.75
replace het=0 if maxaf>=.75
replace het=. if variant==""

save "master_data_file_WHOa_`y'.dta", replace
}


*** run neutral algorithm
******************************************************************************** 
run ../STATA_code_for_neutrals_algorithm.do
******************************************************************************** 

foreach y in complete noNICD {

use "master_data_file_WHOa_`y'.dta", clear
merge m:1 drug variant using "algorithmic_neutrals_`y'.dta"
drop if _m==2
drop _m 


*Identify variants that we don't plan to count when identifying SOLOs, and mask them:
gen byte mask=1 if neutral=="true" | inlist(effect,"synonymous_variant","initiator_codon_variant","stop_retained_variant")
replace mask=0 if mask!=1
replace mask=0 if effect=="missing" 
* As we're focussing on tier 1 first, we mask tier 2 variants for now too
replace mask=1 if tier==2

order sample_id variant phenotype mask  maxaf

save "master_data_file_WHO_`y'.dta", replace

use "master_data_file_pheno_`y'.dta", clear
replace category_phenotype="ALL"

* merge in neutral calls:
merge m:1 drug variant using "algorithmic_neutrals_`y'.dta"
drop if _m==2
drop _m 


*Identify variants that we don't plan to count when identifying SOLOs, and mask them:
gen byte mask=1 if neutral=="true" | inlist(effect,"synonymous_variant","initiator_codon_variant","stop_retained_variant")
replace mask=0 if mask!=1
replace mask=0 if effect=="missing" 
* As we're focussing on tier 1 first, we mask tier 2 variants for now too
replace mask=1 if tier==2

gen byte het=1 if maxaf<.75
replace het=0 if maxaf>=.75
replace het=. if variant==""
replace het=1 if variant=="missing"

order sample_id variant phenotype mask  maxaf

save "master_data_file_ALL_`y'.dta", replace

foreach a in CC CCATU {

use "master_data_file_`a'_`y'.dta", clear

* merge in neutral calls:
merge m:1 drug variant using "algorithmic_neutrals_`y'.dta"
drop if _m==2
drop _m 

*Identify variants that we don't plan to count when identifying SOLOs, and mask them:
gen byte mask=1 if neutral=="true" | inlist(effect,"synonymous_variant","initiator_codon_variant","stop_retained_variant")
replace mask=0 if mask!=1
replace mask=0 if effect=="missing" 
* As we're focussing on tier 1 first, we mask tier 2 variants for now too
replace mask=1 if tier==2

gen byte het=1 if maxaf<.75
replace het=0 if maxaf>=.75
replace het=. if variant==""

order sample_id variant phenotype mask  maxaf

save "master_data_file_`a'_`y'.dta", replace
}



* generate a loop to first go over just WHO phenoptypes, and then go over all phenotypes, and then the CC and CCATU phenotypes.
foreach c in WHO ALL CC CCATU {

use "master_data_file_`c'_`y'.dta", clear
		
* exclude samples with poor qc, as defined above
merge m:1 sample_id using "excluded_after_qc.dta"
keep if _m==1
drop _m	
		
* Quantify the number of non-masked variants per sample and characterise them according to rules of the algorithm (1st pass)
gen byte MUT=1 if mask==0 & variant!=""
replace MUT=0 if MUT!=1

bysort sample_id drug: egen mut=sum(MUT)

gen SOLO=1 if mut==1 & MUT==1 & het==0 & variant!="missing"
replace SOLO=0 if SOLO!=1

gen byte res=1 if mut==1 & MUT==1 & pheno=="R" & het==0 & variant!="missing"
replace res=0 if res!=1
bysort drug variant: egen Res=max(res)
gen character="R" if Res==1 
drop res Res
replace charact="S" if mask==1

* determine the character to be S if always S
bysort drug variant pheno het: gen n=_N
bysort drug variant het: gen m=_N
replace charac="aS" if n==m & pheno=="S" & character=="" & mask==0 & het==0 & variant!="missing" // where 'aS' mean 'algorithmically determined S'

* determine the character to be S if always S as a solo
gen byte s1=1 if MUT==1 & mut==1 & pheno=="S" & het==0 & variant!="missing"
replace s1=0 if s1!=1
gen byte s2=1 if MUT==1 & mut==1 & het==0 & variant!="missing"
replace s2=0 if s2!=1
bysort drug variant het: egen S1=sum(s1)
bysort drug variant het: egen S2=sum(s2)
compress
replace character="aS" if S1==S2 & S1!=0 & character=="" & pheno=="S" & mask==0 & het==0 & variant!="missing"

* anything identified as aS by either rule ('always S' or 'always S when SOLO') as a variant can also be inferred to be S where maxaf<.75, so change these to 'aS' as well: 
gen byte t=1 if character=="aS"
replace t=0 if t!=1
bysort drug variant: egen tt=max(t)
replace character="aS" if tt==1 
drop t tt S1 S2 s1 s2 mut MUT n m

* characterised what's left as 'U'
replace character="U" if character=="" & tier==1 & mask==0 & het==0 & variant!="missing"

* if a variant is U then it should also be a U as a het
gen byte t=1 if character=="U"
replace t=0 if t!=1
bysort drug variant: egen tt=max(t)
replace character="U" if tt==1 
drop t tt

* lable variants that only ever appear as hets as 'U' unless already 'S'
bysort drug variant: gen n=_N
bysort drug variant het: gen m=_N
replace character="U" if het==1 & n==m & tier==1 & character!="S" & variant!="missing"
drop n m

* annotate the characterisation as 'pass 1'
gen byte pass=1 if inlist(character,"R","aS","U")
replace pass=0 if mask==1
compress

*5* Conduct a 2nd pass of the algorithm where we mask 'aS' variants
********************************************************************************
********************************************************************************
replace mask=1 if character=="aS"

* identify new R's:
gen byte MUT=1 if variant!="" & mask==0
replace MUT=0 if MUT!=1
bysort sample_id drug: egen mut=sum(MUT)
compress

gen byte SOLO1=1 if mut==1 & MUT==1 & het==0 & variant!="missing" & character!="R"
replace SOLO1=0 if character=="R" & mut==1 & MUT==1 & het==0 & variant!="missing" & SOLO==0
replace SOLO=1 if SOLO1==1
drop SOLO1

gen byte res=1 if mut==1 & MUT==1 & pheno=="R" & het==0 & variant!="missing"
replace res=0 if res!=1
bysort drug variant: egen Res=max(res)

gen character2="R" if Res==1 
replace pass=2 if character2=="R" & character=="U"
replace character="R" if character2=="R" & inlist(character,"U","aS")

drop character2 res Res 
compress

* Identify new 'aS' based on the rule that 'S if always S as a solo' (those variants that are 'always and only in S isolates' don't need a second pass as will all have been identified in the 1st pass)
gen byte s1=1 if MUT==1 & mut==1 & pheno=="S" & het==0 & variant!="missing"
replace s1=0 if s1!=1
gen byte s2=1 if MUT==1 & mut==1 & het==0 & variant!="missing"
replace s2=0 if s2!=1
bysort drug variant: egen S1=sum(s1)
bysort drug variant: egen S2=sum(s2)
compress
replace pass=2 if S1==S2 & S1!=0 & pheno=="S" & character!="aS" & variant!="missing"
replace character="aS" if S1==S2 & S1!=0 & pheno=="S" & variant!="missing" & het==0 & mask==0
* again, anything identified as 'aS' as a variant can also be inferred to be S where maxaf<.75, so change these to 'aS' as well: 
gen byte t=1 if character=="aS"
replace t=0 if t!=1
bysort drug variant: egen tt=max(t)
replace pass=2 if tt==1 & mask==0 & character!="aS" & variant!="missing"
replace character="aS" if tt==1 
drop s1 s2 S1 S2 t tt mut MUT

save "master_data_file_`c'_results_T1_`y'.dta", replace
 
*6* Identify which samples to explore Tier 2 sequences for, and re-run algorithm for those
* This means everything except those isolates with an R phenotype and an R or U character mutation as those could all be explanations for the phenotype.
********************************************************************************
********************************************************************************
gen byte E=1 if pheno=="R" & inlist(character,"R","U")
replace E=0 if E!=1
bysort sample_id drug: egen e=max(E)
gen byte explanation=1 if e==1
drop e E

* unmask tier 2 mutations
replace mask=0 if tier==2 & neutral!="true" & inlist(effect,"synonymous_variant","initiator_codon_variant","stop_retained_variant")==0
replace character="" if tier==2 & mask==0
replace mask=0 if variant=="missing"

* mask all mutation in isolates that already have an explanation, or where an explanation cannot be ruled out:
replace mask=1 if explanation==1

* Quantify the number of non-masked variants per sample and characterise them according to rules of the algorithm (1st pass)
gen byte MUT=1 if mask==0 & variant!=""
replace MUT=0 if MUT!=1

bysort sample_id drug: egen mut=sum(MUT)

gen byte SOLO1=1 if mut==1 & MUT==1 & het==0 & variant!="missing" & character!="R"
replace SOLO1=0 if character=="R"  & mut==1 & MUT==1 & het==0 & variant!="missing" & SOLO==0
replace SOLO=1 if SOLO1==1
drop SOLO1

gen byte res=1 if mut==1 & MUT==1 & pheno=="R" & het==0 & variant!="missing"
replace res=0 if res!=1
bysort drug variant: egen Res=max(res)

replace character="R" if Res==1 & character==""
drop res Res

* determine the character to be S if always S
bysort drug variant pheno het: gen n=_N
bysort drug variant het: gen m=_N
replace character="aS" if n==m & pheno=="S" & character=="" & mask==0 & het==0 & variant!="missing" 

* determine the character to be S if always S as a solo
gen byte s1=1 if MUT==1 & mut==1 & pheno=="S" & het==0 & variant!="missing"
replace s1=0 if s1!=1
gen byte s2=1 if MUT==1 & mut==1 & het==0 & variant!="missing"
replace s2=0 if s2!=1
bysort drug variant: egen S1=sum(s1)
bysort drug variant: egen S2=sum(s2)
compress
replace character="aS" if S1==S2 & S1!=0 & character=="" & pheno=="S" & mask==0 & het==0 & variant!="missing"

* anything identified as aS by either rule ('always S' or 'always S when SOLO') as a variant can also be inferred to be S where maxaf<.75, so change these to 'aS' as well: 
gen byte t=1 if character=="aS"
replace t=0 if t!=1
bysort drug variant: egen tt=max(t)
replace character="aS" if tt==1 
drop t tt S1 S2 s1 s2 mut MUT n m

* characterised what's left as 'U', including tier 2 variants that are masked in samples where there is a tier 1 'explanation'
replace chara="U" if chara=="" & het==0 & variant!="missing" 

* if a variant is U then it should also be a U as a het
gen byte t=1 if character=="U"
replace t=0 if t!=1
bysort drug variant: egen tt=max(t)
replace character="U" if tt==1 
drop t tt 

* lable variants that only ever appear as hets as 'U' unless already 'S'
bysort drug variant: gen n=_N
bysort drug variant het: gen m=_N
replace character="U" if het==1 & n==m & tier==2 & character!="S" & variant!="missing"
drop n m

drop mask explanation

save "master_data_file_`c'_results_T1andT2_`y'.dta", replace


*7* Run this again for synonymous mutation in tier 1. We'll call it tier 3. Only identify 'R' mutations among these. The others will be set back to 'S'.
********************************************************************************
********************************************************************************
use "master_data_file_`c'_results_T1andT2_`y'.dta", clear

gen byte E=1 if pheno=="R" & inlist(character,"R","U")
replace E=0 if E!=1
bysort sample_id drug: egen e=max(E)
gen byte explanation=1 if e==1
drop e E

* Mask the following mutations: all tier 2; all neutral mutations; aS mutations in tier 1. Unmask synonymous mutations in tier 1 that aren't neutral.
gen byte mask=1 if neutral=="true" | character=="aS"
replace mask=0 if mask!=1

* As we're focussing on tier 1 first, we mask tier 2 variants for now too
replace mask=1 if tier==2

* mask all mutation in isolates that already have an explanation, or where an explanation cannot be ruled out:
replace mask=1 if explanation==1

* Quantify the number of non-masked variants per sample and characterise them according to rules of the algorithm (1st pass)
gen byte MUT=1 if mask==0 & variant!=""
replace MUT=0 if MUT!=1

bysort sample_id drug: egen mut=sum(MUT)

gen byte res=1 if mut==1 & MUT==1 & pheno=="R" & het==0 & variant!="missing"
replace res=0 if res!=1
bysort drug variant: egen Res=max(res)
replace character="sR" if Res==1 & inlist(effect,"synonymous_variant","initiator_codon_variant","stop_retained_variant")
drop res Res

drop mask explanation mut MUT

save "master_data_file_`c'_results_T1andT2andT3_`y'.dta", replace


*8* Calculate the number of times each variant is seen in an isolate with an R phenotype over the total number of isolates in which each mutation is seen (R/tot).
********************************************************************************
********************************************************************************
use "master_data_file_`c'_results_T1andT2andT3_`y'.dta", clear

* For each variant we want the number of times it is seen in a sample with an R phenotype ('totR') and the total number of isolates in which it is seen ('tot').
bysort drug variant pheno het: gen TotR=_N
replace TotR=0 if pheno=="S" | het==1
bysort drug variant het: egen totR=max(TotR)
bysort drug variant het: gen tot=_N
replace tot=0 if het==1 | variant=="" 
replace totR=0 if variant==""
drop TotR

* For each variant we want the number of times it is seen as a SOLO in a sample with an R phenotype, and the total number of times it is seen as a SOLO. For R and aS mutations this is straightforward. For S non-synonymous neutral mutations this involves unmasking those S mutations in tier 1. For tier 2 mutations this involves unmasking tier two mutations in those samples without an 'explanation' (i.e. without an R, U or het). For synonymous mutations this involves unmasking synonymous mutations in tier 1. 

* Mask tier 2 mutations, neutral and synonymous variants to identify SOLOs among tier 1 mutations
gen byte mask=1 if neutral=="true" | inlist(effect,"synonymous_variant","initiator_codon_variant","stop_retained_variant")
replace mask=0 if mask!=1
replace mask=1 if tier==2

* Number of RSOLOs and TotSOLOs for 'R' and 'aS' variants:
bysort drug variant pheno SOLO: gen RSOLO=_N
replace RSOLO=0 if pheno=="S" | SOLO==0 
bysort drug variant: egen Rsolo=max(RSOLO)
drop RSOLO
bysort drug variant SOLO: gen TotSOLO=_N
replace TotSOLO=0 if SOLO==0 
bysort drug variant: egen totsolo=max(TotSOLO)
drop TotSOLO
replace Rsolo=0 if variant=="" | effect=="missing"
replace totsolo=0 if variant=="" | effect=="missing"

* Number of RSOLOs and TotSOLOs for non-synonymous neutral ('S') variants in tier 1 (i.e. neutral variants):
replace mask=0 if tier==1 & inlist(effect,"synonymous_variant","initiator_codon_variant","stop_retained_variant")==0 & neutral=="true"

gen n=1 if mask==0
replace n=0 if n!=1
bysort sample_id drug: egen sum_n=sum(n)
replace SOLO=1 if sum_n==1 & n==1 
replace SOLO=0 if SOLO!=1 
replace SOLO=0 if het==1 | effect=="missing"
drop n sum_n

bysort drug variant pheno SOLO: gen RSOLO=_N
replace RSOLO=0 if pheno=="S" | SOLO==0
bysort drug variant: egen Rsolo1=max(RSOLO)
drop RSOLO
bysort drug variant SOLO: gen TotSOLO=_N
replace TotSOLO=0 if SOLO==0
bysort drug variant: egen totsolo1=max(TotSOLO)
drop TotSOLO 
replace Rsolo=Rsolo1 if character=="S" & inlist(effect,"synonymous_variant","initiator_codon_variant","stop_retained_variant")==0 & neutral=="true" 
replace totsolo=totsolo1 if character=="S" & inlist(effect,"synonymous_variant","initiator_codon_variant","stop_retained_variant")==0 & neutral=="true" 
replace mask=1 if tier==1 & inlist(effect,"synonymous_variant","initiator_codon_variant","stop_retained_variant")==0 & neutral=="true"
drop Rsolo1 totsolo1

* Number of RSOLOs and TotSOLOs tier 2 variants:
gen byte E=1 if pheno=="R" & inlist(character,"R","U") & tier==1
replace E=0 if E!=1
bysort sample_id drug: egen e=max(E)
gen byte explanation=1 if e==1
drop e E mask

* mask all mutations in isolates that already have an explanation, or where an explanation cannot be ruled out:
gen byte mask=1 if neutral=="true" | inlist(effect,"synonymous_variant","initiator_codon_variant","stop_retained_variant")
replace mask=0 if mask!=1
replace mask=1 if explanation==1
replace mask=1 if tier==1 & character=="aS"

gen n=1 if mask==0
replace n=0 if n!=1
bysort sample_id drug: egen sum_n=sum(n)
replace SOLO=1 if sum_n==1 & n==1
replace SOLO=0 if SOLO!=1
replace SOLO=0 if het==1 | effect=="missing"

replace SOLO=0 if explanation==1

drop n sum_n

bysort drug variant pheno SOLO: gen RSOLO=_N
replace RSOLO=0 if pheno=="S" | SOLO==0
bysort drug variant: egen Rsolo1=max(RSOLO)
drop RSOLO
bysort drug variant SOLO: gen TotSOLO=_N
replace TotSOLO=0 if SOLO==0
bysort drug variant: egen totsolo1=max(TotSOLO)
drop TotSOLO
replace Rsolo=Rsolo1 if tier==2
replace totsolo=totsolo1 if tier==2
drop Rsolo1 totsolo1 mask

* Number of RSOLOs and TotSOLOs for synonymous ('S') variants in tier 1:
drop explanation
gen byte E=1 if pheno=="R" & inlist(character,"R","U")
replace E=0 if E!=1
bysort sample_id drug: egen e=max(E)
gen byte explanation=1 if e==1
drop e E

gen byte mask=1 if tier==2 | neutral=="true" | explanation==1 | character=="aS"
replace mask=0 if mask!=1

gen n=1 if mask==0
replace n=0 if n!=1
bysort sample_id drug: egen sum_n=sum(n)
replace SOLO=1 if sum_n==1 & n==1
replace SOLO=0 if SOLO!=1
replace SOLO=0 if het==1 | effect=="missing"
drop n sum_n

bysort drug variant pheno SOLO: gen RSOLO=_N
replace RSOLO=0 if pheno=="S" | SOLO==0
bysort drug variant: egen Rsolo1=max(RSOLO)
drop RSOLO
bysort drug variant SOLO: gen TotSOLO=_N
replace TotSOLO=0 if SOLO==0
bysort drug variant: egen totsolo1=max(TotSOLO)
drop TotSOLO
replace Rsolo=Rsolo1 if  inlist(effect,"synonymous_variant","initiator_codon_variant","stop_retained_variant") & tier==1
replace totsolo=totsolo1 if  inlist(effect,"synonymous_variant","initiator_codon_variant","stop_retained_variant") & tier==1
replace mask=1 if inlist(effect,"synonymous_variant","initiator_codon_variant","stop_retained_variant")
drop Rsolo1 totsolo1 mask SOLO

* generate the denominator for R and S for each drug (the total number of R and S isolates for each drug)
bysort sample_id drug: gen m=_n 
replace m=0 if m!=1
replace m=0 if pheno=="S"
bysort drug: egen R_denominator=sum(m)
drop m

bysort sample_id drug: gen m=_n 
replace m=0 if m!=1
replace m=0 if pheno=="R"
bysort drug: egen S_denominator=sum(m)
drop m

replace Rsolo=0 if variant=="missing"
replace totsolo=0 if variant=="missing"
replace tot=0 if variant=="missing"
replace totR=0 if variant=="missing"

save "master_data_file_`c'_results_T1andT2andT3_totals_`y'.dta", replace

*9* Do separate calculations for katG S315T and rpoB S450L as the QC rule above excluded samples with these mutations and S phenoypes, which has left them with an artificial 100% specificity
********************************************************************************
********************************************************************************

use "master_data_file_`c'_`y'.dta", clear

keep if inlist(drug,"Isoniazid","Rifampicin")

* Quantify the number of non-masked variants per sample and characterise them according to rules of the algorithm (1st pass)
gen byte MUT=1 if mask==0 & variant!=""
replace MUT=0 if MUT!=1

bysort sample_id drug: egen mut=sum(MUT)

gen byte res=1 if mut==1 & MUT==1 & pheno=="R" & het==0 & variant!="missing"
replace res=0 if res!=1
bysort drug variant: egen Res=max(res)
gen character="R" if Res==1 
drop res Res

* For each variant we want the number of times it is seen in a sample with an R phenotype ('totR') and the total number of isolates in which it is seen ('tot').
bysort drug variant pheno het: gen TotR=_N
replace TotR=0 if pheno=="S" | het==1 
bysort drug variant het: egen totR1=max(TotR)
bysort drug variant het: gen tot1=_N
replace tot1=0 if het==1 | variant=="" 
replace totR1=0 if variant==""
drop TotR

* For each variant we want the number of times it is seen as a SOLO in a sample with an R phenotype, and the total number of times it is seen as a SOLO. 
gen n=1 if mask==0
replace n=0 if n!=1
bysort sample_id drug: egen sum_n=sum(n)
gen SOLO=1 if sum_n==1 & n==1
replace SOLO=0 if SOLO!=1 
replace SOLO=0 if het==1 | effect=="missing"
drop n sum_n

* Number of RSOLOs and TotSOLOs for 'R' and 'aS' variants:
bysort drug variant pheno SOLO: gen RSOLO=_N
replace RSOLO=0 if pheno=="S" | SOLO==0 
bysort drug variant: egen Rsolo1=max(RSOLO)
drop RSOLO
bysort drug variant SOLO: gen TotSOLO=_N
replace TotSOLO=0 if SOLO==0
bysort drug variant: egen totsolo1=max(TotSOLO)
drop TotSOLO
replace Rsolo1=0 if variant=="" 
replace totsolo1=0 if variant=="" 

* generate the denominator for R and S for each drug (the total number of R and S isolates for each drug)
bysort sample_id drug: gen m=_n 
replace m=0 if m!=1
replace m=0 if pheno=="S"
bysort drug: egen R_denominator1=sum(m)
drop m

bysort sample_id drug: gen m=_n 
replace m=0 if m!=1
replace m=0 if pheno=="R"
bysort drug: egen S_denominator1=sum(m)
drop m


keep if inlist(variant,"katG_p.Ser315Thr","rpoB_p.Ser450Leu")
keep drug variant totR1 tot1 Rsolo1 totsolo1 R_denominator1 S_denominator1
bysort variant: gen n=_n
keep if n==1
drop n
save "katG_rpoB_`c'.dta", replace

merge 1:m drug variant using "master_data_file_`c'_results_T1andT2andT3_totals_`y'.dta"
replace totR=totR1 if _m==3
replace tot=tot1 if _m==3
replace Rsolo=Rsolo1 if _m==3
replace totsolo=totsolo1 if _m==3
replace R_denominator=R_denominator1 if _m==3
replace S_denominator=S_denominator1 if _m==3

drop totR1 tot1 Rsolo1 totsolo1 R_denominator1 S_denominator1 _m

save "master_data_file_`c'_results_T1andT2andT3_totals1_`y'.dta", replace

}
}


*10* Create pooled entry for frameshifts and premature stop codons to inform the expert rules
********************************************************************************
********************************************************************************

foreach y in complete noNICD {

foreach c in WHO ALL CC CCATU {

use "master_data_file_`c'_results_T1andT2andT3_totals1_`y'.dta", clear

preserve
gen byte lof=1 if inlist(effect,"stop_gained","start_lost","frameshift","feature_ablation")

* annotate individual variants as by different types of loss of function (lof):  
replace variant=gene+"_lof" if lof==1 & inlist(effect,"stop_gained","start_lost","frameshift","feature_ablation")

drop lof

keep if variant==gene+"_lof" & het==0

bysort sample_id drug gene: gen n=_n
keep if n==1
drop n	
replace het=0
replace neutral=""
replace maxaf=.
replace position=""
replace character=""
replace pass=.

save "pooled_lof_`c'_`y'.dta", replace
restore 

append using "pooled_lof_`c'_`y'.dta"

gen byte mask=1 if neutral=="true" | inlist(effect,"synonymous_variant","initiator_codon_variant","stop_retained_variant")
replace mask=0 if mask!=1
replace mask=1 if tier==2

gen byte lof=1 if inlist(effect,"stop_gained","start_lost","frameshift","feature_ablation")

replace mask=1 if variant!=gene+"_lof" & lof==1
replace character="" if inlist(character,"R","aS","U")

gen byte MUT=1 if variant!="" & mask==0
replace MUT=0 if MUT!=1
compress
bysort sample_id drug: egen mut=sum(MUT)

gen byte res=1 if mut==1 & MUT==1 & pheno=="R" & het==0 & variant!="missing"
replace res=0 if res!=1

bysort drug variant: egen Res=max(res)
compress
replace character="R" if Res==1 & tier==1

drop res Res

* determine the character to be S if always S
bysort drug variant pheno het: gen n=_N
bysort drug variant het: gen m=_N
replace charac="aS" if n==m & pheno=="S" & character!="S" & mask==0 & het==0 & variant != "missing"


* determine the character to be S if always S as a solo
gen byte s1=1 if MUT==1 & mut==1 & pheno=="S" & het==0 & variant!="missing"
replace s1=0 if s1!=1

gen byte s2=1 if MUT==1 & mut==1 & het==0 & variant!="missing"
replace s2=0 if s2!=1
bysort drug variant : egen S1=sum(s1)
bysort drug variant : egen S2=sum(s2)
compress
replace character="aS" if S1==S2 & S1!=0 & pheno=="S" & mask==0 & character!="S" 


* anything identified as 'aS' by either rule ('always S' or 'always S when SOLO') as a variant can also be inferred to be S where maxaf<.75, so change these to 'aS' as well: 
gen byte t=1 if character=="aS"
replace t=0 if t!=1
bysort drug variant: egen tt=max(t)
replace character="aS" if tt==1
drop t tt S1 S2 s1 s2 mut MUT n m

* characterised what's left as 'U'
replace character="U" if character=="" & tier==1 & mask==0 

* quantify Rtot/tot and Rsolo/totsolo for the pooled frameshifts and premature stop codons:
bysort drug variant pheno: gen TotR1=_N
replace TotR1=0 if pheno=="S" 
bysort drug variant : egen totR1=max(TotR)
bysort drug variant : gen tot1=_N
replace tot=tot1 
replace totR=totR1 

drop TotR1 tot1 totR1

* number of unmasked mutations as SOLO per sample / drug combination
gen n=1 if mask==0
replace n=0 if n!=1
bysort sample_id drug: egen sum_n=sum(n)
gen SOLO=1 if sum_n==1 & n==1
replace SOLO=0 if SOLO!=1
replace SOLO=0 if het==1
drop n sum_n

* Number of RSOLOs and TotSOLOs for 'R' and 'aS' variants:
bysort drug variant pheno SOLO: gen RSOLO=_N
replace RSOLO=0 if pheno=="S" | SOLO==0 
bysort drug variant: egen Rsolo1=max(RSOLO)
drop RSOLO
bysort drug variant SOLO: gen TotSOLO=_N
replace TotSOLO=0 if SOLO==0 
bysort drug variant: egen totsolo1=max(TotSOLO)

replace Rsolo=Rsolo1 
replace totsolo=totsolo1
drop TotSOLO Rsolo1 totsolo1 SOLO

* repeat for Tier 2
drop explanation
gen byte E=1 if pheno=="R" & inlist(character,"R","U") & tier==1
replace E=0 if E!=1
bysort sample_id drug: egen e=max(E)
gen byte explanation=1 if e==1
drop e E

* unmask tier 2 mutations
replace mask=0 if tier==2 & neutral!="true" & inlist(effect,"synonymous_variant","initiator_codon_variant","stop_retained_variant")==0

replace mask=1 if variant!=gene+"_lof" & lof==1

replace character="" if tier==2 & mask==0 & inlist(charact,"R","aS","U")
replace mask=1 if explanation==1 | character=="aS"

gen byte MUT=1 if variant!="" & mask==0
replace MUT=0 if MUT!=1
compress
bysort sample_id drug: egen mut=sum(MUT)

gen byte res=1 if mut==1 & MUT==1 & pheno=="R" & het==0 & variant!="missing"
replace res=0 if res!=1
bysort drug variant: egen Res=max(res)
compress
replace character="R" if Res==1 & tier==2 

drop res Res

* determine the character to be S if always S
bysort drug variant pheno het: gen n=_N
bysort drug variant het: gen m=_N
replace charac="aS" if n==m & pheno=="S" & character!="S" & mask == 0 & het == 0 & variant != "missing"


* determine the character to be S if always S as a solo
gen byte s1=1 if MUT==1 & mut==1 & pheno=="S"  & het==0 & variant!="missing"
replace s1=0 if s1!=1

gen byte s2=1 if MUT==1 & mut==1 & het==0 & variant!="missing"

replace s2=0 if s2!=1
bysort drug variant het: egen S1=sum(s1)
bysort drug variant het: egen S2=sum(s2)
compress
replace character="aS" if S1==S2 & S1!=0 & pheno=="S" & mask==0 & tier==2

* anything identified as 'aS' by either rule ('always S' or 'always S when SOLO') as a variant can also be inferred to be S where maxaf<.75, so change these to 'aS' as well: 
gen byte t=1 if character=="aS"
replace t=0 if t!=1
bysort drug variant: egen tt=max(t)
replace character="aS" if tt==1
drop t tt S1 S2 s1 s2 mut MUT n m

* characterised what's left as 'U'
replace character="U" if character=="" & tier==2


* quantify Rtot/tot and Rsolo/totsolo for the pooled frameshifts
bysort drug variant pheno : gen TotR1=_N
replace TotR=0 if pheno=="S" | het==1 
bysort drug variant : egen totR1=max(TotR)
bysort drug variant : gen tot1=_N
replace tot=tot1 if tier==2
replace totR=totR1 if tier==2

drop TotR1 tot1 totR1 

* number of unmasked mutations as SOLO per sample / drug combination
gen n=1 if mask==0
replace n=0 if n!=1
bysort sample_id drug: egen sum_n=sum(n)
gen SOLO=1 if sum_n==1 & n==1
replace SOLO=0 if SOLO!=1
drop n sum_n

* Number of RSOLOs and TotSOLOs for 'R' and 'aS' variants:
bysort drug variant pheno SOLO: gen RSOLO=_N
replace RSOLO=0 if pheno=="S" | SOLO==0 
bysort drug variant: egen Rsolo1=max(RSOLO)
drop RSOLO
bysort drug variant SOLO: gen TotSOLO=_N
replace TotSOLO=0 if SOLO==0 
bysort drug variant: egen totsolo1=max(TotSOLO)

replace Rsolo=Rsolo1 if  tier==2
replace totsolo=totsolo1 if tier==2
drop TotSOLO Rsolo1 totsolo1 mask SOLO explanation
compress

keep if variant==gene+"_lof"

save "master_data_file_`c'_pooled_lof_`y'.dta", replace

use  "master_data_file_`c'_results_T1andT2andT3_totals1_`y'.dta"

append using "master_data_file_`c'_pooled_lof_`y'.dta"

save "master_data_file_`c'_results_T1andT2andT3_totals2_`y'.dta", replace

}

use "master_data_file_ALL_results_T1andT2andT3_totals2_`y'.dta", clear

append using "master_data_file_WHO_results_T1andT2andT3_totals2_`y'.dta"
compress

save "master_data_file_results_T1andT2andT3_totals2_`y'.dta", replace

}


* 11. Compute the stats
********************************************************************************
********************************************************************************

set more off, permanently

foreach y in complete noNICD {
	
foreach a in results CC_results CCATU_results {

use "master_data_file_`a'_T1andT2andT3_totals2_`y'.dta", clear


// first prepare the data 
drop if inlist(variant,"missing","") | het==1 
keep drug variant tier tot totR Rsolo totsolo S_denominator R_denominator character pass neutral effect category position
bysort drug variant category: gen n=_n

sort drug variant category n
by drug variant category n: assert _N==1
by drug variant category: assert n==1 if _n==1  
keep if n==1
drop n

* this quantifies the 'k' for the Benjamini Hochberg (BH) correction:
* identify the literature mutations:
merge m:1 drug variant using "detailed_neutral_sets_`y'.dta"
drop if _m==2
gen v=1 
replace v=0 if  ((literature==1 | v1==1 ) & _m==3) | variant=="" | inlist(effect,"synonymous_variant","initiator_codon_variant","stop_retained_variant")
drop _m setA setB setD1 setD2 v1 literature setC n
bysort drug category_phenotype: egen k=sum(v)
drop v

expand 4
bysort drug variant category: gen m=_n
gen pheno=0 if m<3
replace pheno=1 if pheno!=0
gen mut=1 if inlist(m,1,3)
replace mut=0 if inlist(m,2,4)
bysort drug variant category: assert  _N==4
compress


preserve
* need to drop isolates with the mutation and another accompanying mutation from the denominator to facilitate the solo analysis. Therefore subtract from the denominator the difference between the number of times a mutation is seen with a phenotype and the number of times it is seen as a solo with that phenotype, as we'll leave those samples out. 
gen Rsolo_denominator=R_denominator-(totR-Rsolo)
* Ssolo_denominator becomes the number of S isolates, less the difference between mut_S and mutSOLO_S (i.e. less those that have the mut but not as a solo)
gen Ssolo_denominator=S_denominator-((tot-totR)-(totsolo-Rsolo))
gen n=Rsolo if pheno==1 & mut==1
replace n=Rsolo_denominator-Rsolo if pheno==1 & mut==0
replace n=totsolo-Rsolo if pheno==0 & mut==1
replace n=Ssolo_denominator-(totsolo-Rsolo) if pheno==0 & mut==0
gen mut_S=tot-totR if pheno==0 & mut==1
rename S_denominator Sall_denominator
drop R_denominator 
* drop variants that are never solo:
drop if totsolo==0
* drop those where there are no R isolates other than those containing this variant, and where those are never solo (so cannot be assessed):
drop if Rsolo_denominator==0
save "master_stats_solo_`a'_`y'.dta", replace
restore

gen n=totR if pheno==1 & mut==1
replace n=R_denominator-totR if pheno==1 & mut==0
replace n=tot-totR if pheno==0 & mut==1
replace n=S_denominator-(tot-totR) if pheno==0 & mut==0
gen mut_S=tot-totR
gen mut_soloS=totsolo-Rsolo

rename R_denominator Rall_denominator 
rename S_denominator Sall_denominator 

save "master_stats_all_`a'_`y'.dta", replace


** Perform Fischer's exact test on all drug-variant combinations, compute OR, sensitivity, specificity, and all CIs. 
set more off, permanently
foreach i in all solo  {

use "master_stats_`i'_`a'_`y'.dta", clear

sort drug variant category

* Create variables (mut_S, mut_R, nomut_S, nomut_R) and for SOLOs (mut_soloS, mut_soloR, nomut_S, nomut_R)

foreach m in 1 2 3 4 {
gen N`m'=n if m==`m'
bysort variant drug category: egen n`m'=max(N`m')
drop N`m'
}

rename totR r
rename n1 mut_`i'S
rename n2 nomut_S
rename n3 mut_`i'R
rename n4 nomut_R
	
* for the PPV we will need the total number of mut_S for both 'solos' and 'all', so generate a separate variable for that here
bysort drug variant category: egen mut_s=max(mut_S)
drop mut_S
gen allmut_conditional_`i'=mut_s+mut_`i'R 
drop if allmut_conditional_`i'==0 
gen allmut_`i'=mut_`i'S+mut_`i'R 
compress

* generate a group identifier by drug and variant
egen drugvar = group(drug variant category)

bysort drug variant category: assert  _N==4

sort drugvar m

preserve

keep pheno mut n drugvar R`i'_denominator S`i'_denominator  allmut_`i'  allmut_conditional_`i' 

gen ph_mu="sr" if pheno ==0 & mut==1
replace ph_mu="ss" if pheno ==0 & mut==0
replace ph_mu="rs" if pheno ==1 & mut==0
replace ph_mu="rr" if pheno ==1 & mut==1
drop pheno mut 
reshape wide n, i(drugvar ) j(ph_mu, string  )

bysort nrr nrs nsr nss : gen z=_n
tab z if z==1
keep if z==1
drop z drugvar 

gen drugvar=_n

reshape long n, i(drugvar ) j(q, string)

gen byte pheno=1 if substr(q,1,1)=="r"
replace pheno=0 if substr(q,1,1)=="s"

gen byte mut=1 if substr(q,2,1)=="r"
replace mut=0 if substr(q,2,1)=="s"

drop q

gen double or_`i'_p=.
compress

* this leaves behind an r(max), the total number of groups:
su drugvar, meanonly

* loop through each group from 1 to r(max) to perform a Fischer's exact test and replace the p-value with the output
quietly forvalues q = 1/`r(max)' {
	
	tabulate pheno mut [w=n] if drugvar==`q', exact
	replace or_`i'_p=r(p_exact) if drugvar==`q'
}

gen ph_mu="sr" if pheno ==0 & mut==1
replace ph_mu="ss" if pheno ==0 & mut==0
replace ph_mu="rs" if pheno ==1 & mut==0
replace ph_mu="rr" if pheno ==1 & mut==1
drop pheno mut 
reshape wide n, i(drugvar ) j(ph_mu, string  )

rename nrr mut_`i'R
rename nrs nomut_R
rename nsr mut_`i'S
rename nss nomut_S


* generate the odds ratios
gen double or_`i'=.
replace or_`i'= (mut_`i'R * nomut_S)/(mut_`i'S * nomut_R)

* and their exact (binomial) confidence intervals
gen or_`i'_ub=exp(ln(or_`i')+1.96*sqrt((1/mut_`i'R)+(1/mut_`i'S)+(1/nomut_R)+(1/nomut_S))) 
gen or_`i'_lb=exp(ln(or_`i')-1.96*sqrt((1/mut_`i'R)+(1/mut_`i'S)+(1/nomut_R)+(1/nomut_S))) 

drop drugvar


*generate sensitity and specifity, and their exact CIs.
	
gen sens_`i'=.
gen sens_`i'_lb=.
gen sens_`i'_ub=.

quietly forval q=1/`=_N' {
	cii proportions `=R`i'_denominator[`q']'  `=mut_`i'R[`q']' 
	replace sens_`i'=r(proportion) in `q'
	replace sens_`i'_ub=r(ub) in `q'
	replace sens_`i'_lb=r(lb) in `q'
	compress
}

gen spec_`i'=.
gen spec_`i'_lb=.
gen spec_`i'_ub=.

quietly forval q=1/`=_N' {
	cii proportions `=S`i'_denominator[`q']'  `=nomut_S[`q']' 
	replace spec_`i'=r(proportion) in `q'
	replace spec_`i'_ub=r(ub) in `q'
	replace spec_`i'_lb=r(lb) in `q'
	compress
}


gen ppv_`i'=.
gen ppv_`i'_lb=.
gen ppv_`i'_ub=.

quietly forval q=1/`=_N' {
	cii proportions `=allmut_`i'[`q']' `=mut_`i'R[`q']' 
	replace ppv_`i'=r(proportion) in `q'
	replace ppv_`i'_ub=r(ub) in `q'
	replace ppv_`i'_lb=r(lb) in `q'
	compress
}
	
* generate PPV and CIs for 'PPV | SOLO'
gen ppv_conditional_`i'=.
gen ppv_conditional_`i'_lb=.
gen ppv_conditional_`i'_ub=.

quietly forval q=1/`=_N' {
	cii proportions `=allmut_conditional_`i'[`q']' `=mut_`i'R[`q']' 
	replace ppv_conditional_`i'=r(proportion) in `q'
	replace ppv_conditional_`i'_ub=r(ub) in `q'
	replace ppv_conditional_`i'_lb=r(lb) in `q'
	compress	
}

compress

save "stats_`i'_`a'_`y'.dta", replace

restore

* reduce to one line per drug / mutation / category
bysort drug variant category: gen byte nn=_n
keep if nn==1
drop nn mut pheno m 

merge m:1 mut_`i'R nomut_R mut_`i'S nomut_S R`i'_denominator S`i'_denominator  allmut_`i'  allmut_conditional_`i' using "stats_`i'_`a'_`y'.dta"
drop _m allmut_*
compress


merge m:1 drug variant using "detailed_neutral_sets_`y'.dta"
drop if _m==2
gen z=1 if  ((literature==1 | v1==1 ) & _m==3) | variant=="" | inlist(effect,"synonymous_variant","initiator_codon_variant","stop_retained_variant")
drop _m setA setB setD1 setD2 v1 literature setC n
preserve
drop if z!=1
drop z
save "neutrals_not_to_be_counted_`i'_`a'_`y'.dta", replace
restore

drop if z==1
drop z

gsort drug category or_`i'_p - or_`i'
bysort drug category : gen i=_n
compress
gen or_`i'_p_bh=(i/k)*0.05
gen bh_thr_`i'=or_`i'_p if or_`i'_p_bh>or_`i'_p

replace bh_thr_`i' =or_`i'_p_bh if i==1 & bh_thr_`i' ==.
compress
bysort drug category: egen or_`i'_bh_threshold=max(bh_thr_`i')

gen or_`i'_bh_significant=1 if or_`i'_p<=or_`i'_bh_threshold 
replace or_`i'_bh_significant=0 if or_`i'_bh_significant!=1 

order drug variant category

append using "neutrals_not_to_be_counted_`i'_`a'_`y'.dta"

compress

save "master_stats_`i'_`a'2_`y'.dta", replace	
	
}
}


use "master_stats_solo_results2_`y'.dta", clear
drop effect neutral tier character pass r tot Rsolo totsolo drugvar nomut_S nomut_R mut_s k 

merge 1:1 drug variant category using "master_stats_all_results2_`y'.dta"
drop _m
compress

save "master_stats_primary_output_`y'.dta", replace

use "master_stats_solo_CC_results2_`y'.dta", clear
drop effect neutral tier character pass r tot Rsolo totsolo drugvar nomut_S nomut_R mut_s k 

merge 1:1 drug variant category using "master_stats_all_CC_results2_`y'.dta"
drop _m
compress

save "master_stats_CC_`y'.dta", replace

use "master_stats_solo_CCATU_results2_`y'.dta", clear
drop effect neutral tier character pass r tot Rsolo totsolo drugvar nomut_S nomut_R mut_s k 

merge 1:1 drug variant category using "master_stats_all_CCATU_results2_`y'.dta"
drop _m
compress

save "master_stats_CCATU_`y'.dta", replace

}



* 12. Implement the grading rules
********************************************************************************
********************************************************************************


* prep the genotypes that don't have an accompanying phenotype for merging in below: 
import delimited "run-1682403096096-part-r-00000.csv", varnames(1) clear 
gen variant=resolved_symbol+"_"+variant_category
rename drug_name drug
drop if maxaf<.75
keep sample_id variant
bysort sample_id variant: gen n=_n
keep if n==1
drop n
save "orphan_genotypes.dta", replace

* generate counts of all variants:
use "master_data_file_pheno_complete.dta", clear
drop if maxaf<.75
keep sample_id variant 
bysort sample_id variant : gen n=_n
keep if n==1
drop n
append using "orphan_genotypes.dta"
bysort sample_id variant : gen n=_n
keep if n==1
keep sample_id variant 
bysort variant : gen variant_count=_N
bysort variant : gen n=_n
keep if n==1
keep variant variant_count 
save "total_variant_counts_including_orphans.dta", replace

foreach y in complete noNICD {

foreach a in primary_output CC CCATU {

use "master_stats_`a'_`y'.dta", clear

keep drug tier variant pass mut_soloR mut_soloS mut_allS nomut_S mut_allR nomut_R ppv_all ppv_all_lb ppv_all_ub ppv_solo ppv_solo_lb ppv_solo_ub ppv_conditional_solo ppv_conditional_solo_lb ppv_conditional_solo_ub  or_all or_all_lb or_all_ub or_solo or_solo_lb or_solo_ub or_solo_bh_signi neutral category effect position   

rename pass algorithm_pass
replace mut_soloR=0 if mut_soloR==.
replace mut_soloS=0 if mut_soloS==.
gen Present_SOLO_SR=mut_soloR+mut_soloS
rename mut_soloR Present_SOLO_R
rename mut_allS Present_S
rename nomut_S Absent_S
rename mut_allR Present_R
rename nomut_R Absent_R
rename ppv_all PPV
rename ppv_all_lb PPV_lb
rename ppv_all_ub PPV_ub
rename ppv_solo PPV_SOLO
rename ppv_solo_lb PPV_SOLO_lb
rename ppv_solo_ub PPV_SOLO_ub
rename ppv_conditional_solo PPV_conditional_SOLO
rename ppv_conditional_solo_lb PPV_conditional_SOLO_lb
rename ppv_conditional_solo_ub PPV_conditional_SOLO_ub
rename or_all OR
rename or_all_lb OR_lb
rename or_all_ub OR_ub
rename or_solo OR_SOLO
rename or_solo_lb OR_SOLO_lb
rename or_solo_ub OR_SOLO_ub
rename or_solo_bh_signi OR_SOLO_FE_sig
rename neutral Neutral_masked
rename category Datasets

* correct the labelling of effects:
replace effect="frameshift" if strpos(variant,"_lof")

split variant,p("_")
rename variant1 gene
drop variant2 variant3

*INITIAL GRADING:
gen Initial_Confidence_Grading="3) Uncertain significance"

* Define NotAwR by WHO dataset only (neutral by the neutral algorithm)
replace Initial_Confidence_Grading="5) Not assoc w R" if Neutral_masked=="true"

* assoc w R by thresholds
replace Initial_Confidence_Grading="1) Assoc w R" if Present_SOLO_SR>=5 & PPV_conditional_SOLO_lb>=.25 & OR_SOLO>1 & OR_SOLO_FE_sig==1 & Present_SOLO_SR!=. & PPV_conditional_SOLO_lb!=. 

* For the uncertain ones, change to neutral – interim by thresholds (relaxed thresholds pncA only)
replace Initial_Confidence_Grading="4) Not assoc w R - Interim" if PPV_SOLO<.4 & PPV_SOLO_ub<.75 & strpos(variant,"pncA") & Initial_Confidence_Grading=="3) Uncertain significance" & Datasets=="WHO"

* For the uncertain ones, assoc w R – interim by thresholds (relaxed thresholds pncA only)
replace Initial_Confidence_Grading="2) Assoc w R - Interim" if Initial_Confidence_Grading=="3) Uncertain significance" & Present_SOLO_R>=2 & PPV>=.5 & strpos(variant,"pncA")


* FINAL GRADING
gen Final_Confidence_Grading=Initial_Confidence_Grading 

merge m:1 drug variant using "v1_grades.dta"
drop if _m==2
drop _m
replace final_grading_v1="" if inlist(effect,"inframe_deletion","inframe_insertion") & final_grading_v1=="2) Assoc w R - Interim"
gen Additionalgradingcriteria=""

*Downgrade to NotAwR-interim variants classified as neutral (NotAwR) by SetC only (i.e. setD1-2 = 0) but within that only by ver1 guidance or literature guidance (specify which one applies):
merge m:1 drug variant using "detailed_neutral_sets_`y'.dta"
drop if _m==2
drop _m n

replace Final_Confidence_Grading="4) Not assoc w R - Interim" if Initial_Confidence_Grading=="5) Not assoc w R" & setC==1 & setD1!=1 & setD2!=1 & setA!=1 & setB!=1 & literature!=1 & v1==1
replace Additionalgradingcriteria="Neutrals algorithm version 1 only" if Initial_Confidence_Grading=="5) Not assoc w R" & setC==1 & setD1!=1 & setD2!=1 & setA!=1 & setB!=1 & literature!=1 & v1==1

replace Final_Confidence_Grading="4) Not assoc w R - Interim" if Initial_Confidence_Grading=="5) Not assoc w R" & setC==1 & setD1!=1 & setD2!=1 & setA!=1 & setB!=1 & literature==1 & v1!=1
replace Additionalgradingcriteria="Neutrals algorithm literature only" if Initial_Confidence_Grading=="5) Not assoc w R" & setC==1 & setD1!=1 & setD2!=1 & setA!=1 & setB!=1 & literature==1 & v1!=1

*Downgrade to NotAwR-interim all silent variants that are currently Uncertain Significance
replace Final_Confidence_Grading="4) Not assoc w R - Interim" if Final_Confidence_Grading=="3) Uncertain significance" & inlist(effect,"synonymous_variant","initiator_codon_variant","stop_retained_variant")
replace Additionalgradingcriteria="Silent mutation" if Initial_Confidence_Grading=="3) Uncertain significance" & inlist(effect,"synonymous_variant","initiator_codon_variant","stop_retained_variant") & Final_Confidence_Grading=="4) Not assoc w R - Interim"

* Downgrade to AwR-interim variants classified as AwR by ALL dataset only (except for rpoB borderline mutations)
gen z=1 if Final_Confidence_Grading=="1) Assoc w R" & Datasets=="ALL"
replace z=1 if Final_Confidence_Grading=="1) Assoc w R" & Datasets=="WHO"
replace z=0 if z!=1
bysort drug variant: egen zz=sum(z)
gen ALLonly=1 if zz==1 & z==1 & Datasets=="ALL"
replace Final_Confidence_Grading="2) Assoc w R - Interim" if Initial_Confidence_Grading=="1) Assoc w R" & ALLonly==1 & inlist(variant,"rpoB_p.Leu430Pro","rpoB_p.Asp435Tyr","rpoB_p.His445Leu","rpoB_p.His445Asn","rpoB_p.His445Ser","rpoB_p.Leu452Pro","rpoB_p.Ile491Phe")==0
replace Additionalgradingcriteria="ALL dataset only" if Initial_Confidence_Grading=="1) Assoc w R" & ALLonly==1 & inlist(variant,"rpoB_p.Leu430Pro","rpoB_p.Asp435Tyr","rpoB_p.His445Leu","rpoB_p.His445Asn","rpoB_p.His445Ser","rpoB_p.Leu452Pro","rpoB_p.Ile491Phe")==0
drop z zz ALLonly

* Downgrade to AwR-interim PZA_pncA variants classified as AwR by ALL dataset but classified as AwR-interim by WHO dataset
gen z=1 if Final_Confidence_Grading=="1) Assoc w R" & Datasets=="ALL"
replace z=1 if Final_Confidence_Grading=="2) Assoc w R - Interim" & Datasets=="WHO"
replace z=0 if z!=1
bysort drug variant: egen zz=sum(z)
replace Final_Confidence_Grading="2) Assoc w R - Interim" if zz==2 & strpos(variant,"pncA") & effect!="upstream_gene_variant"
replace Additionalgradingcriteria="Downgraded due to WHO dataset" if zz==2 & strpos(variant,"pncA") 
drop z zz

* rpoB borderline mutations move grade up to "1) Assoc w R"
replace Final_Confidence_Grading="1) Assoc w R" if inlist(variant,"rpoB_p.Leu430Pro","rpoB_p.Asp435Tyr","rpoB_p.His445Leu","rpoB_p.His445Asn","rpoB_p.His445Ser","rpoB_p.Leu452Pro","rpoB_p.Ile491Phe") & inlist(Initial_Confidence_Grading,"2) Assoc w R - Interim","3) Uncertain significance")
replace Additionalgradingcriteria="Borderline" if inlist(variant,"rpoB_p.Leu430Pro","rpoB_p.Asp435Tyr","rpoB_p.His445Leu","rpoB_p.His445Asn","rpoB_p.His445Ser","rpoB_p.Leu452Pro","rpoB_p.Ile491Phe") & inlist(Initial_Confidence_Grading,"1) Assoc w R","2) Assoc w R - Interim","3) Uncertain significance")

* Upgrade to AwR-interim any non-silent variant in RIF_rpoB_RRDR. Identify RRDR and implement expert rule around that 
gen aa = regexs(0) if(regexm(variant, "[0-9][0-9][0-9]"))
destring aa, replace
gen rrdr=1 if strpos(variant,"rpoB_p.") & aa>=426 & aa<=452 & effect!="synonymous_variant"
drop aa
* and identify those indels where the first aa in the nomenclature falls outside of the rrdr but the second one falls within it (e.g. rpoB_p.Phe425_Gly426del)
gen aa = regexs(1) if(regexm(variant,"([0-9][0-9][0-9])[a-zA-Z]*$"))
destring aa, replace
replace rrdr=1 if strpos(variant,"rpoB_p.") & aa>=426 & aa<=452
replace Final_Confidence_Grading="2) Assoc w R - Interim" if Initial_Confidence_Grading=="3) Uncertain significance" & rrdr==1 & Additionalgradingcriteria!="Borderline"
replace Additionalgradingcriteria="RRDR" if Initial_Confidence_Grading=="3) Uncertain significance" & rrdr==1 & Additionalgradingcriteria!="Borderline"
drop rrdr aa

* Upgrade to AwR-interim based on "allelic exchange" evidence
merge m:1 drug variant using "allelic_exchanges.dta"

replace Final_Confidence_Grading="2) Assoc w R - Interim" if Initial_Confidence_Grading=="3) Uncertain significance" & _m==3 
replace Final_Confidence_Grading="2) Assoc w R - Interim" if _m==2 
replace Initial_Confidence_Grading="NA" if _m==2 


replace Additionalgradingcriteria="Evidence from allelic exchange experiments" if _m==3 & Initial_Confidence_Grading=="3) Uncertain significance" 
replace Additionalgradingcriteria="Evidence from allelic exchange experiments" if _m==2
drop _m

* if pooled LoF = AwR, then any LoF "Uncertain" mutation in the same drug_gene should be upgraded to AwR-interim 
gen z=1 if strpos(variant,"_lof") & inlist(Final_Confidence_Grading,"1) Assoc w R","2) Assoc w R - Interim") 
replace z=0 if z!=1
bysort drug gene: egen zz=max(z)
replace Final_Confidence_Grading="2) Assoc w R - Interim" if Initial_Confidence_Grading=="3) Uncertain significance" & inlist(effect,"feature_ablation","frameshift","start_lost","stop_gained") & strpos(variant,"_lof")==0 & zz==1
replace Additionalgradingcriteria="Indel frameshift or premature stop codon (LoF)" if Initial_Confidence_Grading=="3) Uncertain significance" & inlist(effect,"feature_ablation","frameshift","start_lost","stop_gained") & strpos(variant,"_lof")==0 & zz==1
gen comment="Can only confer resistance if genetically linked to a functional MmpL5" if (Additionalgradingcriteria=="Indel frameshift or premature stop codon (LoF)" | strpos(variant,"_lof")) & strpos(variant,"Rv0678")
drop z zz

* Upgrade to AwR-interim based on "recognized as DR marker" WHO-endorsed assay
merge m:1 drug variant using "assay_mutations.dta"

replace Final_Confidence_Grading="2) Assoc w R - Interim" if Initial_Confidence_Grading=="3) Uncertain significance" & _m==3 & Additionalgradingcriteria!="Borderline"
replace Final_Confidence_Grading="2) Assoc w R - Interim" if _m==2 

replace Additionalgradingcriteria="Recognized as DR marker through WHO-endorsed assay" if _m==3 & Initial_Confidence_Grading=="3) Uncertain significance" & inlist(Additionalgradingcriteria,"Borderline","RRDR","Indel frameshift or premature stop codon (LoF)")==0
replace Additionalgradingcriteria="Recognized as DR marker through WHO-endorsed assay" if _m==2 & inlist(Additionalgradingcriteria,"Borderline","RRDR","Indel frameshift or premature stop codon (LoF)")==0
drop _m

* Upgrade to AwR-interim based on Literature evidence (PMID 32571824) (applies to Pyrazinamide_pncA_p.Ile31Thr only)
replace Final_Confidence_Grading="2) Assoc w R - Interim" if  variant=="pncA_p.Ile31Thr" & Initial_Confidence_Grading=="3) Uncertain significance"
replace Additionalgradingcriteria="Literature evidence (PMID 32571824)" if  variant=="pncA_p.Ile31Thr" & Initial_Confidence_Grading=="3) Uncertain significance"

replace Final_Confidence_Grading="6) Manual check" if Initial_Confidence_Grading=="3) Uncertain significance" & inlist(final_grading_v1,"1) Assoc w R","2) Assoc w R - Interim") & inlist(Additionalgradingcriteria,"Borderline","RRDR","Recognized as DR marker through WHO-endorsed assay","Indel frameshift or premature stop codon (LoF)")==0 & strpos(Additionalgradingcriteria,"Literature evidence")==0
replace Additionalgradingcriteria="Evidence from catalogue version 1" if Initial_Confidence_Grading=="3) Uncertain significance" & inlist(final_grading_v1,"1) Assoc w R","2) Assoc w R - Interim") & inlist(Additionalgradingcriteria,"Borderline","RRDR","Indel frameshift or premature stop codon (LoF)","Recognized as DR marker through WHO-endorsed assay")==0 & strpos(Additionalgradingcriteria,"Literature evidence")==0

* Upgrade to AwR-interim based on FQ cross-resistance (LEV-MFX on gyrA and gyrB AwR and AwR-interim, bilateral)
gen fqm=1 if drug=="Moxifloxacin" & inlist(Final_Confidence_Grading,"1) Assoc w R","2) Assoc w R - Interim")
replace fqm=0 if fqm!=1
bysort variant: egen fqmox=max(fqm)
replace Final_Confidence_Grading="2) Assoc w R - Interim" if drug=="Levofloxacin" & fqmox==1 & Initial_Confidence_Grading=="3) Uncertain significance"
replace Additionalgradingcriteria="FQ cross-resistance" if drug=="Levofloxacin" & fqmox==1 & Initial_Confidence_Grading=="3) Uncertain significance"

gen fql=1 if drug=="Levofloxacin" & inlist(Final_Confidence_Grading,"1) Assoc w R","2) Assoc w R - Interim")
replace fql=0 if fql!=1
bysort variant: egen fqlev=max(fql)
replace Final_Confidence_Grading="2) Assoc w R - Interim" if drug=="Moxifloxacin" & fqlev==1 & Initial_Confidence_Grading=="3) Uncertain significance"
replace Additionalgradingcriteria="FQ cross-resistance" if drug=="Moxifloxacin" & fqlev==1 & Initial_Confidence_Grading=="3) Uncertain significance" & Additionalgradingcriteria!="WHO-endorsed assay"

drop fql fqm fqlev fqmox

* Upgrade to AwR-interim based on INH-ETH cross-resistance (INH-ETH on inhA AwR and AwR-interim, bilateral)
gen inha=1 if drug=="Isoniazid" & inlist(Final_Confidence_Grading,"1) Assoc w R","2) Assoc w R - Interim") & inlist(gene,"inhA","fabG1")
replace inha=0 if inha!=1
bysort variant: egen inhamax=max(inha)
replace Final_Confidence_Grading="2) Assoc w R - Interim" if drug=="Ethionamide" & inhamax==1 & Initial_Confidence_Grading=="3) Uncertain significance"
replace Additionalgradingcriteria="INH/ETH cross-resistance" if drug=="Ethionamide" & inhamax==1 & Initial_Confidence_Grading=="3) Uncertain significance" & Additionalgradingcriteria!="Recognized as DR marker through WHO-endorsed assay"

drop inha inhamax

gen inha=1 if drug=="Ethionamide" & inlist(Final_Confidence_Grading,"1) Assoc w R","2) Assoc w R - Interim") & inlist(gene,"inhA","fabG1")
replace inha=0 if inha!=1
bysort variant: egen inhamax=max(inha)
replace Final_Confidence_Grading="2) Assoc w R - Interim" if drug=="Isoniazid" & inhamax==1 & Initial_Confidence_Grading=="3) Uncertain significance"
replace Additionalgradingcriteria="INH/ETH cross-resistance" if drug=="Isoniazid" & inhamax==1 & Initial_Confidence_Grading=="3) Uncertain significance" & Additionalgradingcriteria!="Recognized as DR marker through WHO-endorsed assay"

drop inha inhamax

replace comment="Low-level resistance (multiple, genetically linked low-level resistance mutations are additive and confer high-level resistance)" if inlist(variant,"inhA_c.-770T>G","inhA_c.-779G>T") & drug=="Isoniazid"

* Upgrade to AwR-interim based on BDQ-CFZ cross-resistance (BDQ-CFZ on Rv0678 and pepQ AwR and AwR-interim, BDQ over CFZ only)
gen bdq=1 if drug=="Bedaquiline" & inlist(Final_Confidence_Grading,"1) Assoc w R","2) Assoc w R - Interim") & inlist(gene,"Rv0678","pepQ")
replace bdq=0 if bdq!=1
bysort variant: egen bdqmax=max(bdq)
replace Final_Confidence_Grading="2) Assoc w R - Interim" if drug=="Clofazimine" & bdqmax==1 & Initial_Confidence_Grading=="3) Uncertain significance"
replace Additionalgradingcriteria="BDQ/CFZ cross-resistance" if drug=="Clofazimine" & bdqmax==1 & Initial_Confidence_Grading=="3) Uncertain significance" & Additionalgradingcriteria!="Indel frameshift or premature stop codon (LoF)"

drop bdq bdqmax gene

* Downgrade to AwRi BDQ mutations driven by dataset that only included R isolates, potentially inflating PPV:
replace Final_Confidence_Grading="2) Assoc w R - Interim" if drug=="Bedaquiline" & inlist(variant,"Rv0678_p.Asn98fs","Rv0678_p.Cys46Arg","Rv0678_p.Cys46fs","Rv0678_p.Gln51fs","Rv0678_p.Ile67Ser","Rv0678_p.Leu142fs","Rv0678_p.Met146Thr","Rv0678_p.Pro48fs")
replace Additionalgradingcriteria="Potentially inflated PPV" if drug=="Bedaquiline" & inlist(variant,"Rv0678_p.Asn98fs","Rv0678_p.Cys46Arg","Rv0678_p.Cys46fs","Rv0678_p.Gln51fs","Rv0678_p.Ile67Ser","Rv0678_p.Leu142fs","Rv0678_p.Met146Thr","Rv0678_p.Pro48fs")
replace comment="Can only confer resistance if genetically linked to a functional MmpL5; Includes data from one site that only submitted resistant strains, which may have inflated the PPV" if drug=="Bedaquiline" & inlist(variant,"Rv0678_p.Asn98fs","Rv0678_p.Cys46Arg","Rv0678_p.Cys46fs","Rv0678_p.Gln51fs","Rv0678_p.Ile67Ser","Rv0678_p.Leu142fs","Rv0678_p.Met146Thr","Rv0678_p.Pro48fs")

* Specify PMIDs that lead to a downgrade:
replace Additionalgradingcriteria="Literature evidence (PMID 28031270; 34503982)" if drug=="Bedaquiline" & setA==0 & setB==0 & setC==1
replace Additionalgradingcriteria="Literature evidence (PMID 28137812)" if inlist(drug,"Levofloxacin","Moxifloxacin") & setA==0 & setB==0 & setC==1
replace Additionalgradingcriteria="Literature evidence (PMID 32143680)" if inlist(drug,"Ethambutol","Isoniazid","Pyrazinamide","Rifampicin") & setA==0 & setB==0 & setC==1
replace Additionalgradingcriteria="Previous WHO guidance" if setA==0 & setB==0 & setC==0 & setD1==0 & setD2==0 & v1==1 & literature==0
replace Final_Confidence_Grading="4) Not assoc w R - Interim" if Initial_Confidence_Grading=="5) Not assoc w R" & setA==0 & setB==0 & setC==0 & setD1==0 & setD2==0 & v1==1 & literature==0

* Specify PMIDs that lead to an upgrade:
replace Final_Confidence_Grading="4) Not assoc w R - Interim" if inlist(variant,"mmpL5_p.Ile948Val","Rv1979c_c.-129A>G") & drug=="Bedaquiline"
replace Additionalgradingcriteria="Literature evidence (PMID 28031270; 34503982)" if inlist(variant,"mmpL5_p.Ile948Val","Rv1979c_c.-129A>G") & drug=="Bedaquiline"

replace Final_Confidence_Grading="4) Not assoc w R - Interim" if inlist(variant,"embA_p.Ala813Gly","embA_p.Pro639Ser","embB_p.Gly156Cys","embB_p.Val668Ile","embC_c.-20A>C") & drug=="Ethambutol"
replace Additionalgradingcriteria="Literature evidence (PMID 32143680)" if inlist(variant,"embA_p.Ala813Gly","embA_p.Pro639Ser","embB_p.Gly156Cys","embB_p.Val668Ile","embC_c.-20A>C") & drug=="Ethambutol"

replace Final_Confidence_Grading="4) Not assoc w R - Interim" if inlist(variant,"katG_c.-354C>T","mshA_p.Asp218Ala","ndh_c.-70G>T","ndh_p.Arg268His","ndh_p.Gly313Arg") & drug=="Isoniazid"
replace Additionalgradingcriteria="Literature evidence (PMID 32143680)" if inlist(variant,"katG_c.-354C>T","mshA_p.Asp218Ala","ndh_c.-70G>T","ndh_p.Arg268His","ndh_p.Gly313Arg") & drug=="Isoniazid"

* mark mutations that are AwR/AwRi by ALL, and notAWR/notAwRi by WHO as 'for manual grading'
gen manual=1 if inlist(Initial_Confidence_Grading,"1) Assoc w R","2) Assoc w R - Interim") & Datasets=="ALL"
replace manual=1 if inlist(Initial_Confidence_Grading,"4) Not assoc w R - Interim","5) Not assoc w R") & Datasets=="WHO"
replace manual=0 if manual!=1
bysort drug variant: egen manual2=sum(manual)
replace Final_Confidence_Grading="6) Manual check" if manual2==2 
drop manual*

drop Additionalgradingcriteria_v1 final_grading_v1 

replace Additionalgradingcriteria=Additionalgradingcriteria+"Algorithm pass 2" if algorithm_pass==2 & Additionalgradingcriteria=="" & inlist(Final_Confidence_Grading,"1) Assoc w R","5) Not assoc w R")
replace Additionalgradingcriteria=Additionalgradingcriteria+"; Algorithm pass 2" if algorithm_pass==2 & Additionalgradingcriteria!="" & inlist(Final_Confidence_Grading,"1) Assoc w R","5) Not assoc w R")
replace Final_Confidence_Grading="2) Assoc w R - Interim" if Final_Confidence_Grading=="1) Assoc w R"  & algorithm_pass==2


replace comment="Consider that this variant might indicate that the isolate is M. canettii, which is intrinsically PZA resistance" if variant=="pncA_c.138A>G"

save "master_stats_`a'_i_`y'.dta", replace

}

use "master_stats_primary_output_i_`y'.dta", clear

replace Datasets ="na" if Datasets ==""

reshape wide mut_soloS Present_SOLO_R OR_SOLO OR_SOLO_ub OR_SOLO_lb PPV_SOLO PPV_SOLO_lb PPV_SOLO_ub PPV_conditional_SOLO PPV_conditional_SOLO_lb PPV_conditional_SOLO_ub OR_SOLO_FE_sig Present_S Absent_S Present_R  Absent_R OR OR_lb OR_ub PPV PPV_lb PPV_ub Present_SOLO_SR Initial_Confidence_Grading Additionalgradingcriteria Final_Confidence_Grading effect algorithm_pass comment , i(drug variant) j(Datasets , string)

drop mut_soloSna Present_SOLO_Rna OR_SOLOna OR_SOLO_ubna OR_SOLO_lbna PPV_SOLOna PPV_SOLO_lbna PPV_SOLO_ubna PPV_conditional_SOLOna PPV_conditional_SOLO_lbna PPV_conditional_SOLO_ubna OR_SOLO_FE_signa effectna algorithm_passna Present_Sna Absent_Sna Present_Rna Absent_Rna ORna OR_ubna OR_lbna PPVna PPV_lbna PPV_ubna Present_SOLO_SRna commentna

* define confidence grading for the variants added on basis of allelic exchange
gen Final_Confidence_Grading=Final_Confidence_Gradingna if Initial_Confidence_Gradingna=="NA" 
gen Additionalgradingcriteria=Additionalgradingcriteriana if Additionalgradingcriteriana=="Evidence from allelic exchange experiments"
gen Initial_Confidence_Grading=Initial_Confidence_Gradingna if Initial_Confidence_Gradingna=="NA"
replace Initial_Confidence_GradingWHO="NA" if Additionalgradingcriteriana=="Evidence from allelic exchange experiments"
replace Initial_Confidence_GradingALL="NA" if Additionalgradingcriteriana=="Evidence from allelic exchange experiments"
gen evidence="NA" if Additionalgradingcriteriana=="Evidence from allelic exchange experiments"
drop Initial_Confidence_Gradingna Additionalgradingcriteriana Final_Confidence_Gradingna		

* define Initial_Confidence_GradingWHO as Uncertain significance if it's blank
replace Initial_Confidence_GradingWHO="3) Uncertain significance" if Initial_Confidence_GradingWHO==""

* identify which results (WHO or ALL) are represented as final result, according to the following rules:

* rule 1
replace evidence="ALL+WHO" if Initial_Confidence_GradingWHO==Initial_Confidence_GradingALL & Initial_Confidence_GradingWHO!="NA"
replace Initial_Confidence_Grading=Initial_Confidence_GradingWHO if Initial_Confidence_GradingWHO==Initial_Confidence_GradingALL & Initial_Confidence_GradingWHO!="NA"
replace Additionalgradingcriteria=AdditionalgradingcriteriaALL if Initial_Confidence_GradingWHO==Initial_Confidence_GradingALL & Initial_Confidence_GradingWHO!="NA"
replace Final_Confidence_Grading=Final_Confidence_GradingALL if Initial_Confidence_GradingWHO==Initial_Confidence_GradingALL & Initial_Confidence_GradingWHO!="NA"

* rule 2
replace evidence="ALL" if Initial_Confidence_GradingWHO=="3) Uncertain significance" & Initial_Confidence_GradingALL!="3) Uncertain significance" & evidence!="NA"
replace Initial_Confidence_Grading=Initial_Confidence_GradingALL if Initial_Confidence_GradingWHO=="3) Uncertain significance" & Initial_Confidence_GradingALL!="3) Uncertain significance" & evidence!="NA"
replace Additionalgradingcriteria=AdditionalgradingcriteriaALL if Initial_Confidence_GradingWHO=="3) Uncertain significance" & Initial_Confidence_GradingALL!="3) Uncertain significance" & evidence!="NA"
replace Final_Confidence_Grading=Final_Confidence_GradingALL if Initial_Confidence_GradingWHO=="3) Uncertain significance" & Initial_Confidence_GradingALL!="3) Uncertain significance" & evidence!="NA"

* rule 3
replace evidence="WHO" if Initial_Confidence_GradingWHO!="3) Uncertain significance" & Initial_Confidence_GradingALL=="3) Uncertain significance" & evidence!="NA"
replace Initial_Confidence_Grading=Initial_Confidence_GradingWHO if Initial_Confidence_GradingWHO!="3) Uncertain significance" & Initial_Confidence_GradingALL=="3) Uncertain significance"
replace Additionalgradingcriteria=AdditionalgradingcriteriaWHO if Initial_Confidence_GradingWHO!="3) Uncertain significance" & Initial_Confidence_GradingALL=="3) Uncertain significance"
replace Final_Confidence_Grading=Final_Confidence_GradingWHO  if Initial_Confidence_GradingWHO!="3) Uncertain significance" & Initial_Confidence_GradingALL=="3) Uncertain significance"

* rule 4
replace evidence="WHO" if Initial_Confidence_GradingWHO=="2) Assoc w R - Interim" & Initial_Confidence_GradingALL=="1) Assoc w R"
replace Initial_Confidence_Grading="2) Assoc w R - Interim"  if  Initial_Confidence_GradingWHO=="2) Assoc w R - Interim" & Initial_Confidence_GradingALL=="1) Assoc w R"
replace Additionalgradingcriteria="Downgraded to interim based on WHO dataset" if  Initial_Confidence_GradingWHO=="2) Assoc w R - Interim" & Initial_Confidence_GradingALL=="1) Assoc w R"
replace Final_Confidence_Grading=Final_Confidence_GradingWHO  if  Initial_Confidence_GradingWHO=="2) Assoc w R - Interim" & Initial_Confidence_GradingALL=="1) Assoc w R"

replace evidence="WHO" if Initial_Confidence_GradingWHO=="1) Assoc w R" & Initial_Confidence_GradingALL=="2) Assoc w R - Interim"
replace Initial_Confidence_Grading="1) Assoc w R"  if  Initial_Confidence_GradingWHO=="1) Assoc w R" & Initial_Confidence_GradingALL=="2) Assoc w R - Interim"
replace Additionalgradingcriteria=AdditionalgradingcriteriaWHO if Initial_Confidence_GradingWHO=="1) Assoc w R" & Initial_Confidence_GradingALL=="2) Assoc w R - Interim"
replace Final_Confidence_Grading=Final_Confidence_GradingWHO if Initial_Confidence_GradingWHO=="1) Assoc w R" & Initial_Confidence_GradingALL=="2) Assoc w R - Interim"

* rule 5 
replace evidence="FLAG for manual check" if Initial_Confidence_GradingWHO=="4) Not assoc w R - Interim" & inlist(Initial_Confidence_GradingALL,"1) Assoc w R","2) Assoc w R - Interim")
replace Initial_Confidence_Grading="6) Manual check" if Initial_Confidence_GradingWHO=="4) Not assoc w R - Interim" & inlist(Initial_Confidence_GradingALL,"1) Assoc w R","2) Assoc w R - Interim")                               
replace Additionalgradingcriteria="FLAG for manual check" if Initial_Confidence_GradingWHO=="4) Not assoc w R - Interim" & inlist(Initial_Confidence_GradingALL,"1) Assoc w R","2) Assoc w R - Interim")
replace Final_Confidence_Grading="6) Manual check" if Initial_Confidence_GradingWHO=="4) Not assoc w R - Interim" & inlist(Initial_Confidence_GradingALL,"1) Assoc w R","2) Assoc w R - Interim")


rename variant variant_category_v2
merge m:1 drug variant_category_v2 using  "new_variant_matched_to_old.dta"

gen Present_in_Catalogue_v1=1 if _m!=1
replace Present_in_Catalogue_v1=0 if Present_in_Catalogue_v1!=1
drop _m variant_category_v1 final_grading_v1

rename variant_category_v2 variant

merge m:1 variant using "total_variant_counts_including_orphans.dta"
drop if _m==2
drop _m

* finalise Additionalgradingcriteria
replace Additionalgradingcriteria="Silent mutation" if inlist(effectALL,"synonymous_variant","initiator_codon_variant","stop_retained_variant") & strpos(Additionalgradingcriteria,"Literature")==0


save "Final_graded_algorithm_catalogue_primary_analysis_`y'_`c(current_date)'.dta", replace

export excel using "output/Final_graded_algorithm_catalogue_primary_analysis_`y'_`c(current_date)'.xlsx", firstrow(variables) replace


*** export results for CC and CCATU:

foreach a in CC CCATU {

use "master_stats_`a'_i_`y'.dta", clear

rename variant variant_category_v2
merge m:1 drug variant_category_v2 using  "new_variant_matched_to_old.dta"

gen Present_in_Catalogue_v1=1 if _m!=1
replace Present_in_Catalogue_v1=0 if Present_in_Catalogue_v1!=1
drop _m variant_category_v1 final_grading_v1

rename variant_category_v2 variant

merge m:1 variant using "total_variant_counts_including_orphans.dta"
drop if _m==2
drop _m

save "Final_graded_algorithm_catalogue_`a'_`y'_`c(current_date)'.dta", replace

export excel using "output/Final_graded_algorithm_catalogue_`a'_`y'_`c(current_date)'.xlsx", firstrow(variables) replace
}
}



