
clear

import delimited using "run-1682402668815-part-r-00000.csv", varnames(1) clear

gen variant_category_v2=resol+"_"+variant_cat
rename description variant_category_v1  
keep variant_category_v2 variant_category_v1 predicted_effect

* old_variant is currently strL rather than str#. The merge won't work with strL, so need to do the following:

generate str name_string = variant_category_v1 
replace variant_category_v1  = ""
compress variant_category_v1 
replace variant_category_v1  = name_string
drop name_string

* this has now transformed old_variant from strL to str#

save "v1_matching.dta", replace


use "v1_matching.dta", clear
drop if inlist(predicted_effect,"initiator_codon_variant","stop_retained_variant","synonymous_variant")

* now there is only one row per old_variant, and this list inlcudes 'start_lost' mutations, which include the V1V mutations identified in "v1_drug_mutation_grade.dta"

save "v1_matching_nonsyn.dta", replace


clear


import excel "List_of_neutral_mutations_from_cat_ver_1_rev20Jan2023.xlsx", sheet("Dataset_FINAL") cellrange(A3) allstring

save "neutral_mutations_catalogue_v1.dta", replace

clear

* Importing correctly the two line header

import excel "List_of_neutral_mutations_from_cat_ver_1_rev20Jan2023.xlsx", cellrange(A2:W2) sheet("Dataset_FINAL") allstring

save "header_first_line.dta", replace

clear

import excel "List_of_neutral_mutations_from_cat_ver_1_rev20Jan2023.xlsx", cellrange(A1:W1) sheet("Dataset_FINAL") allstring

append using "header_first_line.dta"

foreach var of varlist * {
    quietly replace `var'=`var'[1]+`var'[2] in 1
}

drop in 2


append using "neutral_mutations_catalogue_v1.dta"

foreach var of varlist * {
   rename `var' `=strtoname(`var'[1])'
}

drop in 1

replace drug="Ethambutol" if drug=="EMB" 
replace drug="Isoniazid" if drug=="INH"
replace drug="Rifampicin" if drug=="RIF"  
replace drug="Pyrazinamide" if drug=="PZA"  
replace drug="Moxifloxacin" if drug=="MXF"  
replace drug="Levofloxacin" if drug=="LEV"  
replace drug="Ethionamide" if drug=="ETH"  
replace drug="Amikacin" if drug=="AMI"  
replace drug="Capreomycin" if drug=="CAP"  
replace drug="Kanamycin" if drug=="KAN"  
replace drug="Streptomycin" if drug=="STM"  
replace drug="Clofazimine" if drug=="CFZ"  
replace drug="Bedaquiline" if drug=="BDQ"  
replace drug="Linezolid" if drug=="LZD"  
replace drug="Delamanid" if drug=="DLM"  

drop if strpos(FOR_WHO_CATALOGUE_UPDATE_VER__2, "NOT")

keep drug variant

rename variant variant_category_v1

merge m:1 variant_category_v1 using "sacha_v1_matching_nonsyn.dta"
drop if _m==2

drop if _m==1 

drop _m
keep drug variant_category_v2
rename variant_category_v2 variant
save "neutral_mutations_catalogue_v1.dta", replace


clear
import excel "WHO-UCN-GTB-PCI-2021.7-eng.xlsx", sheet("Mutation_catalogue") firstrow
drop if FINALCONFIDENCEGRADING=="combo"
rename variantcommon_name variant
replace drug="Ethambutol" if drug=="EMB" 
replace drug="Isoniazid" if drug=="INH"
replace drug="Rifampicin" if drug=="RIF"  
replace drug="Pyrazinamide" if drug=="PZA"  
replace drug="Moxifloxacin" if drug=="MXF"  
replace drug="Levofloxacin" if drug=="LEV"  
replace drug="Ethionamide" if drug=="ETH"  
replace drug="Amikacin" if drug=="AMI"  
replace drug="Capreomycin" if drug=="CAP"  
replace drug="Kanamycin" if drug=="KAN"  
replace drug="Streptomycin" if drug=="STM"  
replace drug="Clofazimine" if drug=="CFZ"  
replace drug="Bedaquiline" if drug=="BDQ"  
replace drug="Linezolid" if drug=="LZD"  
replace drug="Delamanid" if drug=="DLM"  
drop Genome 
split variant, p(" ")
replace variant=variant1
drop variant1 variant2 variant3
drop if variant==""
rename variant variant_category_v1
rename FINAL final_grading_v1
rename Additionalgradingcriteria Additionalgradingcriteria_v1

replace Additionalgradingcriteria_v1="" if inlist(Additionalgradingcriteria_v1, "WHO-endorsed gDST assay")

save "v1_drug_mutation_grade.dta", replace

clear
import excel "WHO_endorsed_assay_list_variants_24Apr2023.xlsx", firstrow
rename Drug drug
rename Variant variant
replace drug="Ethambutol" if drug=="EMB" 
replace drug="Isoniazid" if drug=="INH"
replace drug="Rifampicin" if drug=="RIF"  
replace drug="Pyrazinamide" if drug=="PZA"  
replace drug="Moxifloxacin" if drug=="MXF"  
replace drug="Levofloxacin" if drug=="LEV"  
replace drug="Ethionamide" if drug=="ETH"  
replace drug="Amikacin" if drug=="AMI"  
replace drug="Capreomycin" if drug=="CAP"  
replace drug="Kanamycin" if drug=="KAN"  
replace drug="Streptomycin" if drug=="STM"  
replace drug="Clofazimine" if drug=="CFZ"  
replace drug="Bedaquiline" if drug=="BDQ"  
replace drug="Linezolid" if drug=="LZD"  
replace drug="Delamanid" if drug=="DLM"  
save "assay_mutations.dta", replace

use "assay_mutations.dta", clear
rename variant variant_category_v2
gen strL Additionalgradingcriteria_v2="WHO-endorsed gDST assay"
save "v2_WHO_gDST_endorsed.dta", replace

clear

* merge v1_drug_mutations_grade.dta with list of non-synonymous mutations. There are multiple rows per old_variant in the v1_drug_mutation_grade.dta, but these are unique when paired with 'drug'
use "v1_drug_mutation_grade.dta", clear
merge m:1 variant_category_v1 using "v1_matching_nonsyn.dta"

drop if strpos(variant_category_v1,"fprA") & _m==1

* there are rpsL promotor mutations more than 234 bp upstream, so drop these (all such mutations that are _m==1 are >234 bp upstream):
drop if strpos(variant_category_v1,"rpsL") & _m==1 & strpos(variant_category_v1 ,"del")==0 & strpos(variant_category_v1 ,"ins")==0

* same for embR where promotor is now only 103bp long
drop if strpos(variant_category_v1,"embR") & _m==1 & strpos(variant_category_v1 ,"del")==0 & strpos(variant_category_v1 ,"ins")==0

* same for rrs where promotor is now only 151bp long
drop if strpos(variant_category_v1,"rrs") & _m==1 & strpos(variant_category_v1 ,"del")==0 & strpos(variant_category_v1 ,"ins")==0

* same for embC where promotor is now only 1982bp long
drop if strpos(variant_category_v1,"embC") & _m==1 & strpos(variant_category_v1 ,"del")==0 & strpos(variant_category_v1 ,"ins")==0

* same for ubiA where promotor is now only 1982bp long, but don't yet drop M1V
drop if strpos(variant_category_v1,"ubiA") & _m==1 & strpos(variant_category_v1 ,"del")==0 & strpos(variant_category_v1 ,"ins")==0 & variant_category_v1!="ubiA_M1V"

* drop all other promotor, indel variants:
drop if (strpos(variant_category_v1,"-") | strpos(variant_category_v1,"del") | strpos(variant_category_v1,"ins")) & _m==1

* drop other initiator_codon_variant (all are uncertain)
drop if regexm(variant_category_v1, "[MV]1[A-Z]") & _merge==1

* drop gid_V110G as this was an error in v1:
drop if variant_category_v1=="gid_V110G" | variant_category_v1=="inhA_T4I"

* manually annotate the remaining pncA mutations (there are 4 of these):
replace variant_category_v2="pncA_p.His71Asp" if variant_category_v1=="pncA_H71D"
replace variant_category_v2="pncA_p.Ile31Thr" if variant_category_v1=="pncA_I31T"
replace variant_category_v2="pncA_p.Leu116Arg" if variant_category_v1=="pncA_L116R"
replace variant_category_v2="pncA_p.Thr135Asn" if variant_category_v1=="pncA_T135N"
* Remove all v2 variants that have more than one match in v1.
bysort drug variant_category_v2: gen n=_N
keep if n==1
drop n

* now 'drug + variant' is unique and can be merged into the main dataset
* save version with old_variant and new variant and v1_confidence grading:
preserve
keep drug variant_category_v1 variant_category_v2 final_grading_v1
save "new_variant_matched_to_old.dta", replace
restore

drop _merge

merge 1:1 variant_category_v2 drug using "v2_WHO_gDST_endorsed.dta"
drop if _m==2

replace Additionalgradingcriteria_v1="WHO-endorsed gDST assay" if inlist(Additionalgradingcriteria_v2, "WHO-endorsed gDST assay")

rename variant_category_v2 variant

preserve
keep drug final_grading_v1 variant INITIALCONFIDENCEGRADING
rename INITIALCONFIDENCEGRADING Initial_Confidence_Grading
rename final_grading_v1 Final_Confidence_Grading
save "v1_grades_for_sens_spec.dta", replace
restore

keep drug Additionalgradingcriteria_v1 final_grading_v1 variant
save "v1_grades.dta", replace

