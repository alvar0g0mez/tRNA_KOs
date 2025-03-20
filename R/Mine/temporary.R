# Repeat initial proteomics analysis from Alexis in my own way and with non-batch corrected data
# TEMPORARY, TO THEN BE INCORPORATED TO THE MAIN CODE


resp_old <- as.data.frame(fread("C:/MyStuff/tRNA_KOs/Data/responsiveness.csv"))
resp_new <- as.data.frame(fread("C:/MyStuff/tRNA_KOs/Data/responsiveness_batchcorrected_onWTs.csv"))


sum(resp_old$nDEP)
sum(resp_new$nDEP)






























