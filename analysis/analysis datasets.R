library(devtools)
load_all()

# Curculioidea #################

#load data
load("curcu data")

#general overview of the dataset
print(curcu_data)

#generate chronospace
cspace_curcu <- chronospace(data_ages = curcu_data, vartype = "non-redundant")

ordin <- plot(cspace_curcu, axes = )
ordin2 <- plot(cspace_curcu, axes = c(3,4))
ordin$factor_A$ordination
ordin$factor_B$ordination
ordin$factor_C$ordination
ordin$factor_D$ordination

cspace_curcu

#Factor C explains 1% of total variation
ordin$factor_C$PC_extremes$bgPC1

#Factor D explains 60% of total variation, almost 30% of non-redundant variaiton
ordin$factor_D$PC_extremes$bgPC1
ordin$factor_D$PC_extremes$bgPC2
ordin2$factor_D$PC_extremes$bgPC3
ordin2$factor_D$PC_extremes$bgPC4

#get first 5 sensitive nodes for each factor, plot for C and D
sn_curcu <- sensitive_nodes(obj = cspace_curcu)
sn_curcu$factor_C
sn_curcu$factor_D

#get LTT plots for each factor, plot for C and D
ltt_curcu <- ltt_sensitivity(curcu_data)
ltt_curcu$factor_C
ltt_curcu$factor_D


# Decapoda #################

#load data
load("deca data")

#general overview of the dataset
print(deca_data)

#generate chronospace
cspace_deca <- chronospace(data_ages = deca_data, vartype = "non-redundant")

ordin <- plot(cspace_deca, axes = )
ordin2 <- plot(cspace_deca, axes = c(3,4))
ordin$factor_A$ordination
ordin$factor_B$ordination
ordin$factor_C$ordination
ordin$factor_D$ordination

cspace_deca

#Factor C explains ~1% of total variation
ordin$factor_C$PC_extremes$bgPC1

#Factor D explains 50% of total variation, 30% of non-redundant variation
ordin$factor_D$PC_extremes$bgPC1
ordin$factor_D$PC_extremes$bgPC2
ordin2$factor_D$PC_extremes$bgPC3
ordin2$factor_D$PC_extremes$bgPC4

#get first 5 sensitive nodes for each factor, plot for C and D
sn_deca <- sensitive_nodes(obj = cspace_deca)
sn_deca$factor_C
sn_deca$factor_D

#get LTT plots for each factor, plot for C and D
ltt_deca <- ltt_sensitivity(deca_data)
ltt_deca$factor_C
ltt_deca$factor_D



# Eukariota #################

#load data
load("euka data")

#general overview of the dataset
print(euka_data)

#generate chronospace
cspace_euka <- chronospace(data_ages = euka_data, vartype = "non-redundant")

ordin <- plot(cspace_euka, axes = )
ordin2 <- plot(cspace_euka, axes = c(3,4))
ordin$factor_A$ordination
ordin$factor_B$ordination
ordin$factor_C$ordination
ordin$factor_D$ordination

cspace_euka

#Factor C explains ~3% of total variation
ordin$factor_C$PC_extremes$bgPC1

#Factor D explains almost 40% of total variation, but only ~4% of non-redundant variation
ordin$factor_D$PC_extremes$bgPC1
ordin$factor_D$PC_extremes$bgPC2
ordin2$factor_D$PC_extremes$bgPC3
ordin2$factor_D$PC_extremes$bgPC4

#get first 5 sensitive nodes for each factor, plot for C and D
sn_euka <- sensitive_nodes(obj = cspace_euka)
sn_euka$factor_C
sn_euka$factor_D

#get LTT plots for each factor, plot for C and D
ltt_euka <- ltt_sensitivity(euka_data)
ltt_euka$factor_C
ltt_euka$factor_D





