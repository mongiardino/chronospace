library(devtools)
load_all()
#devtools::install_github("YuLab-SMU/ggtree")

#names0 <- list.files("C:/Users/pablo/Desktop/Trabajo/Investigacion/Proyectos/chronospace_local/data/curcu")
names0 <- list.files("C:/Users/pablo/OneDrive/Escritorio/chrospace data/curcu")

names1 <- gsub(names0, pattern = ".tre", replacement = "")
type0 <- sapply(strsplit(names1, split = "_"), function(x) {x[2:5]})
type <- list(type0[1,], type0[2,], type0[3,], type0[4,])
sample <- 500

curcu_data <- extract_ages(path="C:/Users/pablo/OneDrive/Escritorio/chrospace data/curcu",
                           type, sample, facnames = c("clock", "prior", "model", "genes"))
#save(curcu_data, file="C:/Users/pablo/Desktop/Trabajo/Investigacion/Proyectos/chronospace_local/data/curcu data")


#names0 <- list.files("C:/Users/pablo/Desktop/Trabajo/Investigacion/Proyectos/chronospace_local/data/deca")
names0 <- list.files("C:/Users/pablo/OneDrive/Escritorio/chrospace data/deca")

names1 <- gsub(names0, pattern = ".tre", replacement = "")
type0 <- sapply(strsplit(names1, split = "_"), function(x) {x[2:5]})
type <- list(type0[1,], type0[2,], type0[3,], type0[4,])
sample <- 500

deca_data <- extract_ages(path="C:/Users/pablo/OneDrive/Escritorio/chrospace data/deca",
                          type, sample, facnames = c("clock", "prior", "model", "genes"))
#save(deca_data, file="C:/Users/pablo/Desktop/Trabajo/Investigacion/Proyectos/chronospace_local/data/deca data")


#names0 <- list.files("C:/Users/pablo/Desktop/Trabajo/Investigacion/Proyectos/chronospace_local/data/euka")
names0 <- list.files("C:/Users/pablo/OneDrive/Escritorio/chrospace data/euka")

names1 <- gsub(names0, pattern = ".tre", replacement = "")
type0 <- sapply(strsplit(names1, split = "_"), function(x) {x[2:5]})
type <- list(type0[1,], type0[2,], type0[3,], type0[4,])
sample <- 500

euka_data <- extract_ages(path="C:/Users/pablo/OneDrive/Escritorio/chrospace data/euka",
                          type, sample, facnames = c("clock", "prior", "model", "genes"))
#save(euka_data, file="C:/Users/pablo/Desktop/Trabajo/Investigacion/Proyectos/chronospace_local/data/euka data")


#########################################################

document()
load_all()
# Curculioidea #################

#load data
#load("C:/Users/pablo/OneDrive/Escritorio/Software/chronospace/analysis/curcu data")

#general overview of the dataset
print(curcu_data)
head(curcu_data$factors)

#generate chronospace
cspace_curcu <- chronospace(data_ages = curcu_data)
cspace_curcu

#plot chronopace
plot(cspace, sdev = .85, output = "ordination", centroids = TRUE, distances = TRUE, ellipses = FALSE)
xxx <- plot(cspace_curcu, sdev = .85, output = "ordination", factor = 4, axes = 1, centroids = TRUE, ellipses = FALSE, distances = TRUE, pt.alpha = 0, dist.width = 1)
xxx$genes$ordination

xxx <- plot(cspace_curcu, output = "ordination", factor = c(4), distances = TRUE, pt.alpha = 1)





#ordin <- plot(cspace_curcu, sdev = 2)

ordin$clock$ordination
ordin$clock$PC_extremes

ordin$prior$ordination
ordin$prior$PC_extremes

ordin$genes$ordination
ordin$genes$PC_extremes$bgPC1


#Factor A explains ~6% of total variation, ~2% of non-redundant variation
ordin$factor_A$PC_extremes$bgPC1

#Factor D explains ~40% of total variation, ~13% of non-redundant variation
ordin$factor_D$PC_extremes$bgPC1
ordin$factor_D$PC_extremes$bgPC2
ordin2$factor_D$PC_extremes$bgPC3
ordin2$factor_D$PC_extremes$bgPC4

#Neither Factor B or C accounts for more than 1% of total variation

#get first 5 sensitive nodes for each factor, plot for A and D
load_all()
sn_curcu2 <- sensitive_nodes(curcu_data)
sn_curcu2$genes
sn_curcu2
sn_curcu$factor_D

#get LTT plots for each factor, plot for A and D
ltt_curcu <- ltt_sensitivity(curcu_data, factor = c(4))
ltt_curcu$genes

# Decapoda #################

#load data
load("C:/Users/pablo/OneDrive/Escritorio/Software/chronospace/analysis/deca data")

#general overview of the dataset
print(deca_data)

#generate chronospace
cspace_deca <- chronospace(data_ages = deca_data)
#cspace_deca <- chronospace2(data_ages = deca_data, facnames = c("clock", "prior", "model", "genes"))


ordin <- plot(cspace_deca, sdev = .8)

ordin$genes$PC_extremes$bgPC1


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
load("C:/Users/pablo/OneDrive/Escritorio/Software/chronospace/analysis/euka data")


#general overview of the dataset
print(euka_data)

#generate chronospace
cspace_euka <- chronospace(data_ages = euka_data)
#cspace_euka <- chronospace2(data_ages = euka_data, facnames = c("clock", "prior", "model", "genes"))

ordin <- plot(cspace_euka)
ordin$genes$PC_extremes$bgPC1

chronospace
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





