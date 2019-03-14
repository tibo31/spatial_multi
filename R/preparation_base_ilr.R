# this file is for preparing the modeling step
require("compositions")
require("sp")
# we load the data basis
load(file("http://www.thibault.laurent.free.fr/code/jrfm/R/BdD_final_dep.RData"))

# we drop the departement of Paris which is not concerned by the election
dep.2015.spdf <- dep.2015.spdf[dep.2015.spdf@data$CODE_DEPT!="75", ]

# we define the ilr compositions
# Important remark : depending on the order of the compositions in the function ilr(), 
# the ilr transformations focus on different comparaisons of compositions  
y_ilr <- ilr(dep.2015.spdf@data[, c("left" , "right", "others")]) 
# y1 <- ilr(dep.2015.spdf@data[, c("others", "left" , "right")])
# y2 <- ilr(dep.2015.spdf@data[, c("left" , "others", "right")])

# we add the variables to the spatial data
dep.2015.spdf@data$ilr_y1 <- y_ilr[,1]
dep.2015.spdf@data$ilr_y2 <- y_ilr[,2]


# define the matrix of transformation
V_age <- rbind(c(-1/2/sqrt(5), -1/2/sqrt(5), -1/2/sqrt(5), -1/2/sqrt(5), 2/sqrt(5)),
              c(-sqrt(3)/6, -sqrt(3)/6, -sqrt(3)/6, sqrt(3/4), 0),
               c(-sqrt(2/3)/2, -sqrt(2/3)/2,  sqrt(2/3), 0, 0),
               c(sqrt(1/2), -sqrt(1/2), 0, 0, 0))
  

V_dip <- rbind(c(-sqrt(3)/6, -sqrt(3)/6, -sqrt(3)/6, sqrt(3/4)),
               c(-sqrt(2/3)/2, -sqrt(2/3)/2,  sqrt(2/3), 0),
               c(sqrt(1/2), -sqrt(1/2), 0, 0))

V_employ <- rbind(c(-1/2/sqrt(5), -1/2/sqrt(5), -1/2/sqrt(5), -1/2/sqrt(5), 2/sqrt(5)),
               c(-sqrt(3)/6, -sqrt(3)/6, -sqrt(3)/6, sqrt(3/4), 0),
               c(-sqrt(2/3)/2, -sqrt(2/3)/2,  sqrt(2/3), 0, 0),
               c(sqrt(1/2), -sqrt(1/2), 0, 0, 0))

# we create the independant variabes 
# age without "mineur"
age <- dep.2015.spdf@data[, c("age_1824", "age_2540", "age_4055", "age_5564", "age_65")]
age_1840 = dep.2015.spdf@data[, "age_1824"] + dep.2015.spdf@data[, "age_2540"]
age_4064 = dep.2015.spdf@data[, "age_4055"] + dep.2015.spdf@data[, "age_5564"]
age_65 = dep.2015.spdf@data[, "age_65"]
age3 = as.data.frame(cbind(age_1840, age_4064, age_65))

#age_ratio <- log(dep.2015.spdf@data[, c("age_1824", "age_2564")]/dep.2015.spdf@data[, "age_65"])
#cor(cbind(y_ilr, age_ratio))

diplome <- dep.2015.spdf@data[, c("no_diplom", "capbep", "bac", "diplom_sup")]
No_CAPBEP = dep.2015.spdf@data[, "no_diplom"] + dep.2015.spdf@data[, "capbep"]
#diplome3 = cbind(dep.2015.spdf@data[, c("no_diplom", "capbep")],BACSUP)
diplome3 = cbind(No_CAPBEP, dep.2015.spdf@data[, c("bac", "diplom_sup")])

# ternaryDiag(diplome3, col = "royalblue", cex = 1, pch = 16)
# barplot(acomp(colMeans(diplome3)),col = 1:3, main = "Diplome",srt=45,
#         names.arg = c("BAC (16.8%)", "SUP (26.7%)", "No_CAPBEP (56.5%)"))
V_dip3 = rbind(c(1/sqrt(6), 1/sqrt(6), -2/sqrt(6)),
               c(1/sqrt(2), -1/sqrt(2),  0))

# euro.2014 <- dep.2015.spdf@data[, c("left_euro.2014", "right_euro.2014", "others_euro.2014")]

employ <- dep.2015.spdf@data[, c("AZ", "BE", "FZ", "GU", "OQ")]

# we create the matrix of the independent variables
x <- cbind(ilr(age, V = t(V_age)), 
           ilr(diplome, V = t(V_dip)), 
           ilr(employ, V = t(V_employ)),
           as.matrix(dep.2015.spdf@data[, c("unemp_rate", "employ_evol", 
                                                "owner_rate", "income_rate", 
                                                "foreign")]))

colnames(x) <- c("age_ilr1", "age_ilr2", "age_ilr3", "age_ilr4", 
                 "diplome_ilr1", "diplome_ilr2", "diplome_ilr3",
                 "employ_ilr1", "employ_ilr2", "employ_ilr3", "employ_ilr4", 
                 "unemp_rate", "employ_evol", "owner_rate", "income_rate", "foreign")


x_df <- data.frame(x)



x1 = cbind(ilr(age3), 
           ilr(diplome, V = t(V_dip)), 
           ilr(employ, V = t(V_employ)),
           as.matrix(dep.2015.spdf@data[, c("unemp_rate", "employ_evol", 
                                            "owner_rate", "income_rate", 
                                            "foreign")]))

colnames(x1) <- c("age3_ilr1", "age3_ilr2",  
                 "diplome_ilr1", "diplome_ilr2", "diplome_ilr3",
                 "employ_ilr1", "employ_ilr2", "employ_ilr3", "employ_ilr4", 
                 "unemp_rate", "employ_evol", "owner_rate", "income_rate", "foreign")

x1_df <- data.frame(x1)

x2 = cbind(ilr(age3), 
           ilr(diplome3, V = t(V_dip3)), 
           ilr(employ, V = t(V_employ)),
           as.matrix(dep.2015.spdf@data[, c("unemp_rate", "employ_evol", 
                                            "owner_rate", "income_rate", 
                                            "foreign")]))

colnames(x2) <- c("age3_ilr1", "age3_ilr2",  
                  "diplome3_ilr1", "diplome3_ilr2",
                  "employ_ilr1", "employ_ilr2", "employ_ilr3", "employ_ilr4", 
                  "unemp_rate", "employ_evol", "owner_rate", "income_rate", "foreign")

x2_df <- data.frame(x2)

rm(age3, age, diplome, x, x2,  age_1840, age_4064, age_65)



