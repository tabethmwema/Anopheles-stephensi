setwd("")
install.packages("see")
install.packages("vioplot")

library(scales)
library(ggplot2)
library(see)
library(dplyr)
library(vioplot)

####Import .csv data
idata = read.table("intervention.csv", header=T, sep=",")
bdata = read.table("bionomics.csv", header=T, sep=",")
head(bdata)
str(bdata) ###check dataframe

#separate into native
native = bdata[bdata$Origin=="n",,drop=FALSE]
head(native)

#separate resting methods
rest = native[native$Method=="r",,drop=FALSE]  #in the bdata dataset, when the method column is "r", put it in the rest dataset
head(rest)

#separate bite methods
bite = native[native$Method=="b",,drop=FALSE]
head(bite)

#separate into invasive
invasive = bdata[bdata$Origin=="i",,drop=FALSE]
head(invasive)

#separate resting methods
rest = invasive[invasive$Method=="r",,drop=FALSE] 
head(rest)

# separate bite methods
bite = invasive[invasive$Method=="b",,drop=FALSE]
head(bite)


#find adult collection only
rest$One.day.collection <- as.numeric(rest$One.day.collection) ####remove spaces and converts to numeric
badata = rest[rest$Mosquito.stage.reclassification!="Larvae",]
badata$logstnum = log(badata$One.day.collection)
badata$Condensed.collection.method[badata$Condensed.collection.method==""] = NA
badata$Condensed.collection.method[badata$Condensed.collection.method=="Aquatic larval collection"] = NA
badata = badata[!is.na(badata$Condensed.collection.method),] # Exclude rows with blanks (NA values) from column Condensed collection method

####Question 1. Do different methods of collecting An. stephensi have different sensitivity?####
#want to do pairwise comparisons of the methods, ending with a pairwise table summarizing results

#A. adult collection only
#use badata

#B. calculate means just to have something to sort by
methods = unique(badata$Condensed.collection.method)
methods = methods[methods!=""]
OUT = NULL
for(m in 1:length(methods)){
  t = badata[badata$Condensed.collection.method==as.character(methods[m]),]
  OUT = rbind(OUT, c(as.character(methods[m]), mean(t$One.day.collection), nrow(t)))
}
methods = OUT
methods = as.data.frame(methods)
colnames(methods) = c("method", "mean.col.rate", "nrecords")
methods = methods[rev(order(methods$mean.col.rate)),]

#C.calculate pairwise regressions, doing all comparisons
# Convert 'mean.col.rate' and 'nrecords' to numeric
methods$mean.col.rate <- as.numeric(methods$mean.col.rate) ###Only run for native bite
methods$nrecords <- as.numeric(methods$nrecords) ###Only run for native bite

opt1 = methods$method
opt2 = methods$method
comps = expand.grid(opt1, opt2)
comps = comps[comps$Var1!=comps$Var2,]
OUT = NULL
for(m in 1:length(unique(comps$Var2))){
  tt = comps[comps$Var2==as.character(unique(comps$Var2))[m],,drop=F]
  if((m)>nrow(tt)){next}
  tt = tt[(m):nrow(tt),,drop=F]
  OUT = rbind(OUT, tt)
}
comps = as.data.frame(OUT)

#D.iterate over each comparison for regression AND permute t value. Save needed output.
OUT = NULL
for(r in 1:nrow(comps)){
  #isolate data of interest
  t = rbind(badata[badata$Condensed.collection.method==as.character(comps$Var2[r]),,drop=F], badata[badata$Condensed.collection.method==as.character(comps$Var1[r]),,drop=F ])
  t = as.data.frame(t)
  
  #run regression
  model = lm(log(One.day.collection)~Condensed.collection.method, data=t)
  
  #save data 
  emp = c(summary(model)$coefficients[1,1], summary(model)$coefficients[2,1], summary(model)$coefficients[1,3], summary(model)$coefficients[2,3], summary(model)$df[2], summary(model)$r.squared) #actual intercept, slope, t value intercept, t value slope, DF, R2
  
  #permute T distribution to calculate p
  permvalues = NULL
  for(i in 1:1000){
    tt = data.frame(One.day.collection=t$One.day.collection, Condensed.collection.method=sample(t$Condensed.collection.method, replace=F, nrow(t)))
    tmodel = lm(log(One.day.collection)~Condensed.collection.method, data=tt)
    permvalues = rbind(permvalues, c(summary(tmodel)$coeff[1,3], summary(tmodel)$coeff[2,3]))
  }
  summary(tmodel)
  
  #calculate p values and save all the emp data as well as newly calculated p values
  OUT = rbind(OUT, c(as.character(comps$Var1[r]), as.character(comps$Var2[r]), emp, (nrow(permvalues[permvalues[,1]>abs(emp[1]),,drop=F])/i), (nrow(permvalues[permvalues[,2]>abs(emp[2]),,drop=F])/i)))
}
colnames(OUT) = c("method1", "method2", "interceptEst", "slopeEst", "interceptT", "slopeT", "df", "R2", "interceptP", "slopeP")
methodcomp = as.data.frame(OUT)
nrow(methodcomp[methodcomp$slopeP<0.05,])

#F. generate a basic a plot
pdata <- data.frame(y = log(badata$One.day.collection), x = badata$Condensed.collection.method)

pdf("invasive_rest1.pdf", width = 8, height = 4)
par(bg = NA)
ggplot(pdata, aes(x = x, y = y, fill = x)) +
  geom_violin() +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5, binwidth = 0.1, fill = "grey50") + 
  theme_minimal() +
  theme(legend.position = "none", panel.grid = element_blank(),
        axis.line = element_line(color = "black")     )  
dev.off()




####Question 2. Does An. stephensi bite more indoors or outdoors?####
#A. adult collection only
#important: start with filtered version of badata from question 1
#find adult collection only
#Filter based on origin and bite methods from question 1
bite$One.day.collection <- as.numeric(bite$One.day.collection) ####remove spaces and converts to numeric
badata = bite[bite$Mosquito.stage.reclassification!="Larvae",]
badata$logstnum = log(badata$One.day.collection)
badata$Condensed.collection.method[badata$Condensed.collection.method==""] = NA
badata$Condensed.collection.method[badata$Condensed.collection.method=="Aquatic larval collection"] = NA
badata = badata[!is.na(badata$Condensed.collection.method),] # Exclude rows with blanks (NA values) from column Condensed collection method


badata$Biting.indoor.or.outdoor[badata$Biting.indoor.or.outdoor==""] = NA
badata$Biting.indoor.or.outdoor = as.factor(badata$Biting.indoor.or.outdoor)
badata$Condensed.collection.method = as.factor(badata$Condensed.collection.method)
badata = badata[!is.na(badata$Biting.indoor.or.outdoor),] # Exclude rows with blanks (NA values) from column biting in/out

#B. calculate means just to have something to sort by
#don't need for this question

#C. extract residuals from collection method analysis
badata.c = badata[,c("One.day.collection","Condensed.collection.method","Biting.indoor.or.outdoor")]
badata.c = badata.c[complete.cases(badata.c),] ####Exclude blank rows from dataframe
badata.c$collresid <- rstandard(lm(log(One.day.collection)~Condensed.collection.method-1, data=badata.c))

#D. calculate regressions with permuted t values to estimate p values
model = lm(collresid~Biting.indoor.or.outdoor+Condensed.collection.method, data=badata.c)
summary(model)

#save data 
df = summary(model)$df 
r2 = summary(model)$r.squared
emp = data.frame(slope = summary(model)$coefficients[ ,1], t = summary(model)$coefficients[ ,3])#, reg.order = rep(c("Aspiration targeting resting mosquitoes", "Hand collection", "Host seeking baited trap", "Mechanical trap", "Passive trap", "Pyrethrum spray catch", "biting"), 1)) #, bite = c(rep("Indoor", nrow()), rep("Outdoor", 6)))
#rownames(emp) = seq(1:nrow(emp))  

#permute T distribution to calculate p
permvalues = NULL
t = badata.c
for(i in 1:1000){
  tt = data.frame(One.day.collection=sample(t$One.day.collection, nrow(t), replace=F), Condensed.collection.method=t$Condensed.collection.method, Biting.indoor.or.outdoor=t$Biting.indoor.or.outdoor)
  tmodel = lm(log(One.day.collection)~Biting.indoor.or.outdoor+Condensed.collection.method, data=tt)
  summary(tmodel)
  permvalues = rbind(permvalues, t(summary(tmodel)$coeff[ ,1]))
}

#calculate p values and save all the emp data as well as newly calculated p values
emp$p = rep(NA, nrow(emp))
for(r in 1:nrow(emp)){
  perms = permvalues[,r]
  emp$p[r] = length(perms[abs(perms)>abs(emp$t[r])])/length(perms)
}

#E. assess output
emp
emp[emp$p<0.05,]

###Collect residuals from the first model
rstandard(model) ####produces the residuals
cmodel = lm(log(One.day.collection)~Condensed.collection.method, data=badata)
badata$residCmodel = rstandard(cmodel)
iomodel = lm(residCmodel~Biting.indoor.or.outdoor, data=badata)
summary(iomodel)
confint(iomodel)

# Create the data frame
pdata <- data.frame(y = log(badata$One.day.collection), 
                    x = badata$Biting.indoor.or.outdoor)

# Plot
par(bg = NA)
ggplot(pdata, aes(x = x, y = y, fill = x)) +
  geom_violin() +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5, binwidth = 0.1, fill = "grey50") + 
  labs(x = "Biting Indoor or Outdoor", y = "Log(One.day.collection)", 
       title = "Violin Plot of Biting Indoor/Outdoor and Log of One.day.collection") +
  scale_fill_manual(values = c("Indoor" = "indianred", "Outdoor" = "cadetblue")) + # Specify colors
  theme_minimal() +
  theme(legend.position = "none", panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        plot.background = element_rect(fill = "white"))

# Save the plot
ggsave("native_restInOut.png", plot = last_plot(), width = 6, height = 4)



####Question 3. Does An. stephensi rest more indoors or outdoors?####
#A. adult collection only
#important: start with filtered version of badata from question 1
#filter based on origin and rest methods from question 1 
rest$One.day.collection <- as.numeric(rest$One.day.collection) ####remove spaces and converts to numeric
badata = rest[rest$Mosquito.stage.reclassification!="Larvae",]
badata$logstnum = log(badata$One.day.collection)
badata$Condensed.collection.method[badata$Condensed.collection.method==""] = NA
badata$Condensed.collection.method[badata$Condensed.collection.method=="Aquatic larval collection"] = NA
badata = badata[!is.na(badata$Condensed.collection.method),] # Exclude rows with blanks (NA values) from column Condensed collection method

badata$Resting.indoor.or.outdoor[badata$Resting.indoor.or.outdoor==""] = NA
badata$Biting.indoor.or.outdoor = as.factor(badata$Biting.indoor.or.outdoor)
badata$Condensed.collection.method = as.factor(badata$Condensed.collection.method)
badata = badata[!is.na(badata$Resting.indoor.or.outdoor),] # Exclude rows with blanks (NA values) from column resting in/out

#B. calculate means just to have something to sort by
#don't need for this question

#C. extract residuals from collection method analysis
badata.c2 = badata[,c("One.day.collection","Condensed.collection.method","Resting.indoor.or.outdoor")]
badata.c2 = badata.c2[complete.cases(badata.c2),]
badata.c2$collresid = rstandard(lm(log(One.day.collection)~Condensed.collection.method-1, data=badata.c2))

#D. calculate regressions with permuted t values to estimate p values
model = lm(collresid~Resting.indoor.or.outdoor+Condensed.collection.method, data=badata.c2)
summary(model)

#save data 
df = summary(model)$df 
r2 = summary(model)$r.squared
emp = data.frame(slope = summary(model)$coefficients[ ,1], t = summary(model)$coefficients[ ,3])#, reg.order = rep(c("Aspiration targeting resting mosquitoes", "Hand collection", "Host seeking baited trap", "Mechanical trap", "Passive trap", "Pyrethrum spray catch", "biting"), 1)) #, bite = c(rep("Indoor", nrow()), rep("Outdoor", 6)))
#rownames(emp) = seq(1:nrow(emp))  

#permute T distribution to calculate p
permvalues = NULL
t = badata.c2
for(i in 1:1000){
  tt = data.frame(One.day.collection=sample(t$One.day.collection, nrow(t), replace=F), Condensed.collection.method=t$Condensed.collection.method, Resting.indoor.or.outdoor=t$Resting.indoor.or.outdoor)
  tmodel = lm(log(One.day.collection)~Resting.indoor.or.outdoor+Condensed.collection.method, data=tt)
  summary(tmodel)
  permvalues = rbind(permvalues, t(summary(tmodel)$coeff[ ,1]))
}

#calculate p values and save all the emp data as well as newly calculated p values
emp$p = rep(NA, nrow(emp))
for(r in 1:nrow(emp)){
  perms = permvalues[,r]
  emp$p[r] = length(perms[abs(perms)>abs(emp$t[r])])/length(perms)
}

#E. assess output
emp
emp[emp$p<0.05,]

###Collect residuals from the first model
rstandard(model) ####produces the residuals
cmodel = lm(log(One.day.collection)~Condensed.collection.method, data=badata)
badata$residCmodel = rstandard(cmodel)
iomodel = lm(residCmodel~Biting.indoor.or.outdoor, data=badata)
summary(iomodel)
confint(iomodel)


#Question 3. What is the most preferred breeding site?
setwd("C:/Users/tzm0087/OneDrive - Auburn University/Documents/SCHOOL/PhD/CLASSES/FIRST YEAR/SECOND SEMESTER/META-ANALYSIS/ANOPHELES STEPHENSI/Data/Unchanged data")
install.packages("dplyr")
library(dplyr)

datumbreed = read.table("An. stephensi spreadsheet-11-24-22-bio.csv", header=T, sep=",")
head(datumbreed)

#Select columns of interest
databreed <- datumbreed %>%
  filter(Mosquito.stage.reclassification == "Larvae") %>%
  select(Origin, Condensed.breeding.site, Mosquito.stage.reclassification:Density.per.day)

databreed <- datumbreed[datumbreed$Mosquito.stage.reclassification == "Larvae", c("Condensed.breeding.site", "Density.per.day", "Condensed.collection.method", "Origin")]
databreed$Condensed.breeding.site[databreed$Condensed.breeding.site==""] = NA
databreed = databreed[!is.na(databreed$Condensed.breeding.site),] # Exclude rows with blanks (NA values) from column condensed breeding site
databreed$`Density per day`[is.na(databreed$`Density per day`)] = NULL
databreed = databreed[complete.cases(databreed[, "Density.per.day"]), ]
databreed <- databreed[databreed$Origin == 'i', ]


###Runs a pairwise regression
breed= glm(log(Density.per.day)~Condensed.breeding.site, data=databreed[databreed$Condensed.breeding.site%in%c("Drinking water reservoirs","Waste water"),], family="gaussian")
summary(breed)
exp(confint(breed)) #CI
exp(coef(breed)) #Coeff
exp(summary(breed)$coefficients[, "Std. Error"]) #SE
coef(breed) / summary(breed)$coefficients[, "Std. Error"] #T values


###Runs permutation test
sim<-replicate(999,summary(glm(sample(log(Density.per.day),replace=F)~Condensed.breeding.site,data=databreed[databreed$Condensed.breeding.site%in%c("Drinking water reservoirs","Waste water"),],family="gaussian"))$coefficients[2,3])
#hist(sim)
sum(abs(sim)>abs(summary(breed)$coefficients[2,3]))/1000 #p-value

#Question 4: Which An. stephensi control interventions are effective for adults and larvae?
setwd("")
####datuminter = read.csv(file.choose()) ####An. stephensi spreadsheet_11_04_22 ####Carbamate
intervention = read.table("intervention.csv", header=T, sep=",")

library(dplyr)
library(tidyverse)
library(ggplot2)
library(dunn.test)

head(intervention)

#Adults
# Filter data for adults and by origin (native "n", invasive "i"), 
adult_data <- intervention %>%
  filter(
    Mosquito.stage.reclassification == "Adult" &
      trimws(Origin) == "i"  # change to "n" for native
  )

# Count the number of rows for each insecticide
insecticide_counts <- adult_data %>%
  group_by(Insecticide.reclassification) %>%
  summarize(count = n())

# Filter data for insecticides with counts >= 2
filtered_data <- adult_data %>%
  filter(Insecticide.reclassification %in% insecticide_counts$Insecticide.reclassification[insecticide_counts$count >= 3])

# Count the number of rows for each unique insecticide
insecticide_counts_filtered <- filtered_data %>%
  count(Insecticide.reclassification)

# Scale mortality rate
filtered_data <- filtered_data %>%
  mutate(scaled_mortality_rate = scale(Mortality.Rate))

# Run Shapiro-Wilk normality test
shapiro_test_result <- shapiro.test(filtered_data$scaled_mortality_rate)

# Print the result
print(shapiro_test_result)

# Run Kruskal-Wallis test
kruskal_test_result <- kruskal.test(scaled_mortality_rate ~ Insecticide.reclassification, data = filtered_data)

# Print the result
print(kruskal_test_result)

# Perform Dunn test for pairwise comparisons
dunn_result <- dunn.test(filtered_data$scaled_mortality_rate, g = filtered_data$Insecticide.reclassification, method = "bonferroni")

# Print the results
print(dunn_result)



# Create a boxplot
ggplot(filtered_data, aes(x = Insecticide.reclassification, y = scaled_mortality_rate, fill = Insecticide.reclassification)) +
  geom_boxplot() +
  labs(x = "Insecticide Reclassification", y = "Scaled Mortality Rate", title = "Boxplot of Scaled Mortality Rate by Insecticide Reclassification") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_discrete(name = "Insecticide Reclassification")


#Larvae
#Biocontrol, larvicides, organic compounds
# Filter the intervention dataset for native/invasive
# Filter data based on the first set of conditions
first_filter <- intervention %>%
  filter(
    trimws(Insecticide.reclassification) %in% c("Biopesticide", "Organophosphate") & 
      trimws(Mosquito.stage.reclassification) == "Larvae" &
      (trimws(Insecticide) == "Aphanius dispar" |
         trimws(Insecticide) == "Aplocheilus blockii" |
         trimws(Insecticide) == "Bacillus sphaericus" |
         trimws(Insecticide) == "Bacillus thuringiensis israelensis" |
         trimws(Insecticide) == "Beauveria bassiana" |
         trimws(Insecticide) == "Gambusia affinis" |
         trimws(Insecticide) == "Poecilia reticulata" |
         trimws(Insecticide) == "Fenthion" |
         trimws(Insecticide) == "Temephos")
  )

# Filter the result of the first filter with the second set of conditions
organophosphate_biopesticide_larvae_data <- first_filter %>%
  filter(
    trimws(Origin) == "n"
  )

#For invasive range
first_filter <- intervention %>%
  filter(
    trimws(Insecticide.reclassification) %in% c("Biopesticide", "Organophosphate") & 
      trimws(Mosquito.stage.reclassification) == "Larvae" &
      (trimws(Insecticide) == "Bacillus thuringiensis israelensis" |
         trimws(Insecticide) == "Temephos")
  )

# Filter the result of the first filter with the second set of conditions
organophosphate_biopesticide_larvae_data <- first_filter %>%
  filter(
    trimws(Origin) == "i"
  )

#Print the filtered dataset
print(organophosphate_biopesticide_larvae_data)


# scale mortality rate 
scaled_mortality <- scale(organophosphate_biopesticide_larvae_data$Mortality.Rate)

# Create a new column with the scaled values
organophosphate_biopesticide_larvae_data$Scaled_Mortality <- scaled_mortality

#test for normality
shapiro_test <- shapiro.test(organophosphate_biopesticide_larvae_data$Scaled_Mortality)
shapiro_test

# Run Kruskal-Wallis test
kruskal_test <- kruskal.test(Scaled_Mortality ~ Insecticide, data = organophosphate_biopesticide_larvae_data)

# Print the results
print(kruskal_test)

# Run Dunn's test for pairwise comparisons
dunn_result <- dunn.test(organophosphate_biopesticide_larvae_data$Scaled_Mortality, g = organophosphate_biopesticide_larvae_data$Insecticide, method = "bonferroni")

# Print the results
print(dunn_result)

#create boxplot
# Define custom colors
custom_colors <- c("Aphanius dispar" = "#F8766D",
                   "Bacillus sphaericus" = "#CD9600",
                   "Bacillus thuringiensis israelensis" = "#ABA300",
                   "Beauveria bassiana" = "#0CB702",  
                   "Gambusia affinis" = "#00B8E7",
                   "Poecilia reticulata" = "#8494FF",
                   "Temephos" = "#AA4499")


# Get unique values of the "Insecticide" column
unique_insecticides <- unique(organophosphate_biopesticide_larvae_data$Insecticide)

# Match the custom colors to the unique insecticides
colors <- custom_colors[organophosphate_biopesticide_larvae_data$Insecticide]

# Create the boxplot
boxplot <- ggplot(organophosphate_biopesticide_larvae_data, aes(x = Insecticide, y = Scaled_Mortality, fill = Insecticide)) +
  geom_boxplot(fill = NA, aes(colour = Insecticide), position = position_dodge(0.8), width = 0.7, coef = 1.5) + # Hollow box with custom colors
  labs(title = "Effect of Different Insecticides on Larval Mortality Rate",
       x = "Insecticide",   # add x-axis label
       y = "Median Mortality Rate") + # add y-axis label
  theme_minimal() +
  theme(legend.position = "none",           # remove legend
        axis.text.x = element_text(angle = 45, hjust = 1),  # rotate x-axis labels
        #panel.border = element_rect(colour = "black", fill = NA),  # add border around plot
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        panel.grid.major = element_blank(),  # remove major grid lines
        panel.grid.minor = element_blank())  # remove minor grid lines)  # add axis lines


print(boxplot)
# Save as PDF
# Save as PDF with increased size
ggsave("Native larvae control.pdf", plot = boxplot, device = "pdf", width = 10, height = 6)


