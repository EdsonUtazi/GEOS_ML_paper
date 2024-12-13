#Set working directory
setwd("//worldpop.files.soton.ac.uk/worldpop/Projects/WP516374_VaxPop/Working/DHS_2018_data/NGKR7ADT")

#-------------Read in data---------------------#
data1 <- read.csv("Vax_NG_measles_proc.csv", header = T)


length(unique(data1$Cluster_no))
max(unique(data1$Cluster_no))

#Find missing clusters in data1
clust <- 1:max(unique(data1$Cluster_no))    
miss <- clust[!clust %in% unique(data1$Cluster_no)] #missing clusters in stata file
miss
clust <- clust[-miss]


#------------------------Data extraction starts here-------------------------------------#
#Summarize data1 - total and by age <9, 9-11, 12-23, 24-35
#n.clust = length(unique(data1$Cluster_no))

#subset to age 0-35
#data1 <- subset(data1, childs_age_in_months <= 35)

n.clust = length(clust)
agebrks = c(0,8,11,23,35)

#Unique cluster nos
#Clust <- unique(data1$Cluster_no)
Clust <- clust

clust.tab.1 = clust.tab.2 = matrix(0, nrow = n.clust, ncol = 4)
for (i in 1:n.clust)
{
  print(i)
  subdat.1 = subset(data1, Cluster_no==Clust[i], select=c(childs_age_in_months, received_measles))
  for (j in 1:(length(agebrks)-1)){
    if (j>1) subdat.2 = subset(subdat.1, childs_age_in_months>agebrks[j] & childs_age_in_months<=agebrks[j+1])
    if (j==1) subdat.2 = subset(subdat.1, childs_age_in_months>=agebrks[j] & childs_age_in_months<=agebrks[j+1])
    print(j)
    clust.tab.1[i,j] = nrow(subdat.2)
    clust.tab.2[i,j] = length(subdat.2[subdat.2[,2]==1,2])
  }
}

#Add column headings
#Totals ind. age groups
colnames(clust.tab.1) = c("0-8", "9-11", "12-23", "24-35")

#Numbers vaccinated ind. age groups
colnames(clust.tab.2) = c("0-8", "9-11", "12-23", "24-35")
head(clust.tab.1); head(clust.tab.2)


#Total counts for 0-35
clust.tab.3 = matrix(0, nrow = n.clust, ncol = 3)
for (i in 1:n.clust) #NOTE 35
{
  sub.dat3 = subset(data1, Cluster_no==Clust[i] & childs_age_in_months <= 35, select = c(Cluster_no, childs_age_in_months,received_measles))
  
  Totchild = nrow(sub.dat3); Numvacc = length(sub.dat3$received_measles[sub.dat3$received_measles==1])
  #print(Totchild); print(Numvacc)
  clust.tab.3[i,]= c(as.numeric(sub.dat3[1,1]), Totchild, Numvacc) 
}
head(clust.tab.3)

colnames(clust.tab.3)=c("cluster", "Totchild", "Numvacc")

all.var=data.frame(1:length(clust.tab.1[,1]))
all.var$Age0_8_Tot = clust.tab.1[,1]; all.var$Age0_8_Vax = clust.tab.2[,1]
all.var$Age9_11_Tot = clust.tab.1[,2]; all.var$Age9_11_Vax = clust.tab.2[,2]
all.var$Age12_23_Tot = clust.tab.1[,3]; all.var$Age12_23_Vax = clust.tab.2[,3]
all.var$Age24_35_Tot = clust.tab.1[,4]; all.var$Age24_35_Vax = clust.tab.2[,4]
all.var$Totchild = clust.tab.3[,2]; all.var$Numvacc = clust.tab.3[,3]
all.var$Age9_35_Tot = all.var$Totchild - all.var$Age0_8_Tot
all.var$Age9_35_Vax = all.var$Numvacc - all.var$Age0_8_Vax
all.var$DHSCLUST = Clust
all.var = all.var[,-1]

head(all.var)
#-----------------------End of data extraction--------------------------------#

#------Write output file-------------#
write.csv(all.var, "NG_Vax_measles.csv")