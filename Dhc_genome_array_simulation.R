strain195 <- c('AAW40468.1','AAW39240.1','AAW39053.1','AAW40373.1','AAW39215.1','AAW40589.1','AAW40575.1','AAW40356.1','AAW39273.1','AAW39060.1','AAW39843.1','AAW39605.1','AAW39229.1','AAW40342.1','AAW39262.1','AAW40361.1','AAW39214.1','AAW39256.1')
CBDB1 <- c('ADC74632.1','CAI82345.1','CAI83513.1','ACZ62529.1','CAI83526.1','ADC74550.1','ADC74673.1','ABQ16712.1','ABQ16703.1','CAI83531.1','CAI82479.1','ADC74676.1','ADC74693.1','CAI83519.1','ADC74643.1','CAI82340.1','AAW39273.1','ABQ16764.1','ADC74641.1','CAI83563.1','ADC74667.1','ADC74627.1','ABQ16695.1','ADC74658.1','ACZ62535.1','CAI83191.1','CAI83570.1','ACZ62486.1','ADC74548.1','CAI83566.1','ADC74651.1')
GT <- c('ADC74632.1','ADC74641.1','ACZ62529.1','AAW39215.1','ACZ62391.1','ABQ16764.1','ADC74693.1','ACZ62486.1','ADC74548.1','ADC74627.1','ADC74550.1','ADC74667.1','ADC74673.1','ADC74658.1','ADC74643.1','ADC74651.1','AAW39273.1','ACZ62535.1','ADC74676.1')
VS <- c('ACZ62407.1','ACZ62529.1','ACZ62474.1','ACZ62428.1','ACZ62391.1','ACZ62439.1','ACZ62362.1','ACZ61255.1','ACZ62501.1','ACZ62364.1','ACZ62492.1','ACZ62426.1','ACZ62463.1','ACZ62443.1','AAW39605.1','ACZ62520.1','ACZ62415.1','ACZ62484.1','ABQ16863.1','ACZ62482.1','ACZ62441.1','ABQ16882.1','ACZ62413.1','ACZ62535.1','ACZ62448.1','ACZ62470.1','AAW40575.1','ACZ62446.1','ACZ62423.1','ACZ62486.1','ACZ62452.1','ACZ61272.1','ACZ62477.1','ACZ62419.1','ACZ61269.1')
BAV1 <- c('ABQ16703.1','ABQ16764.1','ABQ17429.1','ABQ16710.1','ABQ16871.1','ABQ16695.1','ABQ16863.1','ABQ16882.1','ABQ16868.1','ABQ16712.1')
rdh_set <- c('ACZ62407.1','ADC74632.1','ACZ62474.1','ACZ62477.1','ADC74550.1','AAW40589.1','ACZ62529.1','AAW39273.1','ACZ62413.1','ACZ62362.1','AAW39843.1','ADC74673.1','ACZ61255.1','ACZ62501.1','CAI82345.1','AAW39262.1','ACZ62391.1','AAW39214.1','ACZ62364.1','ABQ16712.1','ACZ62428.1','AAW40468.1','ADC74548.1','AAW39053.1','ACZ62426.1','CAI82479.1','CAI83513.1','ACZ62463.1','ADC74643.1','ADC74693.1','CAI83519.1','ACZ62443.1','AAW39605.1','AAW39229.1','ACZ62520.1','CAI82340.1','ACZ62415.1','ACZ62484.1','CAI83526.1','ACZ61272.1','ACZ62448.1','ABQ16703.1','ACZ62441.1','ADC74641.1','AAW39215.1','AAW40356.1','ABQ16710.1','CAI83563.1','CAI83531.1','ABQ16871.1','ABQ16868.1','ADC74627.1','ABQ16695.1','ACZ62482.1','ACZ62439.1','AAW39060.1','ADC74658.1','AAW40373.1','ACZ62535.1','AAW39256.1','AAW39240.1','ACZ61269.1','ADC74667.1','ACZ62470.1','ABQ16764.1','AAW40575.1','CAI83191.1','ACZ62446.1','ACZ62423.1','CAI83570.1','ACZ62486.1','ACZ62452.1','ACZ62492.1','ADC74676.1','AAW40342.1','ABQ16863.1','CAI83566.1','ABQ17429.1','AAW40361.1','ADC74651.1','ACZ62419.1','ABQ16882.1')

cov = 0.34
iterations = 10000
result.matrix = matrix(nrow=iterations,ncol=5)
colnames(result.matrix) <- c("strain195","BAV1","CBDB1","VS","GT")
for(x in 1:iterations)
{
strain_table <- data.frame(matrix(0,nrow=5,ncol=6), row.names=c("strain195","BAV1","CBDB1","VS","GT"))

rdh_strain_list <- list()

rdh_strain_list[['strain195']] <- strain195
rdh_strain_list[['BAV1']] <- BAV1
rdh_strain_list[['CBDB1']] <- CBDB1
rdh_strain_list[['VS']] <- VS
rdh_strain_list[['GT']] <- GT


for(i in 1:5)
{
	for(j in 1:6)
	{
		strain_table[i,j] <- sample(500*c(1:20),1)
		
	}
}


rdh_table <- data.frame(matrix(0, nrow = 82, ncol = 6), row.names = rdh_set)

for(s in c("strain195","BAV1","CBDB1","VS","GT"))
{
for(i in 1:82)
{
	rdh_name <- row.names(rdh_table)[i]
	for(t in 1:6)
	{
		if(is.element(rdh_name, rdh_strain_list[[s]]))
			{
			rdh_table[rdh_name, t] = rdh_table[rdh_name, t] + rnorm(1, mean=strain_table[s,t], sd=cov*strain_table[s,t])
			}
	}
}
}

array_data_table <- data.frame(matrix(0, nrow = 82, ncol = 6), row.names = rdh_set)

for(t in 2:6)
{
	for(r in 1:82)
	{
		array_data_table[r, t] = log(rdh_table[r, t]/rdh_table[r, 1])
	}
}

#plot(hclust(dist(array_data_table)))
#rect.hclust(hclust(dist(array_data_table)), k=5)
model_result <- cutree(hclust(dist(array_data_table)), k=5)

reconstructed_rdh_list <- list()

for(i in 1:5)
{
	reconstructed_rdh_list[[i]] = list()
}

for(i in 1:length(model_result))
{
	accession <- names(model_result)[i]
	strain <- model_result[[i]]
	
	reconstructed_rdh_list[[strain]][length(reconstructed_rdh_list[[strain]]) + 1] <- unlist(accession)
}

strain.list <- c("simstrain1","simstrain2","simstrain3","simstrain4","simstrain5")
reconstructed_rdh_list_2 <- data.frame(row.names=strain.list)

for(i in 1:5)
{
	reconstructed_rdh_list_2[[strain.list[i]]] = list(unlist(reconstructed_rdh_list[[i]]))
}

allstrains.data <- data.frame(matrix(0, nrow = length(rdh_set), ncol=10), row.names = rdh_set)
colnames(allstrains.data) <- c("strain195","BAV1","CBDB1","VS","GT","simstrain1","simstrain2","simstrain3","simstrain4","simstrain5")

for(i in 1:dim(allstrains.data)[1]) if(is.element(row.names(allstrains.data)[i], strain195)) allstrains.data[i,'strain195'] <- 1
for(i in 1:dim(allstrains.data)[1]) if(is.element(row.names(allstrains.data)[i], BAV1)) allstrains.data[i,'BAV1'] <- 1
for(i in 1:dim(allstrains.data)[1]) if(is.element(row.names(allstrains.data)[i], CBDB1)) allstrains.data[i,'CBDB1'] <- 1
for(i in 1:dim(allstrains.data)[1]) if(is.element(row.names(allstrains.data)[i], VS)) allstrains.data[i,'VS'] <- 1
for(i in 1:dim(allstrains.data)[1]) if(is.element(row.names(allstrains.data)[i], GT)) allstrains.data[i,'GT'] <- 1

for(i in 1:dim(allstrains.data)[1]) if(is.element(row.names(allstrains.data)[i], unlist(reconstructed_rdh_list_2[['simstrain1']]))) allstrains.data[i,'simstrain1'] <- 1
for(i in 1:dim(allstrains.data)[1]) if(is.element(row.names(allstrains.data)[i], unlist(reconstructed_rdh_list_2[['simstrain2']]))) allstrains.data[i,'simstrain2'] <- 1
for(i in 1:dim(allstrains.data)[1]) if(is.element(row.names(allstrains.data)[i], unlist(reconstructed_rdh_list_2[['simstrain3']]))) allstrains.data[i,'simstrain3'] <- 1
for(i in 1:dim(allstrains.data)[1]) if(is.element(row.names(allstrains.data)[i], unlist(reconstructed_rdh_list_2[['simstrain4']]))) allstrains.data[i,'simstrain4'] <- 1
for(i in 1:dim(allstrains.data)[1]) if(is.element(row.names(allstrains.data)[i], unlist(reconstructed_rdh_list_2[['simstrain5']]))) allstrains.data[i,'simstrain5'] <- 1

require(vegan)

#plot(hclust(vegdist(t(allstrains.data), method="jaccard"), method="complete"))

strain.matrix <- as.matrix(vegdist(t(allstrains.data), method="jaccard"))

strain.pairs <- data.frame(row.names = c("strain195","BAV1","CBDB1","VS","GT","simstrain1","simstrain2","simstrain3","simstrain4","simstrain5"))

for(realstrain in c("strain195","BAV1","CBDB1","VS","GT"))
{
	distance = 1.0
	for(fakestrain in c("simstrain1","simstrain2","simstrain3","simstrain4","simstrain5"))
	{
		new.distance <- strain.matrix[realstrain,fakestrain]
		if(new.distance < distance)
		{
			distance <- new.distance
			strain.pairs[realstrain,1] <- fakestrain
		}
	}
}

matches = 0
strains.matched = vector()

for(fakestrain in c("simstrain1","simstrain2","simstrain3","simstrain4","simstrain5"))
{
	distance = 1.0
	for(realstrain in c("strain195","BAV1","CBDB1","VS","GT"))
	{
		new.distance <- strain.matrix[fakestrain,realstrain]
		if(new.distance < distance)
		{
			distance <- new.distance
			strain.pairs[fakestrain,1] <- realstrain
		}
	}
	if(fakestrain == strain.pairs[strain.pairs[fakestrain,1],1])
	{
		matches = matches +1
		strains.matched[matches] <- strain.pairs[fakestrain,1]
	}
}

for(realstrain in strains.matched)
{
	is.fraction <- length(intersect(unlist(rdh_strain_list[[realstrain]]), unlist(reconstructed_rdh_list_2[[strain.pairs[realstrain,1]]]))) / length(unlist(rdh_strain_list[[realstrain]]))
	result.matrix[x,][realstrain] <-is.fraction
}

}

result.summary = matrix(nrow=iterations,ncol=5)
colnames(result.summary) <- c("matches","mean","sd","median","mad")


for(x in 1:iterations)
{
	matches = 0
	for(y in 1:5)
	{
		if(is.na(result.matrix[x,y]) == FALSE)
		{
			matches = matches + 1
		}
	}
	result.summary[x,'sd'] <- sd(result.matrix[x,], na.rm =TRUE)
	result.summary[x,'mean'] <- mean(result.matrix[x,], na.rm =TRUE)
	result.summary[x,'mad'] <- mad(result.matrix[x,], na.rm =TRUE)
	result.summary[x,'median'] <- median(result.matrix[x,], na.rm =TRUE)
	result.summary[x,'matches'] <- matches
}