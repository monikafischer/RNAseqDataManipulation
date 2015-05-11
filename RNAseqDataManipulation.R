###### Playing with Abby's PP1 data and Deb Bell-Pederson's ADV1 data #

pp1dat<- read.table("Desktop/ADV1_RNAseq/RNAseqPP1.txt", header=TRUE)
pp1.dt <- data.table(pp1dat)
adv1dat <- read.table("Desktop/ADV1_RNAseq/RNAseqADV1.txt", header=TRUE)
adv1.dt <- data.table(adv1dat)

#read in a table with the annotations!
#read.table can't handel the spaces between words in the annotations, but read.delim can!
GeneAnnot <- read.delim("Desktop/ADV1_RNAseq/AnnotatedGenes.txt", header=TRUE)
GeneAnnot.dt <- data.table(GeneAnnot)

# Matching RNAseq data from two different dataframes based on NCU#
# The two datatables contain slightly different NCU's
# This code matches the data from each experiment to each gene
# Genes that were only represented in one experiment are given an "NA" for the other experiment.
# Add the PP1 datatable to the ADV1 datatable:
adv1pp1.dt <- adv1.dt #make a copy of adv1.dt to turn into adv1.dt + pp1.dt
adv1pp1.dt$Gene_pp1 <- pp1.dt$Gene_pp1[match(adv1pp1.dt$Gene_adv1, pp1.dt$Gene_pp1, nomatch="NA")]
adv1pp1.dt$germWT_FPKM <- pp1.dt$X5hr_WT_FPKM[ match(adv1pp1.dt$Gene_adv1, pp1.dt$Gene_pp1, nomatch="NA")]
adv1pp1.dt$germPP1_FPKM <- pp1.dt$X5hr_PP1_FPKM[ match(adv1pp1.dt$Gene_adv1, pp1.dt$Gene_pp1, nomatch="NA")]
adv1pp1.dt$germ_log2_foldchange <- pp1.dt$log2_fold_change[ match(adv1pp1.dt$Gene_adv1, pp1.dt$Gene_pp1, nomatch="NA")]
#adv1pp1.dt$Broad_Annotation <- GeneAnnot.dt$Broad.annotation[match(adv1pp1.dt$Gene_adv1, pp1.dt$Gene_pp1, nomatch="NA")]
#the above method adds columns to the end of the data.table
#but I want to add the annotations next to the first column:
adv1pp1.df <- data.frame(adv1pp1.dt)
adv1pp1.df <- as.data.frame(append(adv1pp1.df, list(GeneAnnot$Broad.annotation[match(adv1pp1.df$Gene_adv1, GeneAnnot$Gene)]), after=1))
# The above doesn't work, I initially did this by adding the annotations at the end,
# then copying that column to the 2nd position, then deleting the extra column at the end
# ...there has to be a more elegant way of doing this...
colnames(adv1pp1.df)[2] <- "Broad_Annotation" #Change the name of the second column
write.csv(adv1pp1.df, "Desktop/ADV1_RNAseq/ADV1andPP1_RNAseq.csv")


WTonly <- data.table(adv1dat$Gene_adv1, adv1dat$DD_WT_FPKM)
colnames(WTonly) <- c("Gene_adv1", "WTadv1_FPKM")
WTonly$Gene_pp1 <- pp1dat$Gene_pp1[match(WTonly$Gene_adv1, pp1dat$Gene_pp1, nomatch="NA")]
WTonly$WTpp1_FPKM2 <- pp1dat$X5hr_WT_FPKM[match(WTonly$Gene_adv1, pp1dat$Gene_pp1, nomatch="NA")]

#checking that matching actually worked:
setkey(WTonly, Gene_adv1)
setkey(adv1.dt, Gene_adv1)
setkey(pp1.dt, Gene_pp1)
WTonly["NCU08774"] #outputs all values in rows that contain "NCU08774"
adv1.dt["NCU08774"]
pp1.dt["NCU08774"]

#create data.table of only 5-fold differentially expressed genes in adv1 and pp1:
adv1pp1.dt <- data.table(adv1pp1)
adv1_5up <- subset(adv1pp1.dt, DD_log2_fold_change > 2.3)
adv1_5down <- subset(adv1pp1.dt, DD_log2_fold_change < -2.3)
pp1_5up <- subset(adv1pp1.dt, PP1_log2_foldchange < -2.3)
pp1_5down <- subset(adv1pp1.dt, PP1_log2_foldchange > 2.3)
adv1_diff <- rbind(adv1_5up, adv1_5down)
pp1_diff <- rbind(pp1_5up, pp1_5down)
setkey(adv1_diff, Gene_adv1, DD_log2_fold_change)
setkey(pp1_diff, Gene_adv1, PP1_log2_foldchange)
adv1pp1_diff <- merge(adv1_diff, pp1_diff, by.x="Gene_adv1", by.y="Gene_adv1")

####compare differentially expressed genes between adv1 and pp1 experiments###
#creating the pp1 differentially expressed dataset to compare with adv1
pp1.dt <- data.table(pp1dat)
pp15up <- subset(pp1.dt, log2_fold_change < -2.3)
pp15down <- subset(pp1.dt, log2_fold_change > 2.3)
pp1diff <- rbind(pp15up, pp15down)
pp1diff2 <- pp1diff[ X5hr_PP1_FPKM > 0 | X5hr_WT_FPKM > 0] #removes rows where both columns = 0
pp1diff3 <- pp1diff[ X5hr_PP1_FPKM > 0] #removes all rows where PP1_FPKM = 0
pp1diff4 <- pp1diff3[X5hr_WT_FPKM > 0] #removes all rows where WT_FPKM = 0
#pp1diff4 = data.table of all genes with 5fold change up or down greater than 0. (256 genes total)

adv1.dt <-data.table(adv1dat)
#create new datatable that only has significantly differentially expressed genes:
adv1.dt.p <- subset(adv1.dt, DD_p_value < 0.001) #139 genes

#pp1diff4 and adv1.dt.p are the datasets I want to compare
common <- data.table(intersect(pp1diff4$Gene_pp1, adv1.dt.p$Gene_adv1)) #20
pp1only <- data.table(setdiff(pp1diff4$Gene_pp1, adv1.dt.p$Gene_adv1)) #236
adv1only <- data.table(setdiff(adv1.dt.p$Gene_adv1, pp1diff4$Gene_pp1)) #119
#adv1.dt.p = 139 = 119 + 20
#pp1diff4 = 256 = 236 + 20
write.csv(common, "Desktop/commongenes.csv")

#Simple black & white Venn Diagram:
venn(list(first.vector = pp1diff4$Gene_pp1, adv1.dt.p$Gene_adv1))

#pretty venn diagram, but for some reason I can't edit it in in illustrator :( 
venn.diagram(list(pp1=pp1diff4$Gene_pp1,  adv1=adv1.dt.p$Gene_adv1), 
             filename = "Desktop/venn2.eps", #venndiagram doesn't show up in plots window...
             fill = c("green", "purple"), #circle fill color
             col = c("transparent", "transparent"), #circle border color ("transparent" = no boarder)
             cex = 2.5, #font size for the numbers inside the circles
             cat.cex = 2, #font size for the circle labels
             rotation.degree = 0, #90 or 270 = circles on top of eachother, 0 or 180 = circles next to eachother
             main = "Differential expressed genes in adv1 and pp1",
             sub = "this is a subtitle",
             main.cex = 1.5,
             sub.cex = 1
             );


####trying to figure out how to generate p-values for the pp1 data###
#Abby and Wilfried's paper says they did a Fisher's Exact Test...
#First make a new dataframe with just the gene name and FPKM's
getwd() #get working directory
pp1dat<- read.table("Desktop/ADV1_RNAseq/RNAseqPP1.txt", header=TRUE)
pp1FPKM <- data.frame(pp1dat$Gene_pp1, pp1dat$X5hr_PP1_FPKM, pp1dat$X5hr_WT_FPKM)
colnames(pp1FPKM) <- c('gene', 'pp1.FPKM', 'WT.FPKM')

fisher <- fisher.test(pp1FPKM)
#Basic code for one fisher's exact test, but how to do it on each row?...
#A fisher's exact test needs 4 datapoints to build a 2x2 table
#Abby's data only has two data points for each gene...
#How did they do a Fisher's Exact Test?!

fisher.test(pp1FPKM[,c(2,3)])


