library(data.table)
library(stringr)
library(Biostrings)
library(ggplot2)
library(plyr)
library(MASS)


option_list = list(
  make_option(c("-n", "--no.sample"), type="double", default=NULL, help="sample's no", metavar="character"),
  make_option(c("-w", "--second.allele"), type="character", default=NULL, help="whether to select the second allele, Y or N", metavar="character"),
  make_option(c("-m", "--methods"), type="character", default=NULL, help="v-voting, m-max, k-Kendall", metavar="character")
) 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


## parameters
sample.no = opt$no.sample  #
second.allele = opt$second.allele  #
method.lag = opt$methods  #



#######################to combine multiple files with different three different exposure times##################################################
#######################and to produce a file ###########################

filess = "D:\\Jiang_work\\Projects\\COVID19\\data"
setwd(filess)


getdata = function(data){
  
  data_all = data.table()
  
  
  locations = unique(data$location)
  locations = locations[!is.na(locations)]
  for (ii in min(locations):max(locations)){
    sub_d = subset(data,location == ii)
    if (nrow(sub_d)== 16){
      sub_d = subset(sub_d,gene != "")
    }
    
    
    y_idx = unique(sub_d$y)
    
    for (jj in y_idx){
      sub_d_1 = subset(sub_d,y == jj)
      
      
      if (nrow(sub_d_1) == 4){
        sub_d_1 = sub_d_1[order(sub_d_1$base),]
        
        aa = data.frame(ID = jj, location = ii, strand = unique(sub_d_1$strand),A = sub_d_1$mean[1],
                        C = sub_d_1$mean[2], G = sub_d_1$mean[3], T = sub_d_1$mean[4]
                        ,stringsAsFactors = FALSE)
        
        data_all = rbind(data_all,aa)
        
        aa = match(sub_d$y,jj)
        sub_d = sub_d[is.na(aa),]
      }
      
      
    }
    
    if (nrow(sub_d) == 4){
      
      sub_d = sub_d[order(sub_d$base),]
      
      aa = data.frame(ID = jj, location = ii, strand = unique(sub_d$strand),A = sub_d$mean[1],
                      C = sub_d$mean[2], G = sub_d$mean[3], T = sub_d$mean[4]
                      ,stringsAsFactors = FALSE)
      
      data_all = rbind(data_all,aa)
      
    }
    
    
  }
  
  return(data_all)
  
}


d_anno = fread("core0_annot_mod.tsv",header = T,sep = "\t")
d_anno = setDF(d_anno)
d_anno$location = d_anno$start + 12
d_anno$strand = str_sub(d_anno$probe_seq_5up,13,13)
d_anno$strand[d_anno$strand != d_anno$base] = "antisense"
d_anno$strand[d_anno$strand == d_anno$base] = "sense"
d_anno = subset(d_anno,species == "SARS-COV2")

d_anno = d_anno[,c(1,2,5,9,18,19)]
d_anno$ID = paste(d_anno$x,d_anno$y,sep = "-")
d_anno = d_anno[,c(3:7)]




#######################annotation file ############################################


exposes = c("0.5s","1s","4s")


for (expose in exposes){
  
  setwd(paste0(filess,"/",expose))
  
####################### intensty file 0.5s ############################################
  
  filename = list.files()
  filename = filename[grepl(paste0(No.sample,"_"),filename)]
  d_1 = fread(filename,header = T,sep = ",")
  d_1$y = 495 - d_1$y
  d_1$ID = paste(d_1$x,d_1$y,sep = "-")
  d_1 = d_1[,c(2,3,5,14)]
  
####################### merge annation and intensty files############################################
  
  
  d_1 = merge(d_1,d_anno,by.x = "ID",by.y = "ID")
  data_all = getdata(d_1)
  
  
  
  data_all$max = apply(data_all[,c(4:7)], 1, max)
  data_all$idx = apply(data_all[,c(4:7)], 1, which.max)
  
  
  data_all$min = apply(data_all[,c(4:7)], 1, min,na.rm=TRUE)
  
  for (ii in 1:nrow(data_all)){
    data_all[ii,(data_all$idx[ii]+3)] = NA
  }
  
  
  bases = c("A","C","G","T")
  data_all$idx = bases[data_all$idx]
  data_all = setDF(data_all)
  
  
  
  
  data_all$A = (data_all$A )/data_all$max
  data_all$C = (data_all$C )/data_all$max
  data_all$G = (data_all$G )/data_all$max
  data_all$T = (data_all$T )/data_all$max
  
  
  data_all$max_2 = apply(data_all[,c(4:7)], 1, max,na.rm=TRUE)
  data_all$idx_2 = apply(data_all[,c(4:7)], 1, which.max)
  bases = c("A","C","G","T")
  data_all$idx_2 = bases[data_all$idx_2]
  
  
  
  
  
  
  data_all$DIS = data_all$min/data_all$max
  
  setwd(paste0(filess,"/data_all"))
  write.csv(data_all,paste(No.sample,expose,"all.csv",sep = "_"),row.names = FALSE)
  
  
  
  data_all = subset(data_all,DIS <=0.5)
  data_all = setDF(data_all)
  
  
  data_p = data.frame()
  
  for (ii in 4:7){
    aa = data.frame(label = data_all[,3],value = data_all[,ii],stringsAsFactors = FALSE)
    data_p = rbind(data_p,aa)
  }
  
  
  data_p = data_p[!is.na(data_p$value),]
  f1n <- fitdistr(data_p$value,"normal")
  
  
  
  data_all$p = pnorm(data_all$max_2, mean = f1n$estimate[1], sd = f1n$estimate[2] , lower.tail = FALSE, log.p = FALSE)
  
  
  write.csv(data_all,paste(No.sample,expose,"cal_p.csv",sep = "_"),row.names = FALSE)
  
  
}






#######################to calculate the base with maximum value #################################



setwd("D:\\Jiang_work\\Projects\\COVID19\\data\\data_all")


f1 <- function(x) { d <- !duplicated(x) ; data.frame(uniqueValue=x[d], firstIndex=which(d)) }


getmatrix = function(d1,colname){
  aa = matrix(NA,1,6)
  colnames(aa) = colname
  aa = data.frame(aa)
  l = match(d1$ID,colname)
  aa[1,l] = paste0(d1$idx,"(",d1$value,")")
  
  return(aa)
}



filenames = list.files()
filenames = filenames[grepl("_all.csv",filenames)&grepl(sample.no,filenames)]

data_all = data.frame()
for (ii in filenames){
  d1 = fread(ii,header = T,sep = ",")
  d1 = setDF(d1)
  d1$se = ii
  data_all = rbind(data_all,d1)
}


data_all$se = gsub("s_all.csv","",data_all$se)

d1_ref = fread("ref.csv",header = T,sep = ",")
d1_ref = setDF(d1_ref)


data_all = merge(data_all,d1_ref,by.x = "location",by.y = "location")


data_all$A = data_all$A*data_all$max
data_all$C = data_all$C*data_all$max
data_all$G = data_all$G*data_all$max
data_all$T = data_all$T*data_all$max
data_all$max_2 = data_all$max_2*data_all$max



data_all$value = data_all$max/(rowMeans(data_all[,4:7],na.rm = TRUE))


data_all$strand = gsub("antisense","AS",data_all$strand)
data_all$strand = gsub("sense","S",data_all$strand)

data_all$ID = paste0(data_all$se,data_all$strand)


locations = unique(data_all$location)
locations = sort(locations)


colname = unique(data_all$ID)
data_end = data.frame()




for (ii in locations){
  subdata = subset(data_all,location == ii)
  
  un_ID = unique(subdata$ID)
  if (length(un_ID) < nrow(subdata)){
    dd = f1(subdata$ID) 
    d1 = subdata[dd$firstIndex,]
    aa = getmatrix(d1,colname)
    aa$location = ii
    data_end = rbind(data_end,aa)
    subdata = subdata[-dd$firstIndex,]
  }
  
  aa = getmatrix(subdata,colname)
  aa$location = ii
  data_end = rbind(data_end,aa)
  
}


data_end = merge(data_end,d1_ref,by.x = "location",by.y = "location")
colnames(data_end)[8] = "ref"


for (ii in 2:7){
  
  bb = !is.na(data_end[,ii])
  aa = sapply(strsplit(data_end[bb,ii], split='(', fixed=TRUE), function(x) (x[]))
  data_end$base[bb] = aa[1,]
  data_end[bb,ii] = aa[2,]
  data_end[,ii] = as.numeric(gsub(")","",data_end[,ii]))
  colnames(data_end)[ncol(data_end)] = paste0(colnames(data_end)[ii],"base")
}

aa = sapply(strsplit(colnames(data_end)[2], split='_', fixed=TRUE), function(x) (x[1]))

colnames(data_end) = gsub(aa,"X",colnames(data_end))
colnames(data_end) = gsub("_","",colnames(data_end))


data_end = data_end[,c("location","X0.5S","X0.5Sbase","X0.5AS","X0.5ASbase",
                       "X1S","X1Sbase","X1AS","X1ASbase",
                       "X4S", "X4Sbase", "X4AS","X4ASbase" ,"ref" )]



write.csv(data_end,paste0(sample.no,".csv"),row.names = F)



#######################The weighted voting method #################################





data_end = fread(paste0(sample.no,".csv"),header = T,sep = ",")
data_end = setDF(data_end)


n = which(data_end$X0.5Sbase == data_end$X0.5ASbase &
            data_end$X0.5Sbase == data_end$X1Sbase &
            data_end$X0.5Sbase == data_end$X1ASbase &
            data_end$X0.5Sbase == data_end$X4Sbase &
            data_end$X0.5Sbase == data_end$X4ASbase & !is.na(data_end$X0.5Sbase))

d1 = data_end[-n,]

weights = data.frame(X0.5S = sum(d1$X0.5Sbase == d1$ref & !is.na(d1$X0.5Sbase), na.rm = TRUE)/sum(!is.na(d1$X0.5Sbase)),
                     X0.5AS = sum(d1$X0.5ASbase == d1$ref & !is.na(d1$X0.5ASbase), na.rm = TRUE)/sum(!is.na(d1$X0.5ASbase)),
                     X1S = sum(d1$X1Sbase == d1$ref & !is.na(d1$X1Sbase), na.rm = TRUE)/sum(!is.na(d1$X1Sbase)),
                     X1AS = sum(d1$X1ASbase == d1$ref & !is.na(d1$X1ASbase), na.rm = TRUE)/sum(!is.na(d1$X1ASbase)),
                     X4S = sum(d1$X4Sbase == d1$ref & !is.na(d1$X4Sbase), na.rm = TRUE)/sum(!is.na(d1$X4Sbase)),
                     X4AS = sum(d1$X4ASbase == d1$ref & !is.na(d1$X4ASbase), na.rm = TRUE)/sum(!is.na(d1$X4ASbase)),
                     stringsAsFactors = FALSE
)

weights = weights/sum(weights)


data_end = data_end[,-14]

data_all = data.frame()

location = unique(data_end$location)
for (ii in location){
  sub_data = subset(data_end,location == ii)
  bases = unique(unlist(sub_data[,c(3,5,7,9,11,13)]))
  bases = bases[!is.na(bases)]
  matrix_all = matrix(0,1,length(bases))
  colnames(matrix_all) = bases
  
  for (jj in 1:nrow(sub_data)){
    for (kk in c(3,5,7,9,11,13)){
      if (!is.na(sub_data[jj,kk])){
        idx =(match(sub_data[jj,kk],bases))
        value = sub_data[jj,(kk-1)]
        weight = weights[colnames(sub_data)[kk-1]]
        matrix_all[1,idx] = matrix_all[1,idx] + as.numeric(value*weight) 
      }
    }
  }
  base_r = colnames(matrix_all)[which.max(matrix_all)]
  
  
  sub_data = subset(data_end,location == ii)
  sub_data_1 = sub_data[,c(2,4,6,8,10,12)]
  aa = max(apply(sub_data_1, 1, max, na.rm = TRUE))
  
  for (jj in 1:nrow(sub_data)){
    idx = (sub_data[jj,] == aa)
    if (sum(idx,na.rm = TRUE)>0){
      base_m = sub_data[jj,(which(idx)+1)]
    }
  }
  
  data_all = rbind(data_all,data.frame(location = ii, voting = base_r, max = base_m,stringsAsFactors = FALSE))
  
}


write.csv(data_all,paste0("result_",sample.no,".csv"),row.names = F)




##################################The credibility score method ####################################################


library(data.table)

setwd("D:\\Jiang_work\\Projects\\COVID19\\data\\data_all")
sample.no = "30"

filenames = list.files()
filenames = filenames[grepl("_all.csv",filenames)&grepl(sample.no,filenames)]

data_all = data.frame()
for (ii in filenames){
  d1 = fread(ii,header = T,sep = ",")
  d1 = setDF(d1)
  d1$se = ii
  data_all = rbind(data_all,d1)
}

data_all$se = gsub("30_","",data_all$se)
data_all$se = gsub("_all.csv","",data_all$se)




d1_ref = fread("ref.csv",header = T,sep = ",")
d1_ref = setDF(d1_ref)




data_all = merge(data_all,d1_ref,by.x = "location",by.y = "location")


data_all$isref = TRUE
data_all$isref[data_all$idx!=data_all$base] = FALSE

data_all$A = data_all$A*data_all$max
data_all$C = data_all$C*data_all$max
data_all$G = data_all$G*data_all$max
data_all$T = data_all$T*data_all$max
data_all$max_2 = data_all$max_2*data_all$max


data_all$D = data_all$max - data_all$min
data_all$drel = (data_all$max - data_all$max_2)/data_all$D



data_correct = subset(data_all,isref == TRUE)
data_incorrect = subset(data_all,isref == FALSE)

data_all$score = NA

for (ii in 1:nrow(data_all)){
   D1 = data_all$D[ii]
  drel1 = data_all$drel[ii]
  
  Wc = nrow(subset(data_correct,D <= D1 & drel <= drel))
  Wic = nrow(subset(data_incorrect,D >= D1 & drel >= drel))
  
  data_all$score[ii] = Wc/(Wc+Wic+1)
  
  
  
}


locations = as.numeric(unique(data_all$location)) 

seq_result = data.frame()
for (ii in locations){
  sub_data = subset(data_all,location == ii)
  sub_data = sub_data[which.max(sub_data$score),]
  seq_result = rbind(seq_result,data.frame(location = ii, base = sub_data$idx,stringsAsFactors = FALSE))
}


write.csv(seq_result,"seq.csv",row.names = FALSE)















