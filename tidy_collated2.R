file = "/export/share/prkfs2/shared/bioinf-data/Papenfuss_lab/projects/lachlan/PHD/GD_matrices/combined.csv"
file = "/export/share/prkfs2/shared/bioinf-data/Papenfuss_lab/projects/lachlan/PHD/GD_matrices3/precomputed/collated_c32_p128_v1.csv"
file = "/vast/scratch/users/lmcintosh/GD2/GD/collated_128_p128_v3.csv"
# install.packages("remotes")
# remotes::install_github("traversc/trqwe")
# install.packages("data.table")
library("data.table")
setDTthreads(threads = 0)
data = fread(file)
saveRDS(data, "/vast/scratch/users/lmcintosh/GD2/GD/collated_128_p128_v3_uncompressed.Rdata",compress=F)
#data = readRDS("/vast/scratch/users/lmcintosh/GD2/GD/collated_32_p128_v3.Rdata")

data = data[order(data[,2]),]
data = data[order(data[,1]),]
head(data)
dim(data)
path_code = data[,3]
paths = sapply(path_code,function(x) strsplit(x,"G"))
nWGD = sapply(paths, function(x) length(x)-1)
pre = sapply(paths, function(x) as.numeric(x[[1]]))
mid = sapply(paths, function(x){
  if(length(x) >= 2){
    as.numeric(x[[2]])
  } else{
    NA
  }
})
post = sapply(paths, function(x){
  if(length(x) == 3){
    as.numeric(x[[3]])
  } else{
    NA
  }
})
tWGD1 = sapply(paths, function(x){
  if(length(x)>1){
    as.numeric(x[[1]])+1
  } else{
    NA
  }
})
tWGD2 = sapply(paths, function(x){
  if(length(x)>2){
    as.numeric(x[[1]])+as.numeric(x[[2]])+2
  } else{
    NA
  }
})
pup = data[,1]/100
pdown = data[,2]/100
updown = data[,1:2]/100
colnames(updown) <- c("pup","pdown")



CNS = data[,4:ncol(data)]
CNS = data.frame(lapply(CNS,as.numeric))
dim(CNS)
for(col in 1:dim(CNS)[2]){
  print("#")
  print(col)
  those = which(CNS[,col] < 0)
  print(length(those))
  CNS[those,col] = rep(0,length(those))
  those = which(CNS[,col] > 1)
  print(length(those))
  CNS[those,col] = rep(1,length(those))
}

CNS = data.frame(t(apply(CNS,1,function(x) x/sum(x))))
#alldata[which(!complete.cases(CNS)),]

colnames(CNS) <- sapply(0:(ncol(CNS)-1),function(x) as.character(x))
colnames(pup) <- "pup"
colnames(pdown) <- "pdown"
alldata = cbind(pup,pdown,nWGD,pre,mid,post,tWGD1,tWGD2,CNS)
alldata$paths = path_code
head(alldata)

saveRDS(alldata, "/vast/scratch/users/lmcintosh/GD2/GD/collated_128_p128_v3_final_uncompressed.Rdata",compress=F)

# 
# alldata[which(alldata[,9] == 1 & rowSums(alldata[,10:41]) != 0),]
# alldata$zero = alldata[,9]
# alldata[which(alldata[,9] == 1),"zero"] = 1-rowSums(alldata[which(alldata[,9] == 1),10:41])
# alldata[which(alldata$zero != alldata[,9]),"zero"]
# # up to here! fix this!
# 
# parts = strsplit(file, "/")
# name = parts[[1]][length(parts[[1]])]
# name = paste0("WGD_",name)
# #
# # write.csv(alldata,name)
# saveRDS(alldata, file = paste0(name,".Rda"),compress = FALSE)
# #
# # # should save this for later things
# # library("data.table")
# # alldata = fread(name)
# 
# 
# alldata = readRDS(file = paste0(name,".Rda"))
# 
# dim(alldata)
# # rs = rowSums(alldata[,10:(10+31)])
# # alldata[which(abs(rs- 1)>1)]
# 
# get_CN_counts <- function(CNS,threshold= 0.01){
#   major = unlist(CNS[,1])
#   minor = unlist(CNS[,2])
#   count = unlist(CNS[,3])
#   
#   major[which(CNS[,2] > CNS[,1])] = CNS[which(CNS[,2] > CNS[,1]),2]
#   major[which(CNS[,2] < CNS[,1])] = CNS[which(CNS[,2] < CNS[,1]),1]
#   minor[which(CNS[,2] < CNS[,1])] = CNS[which(CNS[,2] < CNS[,1]),2]
#   minor[which(CNS[,2] > CNS[,1])] = CNS[which(CNS[,2] > CNS[,1]),1]
#   
#   FAs = unlist(major[which(minor == 0)])
#   FACounts = unlist(count[which(minor == 0)])
#   
#   FA = t(sapply(unique(FAs), function(x) c(x,sum(FACounts[which(FAs == x)]))))
#   
#   nonFAs = unlist(c(minor,major[which(minor != 0)]))
#   nonFACounts = unlist(c(count,count[which(minor != 0)]))
#   
#   nonFA = t(sapply(unique(nonFAs), function(x) c(x,sum(nonFACounts[which(nonFAs == x)]))))
#   
#   if(ncol(FA) > 0 & ncol(nonFA) > 0){
#     mat = rbind(
#       cbind(FA,rep(1,nrow(FA))),
#       cbind(nonFA,rep(0,nrow(nonFA)))
#     )    
#   } else if(ncol(FA) > 0){
#     mat = cbind(FA,rep(1,nrow(FA)))
#   } else{
#     mat = cbind(nonFA,rep(1,nrow(nonFA)))
#   }
#   
#   
#   mat = mat[order(mat[,2]),]
#   mat = mat[which(mat[,2] > threshold),]
#   
#   mat = as.data.frame(mat)
#   colnames(mat) <- c("copy","freq","FA")
#   return(mat)
# }
# 
# 
# (CNS_example = matrix(c(1,1,2,0,2,2,1,5,4,2,4,2,4,2,4,1,1,2,0,2,2,1,5,4,2,4,2,4,2,4,rep(1,16)),nrow=2,byrow=F))
# (CNS_example <- matrix(c(rep(2,30),rep(4,6),rep(3,10)),nrow=2,byrow=T))
# (CNS_example <- matrix(c(rep(0,20),rep(2,20),rep(4,6)),nrow=2,byrow=T)) # GD happens with greater bayesian prob, doesn't happen ML VERY INTERESTING CASE
# 
# CNS_example  = t(rbind(CNS_example,rep(1,23)))
# (CN_table = get_CN_counts(CNS_example))
# 
# for(j in 1:nrow(CN_table)){
#   col = as.character(CN_table[j,"copy"])
#   print(j)
#   print(col)
#   if(CN_table[j,"FA"]){
#     ll = (log(alldata[,col]) - log(1-alldata[,"zero"]))*CN_table[j,"freq"]
#     ll[which(ll==Inf)] = -Inf
#     ll[which(sapply(ll,is.nan))] = -Inf
#   } else{
#     ll = log(alldata[,col])*CN_table[j,"freq"]
#   }
#   print(ll[order(-ll)])
#   print(alldata[which(ll==Inf),])
#   if(j==1){
#     loglik = ll
#   } else{
#     loglik = ll+loglik
#   }
# }
# 
# #temp = alldata[!complete.cases(alldata),]
# #dim(temp)
# #head(temp)
# alldata$loglik = loglik
# 
# alldata[which(is.na(alldata$loglik)),]
# alldata[which(alldata$pup < 0 | alldata$pdown < 0),]
# head(alldata)
# #alldata[is.na(alldata$loglik) & alldata$pup >0 & alldata$pdown >0,]
# 
# #alldata[which(alldata$loglik == -Inf),"loglik"] = NA
# gooddata = alldata[which(!is.na(alldata$loglik)),]
# dim(gooddata)
# gooddata = gooddata[order(-gooddata$loglik),]
# head(gooddata)
# 
# # CN_table
# # tot = sum(exp(gooddata$loglik),na.rm=T)
# # gooddata$dens = exp(gooddata$loglik)/tot
# 
# # gooddata = gooddata[order(-gooddata$dens),]
# # head(gooddata)
# # alldata[which(alldata$loglik == max(alldata$loglik,na.rm=T)),]
# # alldata[which(alldata$loglik == min(alldata$loglik,na.rm=T)),]
# alldata$count = 1
# alldata$prob = exp(alldata$loglik)
# 
# loglik_marginal = aggregate(alldata[,c("prob","count")],by=list(path = alldata$paths),function(x) sum(x,na.rm=TRUE))
# 
# loglik_marginal = loglik_marginal[order(-loglik_marginal$prob),]
# head(loglik_marginal,50)
# # calculate the posterior expection of the parameters 
# sum(alldata$dens * alldata$pup)
# sum(alldata$dens * alldata$pdown)
# sum(alldata$dens * (1-alldata$pup-alldata$pdown))
# 
# # these 1d plots of probabilits are useless without considering epoch paths
# #ggplot(alldata, aes(pup)) + geom_density(aes(weight=dens))
# #ggplot(alldata, aes(pdown)) + geom_density(aes(weight=dens))
# 
# 
# # geom tile isn't working great, let's do a 2d density plot with the rep function
# alldata$perc <- round(1000*alldata$dens/max(alldata$dens))
# alldata$dens2 <- alldata$dens#/1000
# 
# library(data.table)
# rep.with.dt<-data.table(alldata[which(alldata$nWGD==0),])[,list(pup=rep(pup,perc),
#                                        pdown=rep(pdown,perc),
#                                        pre=rep(pre,perc),
#                                        mid=rep(mid,perc),
#                                        post=rep(post,perc),
#                                        dens2=rep(dens2,perc))]
# ggplot(rep.with.dt, aes(pup, pdown)) + 
#   geom_density_2d_filled() +
#   facet_grid(vars(pre))
# 
# rep.with.dt<-data.table(alldata[which(alldata$nWGD==1),])[,list(pup=rep(pup,perc),
#                                                                 pdown=rep(pdown,perc),
#                                                                 pre=rep(pre,perc),
#                                                                 mid=rep(mid,perc),
#                                                                 post=rep(post,perc))]
# ggplot(rep.with.dt, aes(pup, pdown)) + 
#   geom_density_2d_filled() +
#   #facet_wrap(vars(pre,mid)) +
#   facet_grid(vars(pre),vars(mid))
# 
# rep.with.dt<-data.table(alldata[which(alldata$nWGD==2),])[,list(pup=rep(pup,perc),
#                                                                 pdown=rep(pdown,perc),
#                                                                 pre=rep(pre,perc),
#                                                                 mid=rep(mid,perc),
#                                                                 post=rep(post,perc))]
# 
# ggplot(rep.with.dt, aes(pup, pdown)) + 
#   geom_density_2d_filled() +
#   facet_wrap(vars(pre,mid,post))
# 
# repped_data0 <- data.table(alldata[which(alldata$nWGD==0),])[,
#                                                              list(pup=rep(pup,perc),
#                                                                   pdown=rep(pdown,perc),
#                                                                   pre=rep(pre,perc),
#                                                                   mid=rep(mid,perc),
#                                                                   post=rep(post,perc))]
# 
# g0 <- ggplot(repped_data0, aes(pup, pdown)) + 
#   geom_density_2d() +
#   facet_grid(vars(pre),vars(mid))
# g0
# 
# repped_data1 <- data.table(alldata[which(alldata$nWGD==1),])[,
#                                                              list(pup=rep(pup,perc),
#                                                                   pdown=rep(pdown,perc),
#                                                                   pre=rep(pre,perc),
#                                                                   mid=rep(mid,perc),
#                                                                   post=rep(post,perc))]
# 
# g1 <- ggplot(repped_data1, aes(pup, pdown)) + 
#   stat_density_2d(geom = "polygon", aes(alpha = ..level.., fill = Label))+
#   facet_grid(vars(pre),vars(mid))
# g1
# 
# repped_data <- data.table(alldata)[,list(pup=rep(pup,perc),
#                                           pdown=rep(pdown,perc),
#                                           nWGD=rep(nWGD,perc),
#                                           pre=rep(pre,perc),
#                                           mid=rep(mid,perc),
#                                           post=rep(post,perc))]
# 
# ggplot(repped_data, aes(pup, pdown)) + 
#   stat_density_2d(geom = "polygon", aes(fill = as.factor(nWGD)))+
#   facet_grid(vars(pre),vars(mid))
# 
# 
# temp <- unique(alldata[,c("pre","mid","post")])
# temp <- temp[order(temp$post),]
# temp <- temp[order(temp$mid),]
# temp <- temp[order(temp$pre),]
# 
# # why are the levels going up with dimension when that isn't even true
# 
# # maybe i need to check for duplicate rows...
# temp <- temp[,c("pup","pdown","nWGD")]
# 
# ggplot() + geom_density(aes(temp$dens))
# 
# sum(exp(alldata$loglik)/tot * alldata$nWGD)
# 
# # posterior probability of number of WGD rounds. 
# aggregate(exp(alldata$loglik)/tot, by=list(nWGD=alldata$nWGD), FUN=sum)
# 
# 
# alldata$density <- exp(alldata$loglik)/tot
# 
# seq(-4.60517019,0,0.1/4.60517019)
# 
# 10^(seq(-4,0,0.05))
# 
# get_posterior_plot_2d <- function(alldata){
#   dd <- alldata[which(alldata$nWD == 1),]
#   library(cowplot)
#   cols <- rev(rainbow(7)[-7])
#   low = min(dd$density)
#   high = max(dd$density)
#   g1 <- ggplot(dd, aes(pdown, pup, fill= density)) + 
#     geom_tile() +
#     scale_fill_gradientn(colours = cols,limits=c(low,high))+ facet_grid(vars(pre), vars(post))
#   
#   return(plot_grid(g1))
# }
# 
# get_posterior_plot_2d(alldata)
# 
# 
# 
# 
# 
# 
# 
# # lets see which of these plotting functions we can match up again
# 
# ggplot(alldata)+geom_line(aes(event_prob,density))+ facet_grid(vars(N), vars(M)) + ylim(0,ylimit)
# 
# 
# get_posterior_plot <- function(dd){
#   ylimit = max(dd$density)
#   dd$event_prob <- dd$a
#   library(cowplot)
#   g1 <- ggplot(dd[which(dd$O == -1),])+geom_line(aes(event_prob,density))+ facet_grid(vars(N), vars(M)) + ylim(0,ylimit)
#   g2 <- ggplot(dd[which(dd$O == 0),])+geom_line(aes(event_prob,density))+ facet_grid(vars(N), vars(M)) + ylim(0,ylimit)
#   g3 <- ggplot(dd[which(dd$O == 1),])+geom_line(aes(event_prob,density))+ facet_grid(vars(N), vars(M)) + ylim(0,ylimit)
#   g4 <- ggplot(dd[which(dd$O == 2),])+geom_line(aes(event_prob,density))+ facet_grid(vars(N), vars(M)) + ylim(0,ylimit)
#   return(plot_grid(g1,g2,g3,g4))
# }
# 
# get_log_posterior_plot <- function(dd){
#   #dd2 <- dd[which(dd$res>-Inf),]
#   dd2 <- dd
#   dd2$event_prob <- dd2$a
#   ymax = max(dd2$res)
#   ymin = 2*ymax 
#   library(cowplot)
#   dd2$log_density <- dd2$res
#   g1 <- ggplot(dd2[which(dd2$O == -1),])+geom_line(aes(event_prob,log_density))+ facet_grid(vars(N), vars(M)) + ylim(ymin,ymax)
#   g2 <- ggplot(dd2[which(dd2$O == 0),])+geom_line(aes(event_prob,log_density))+ facet_grid(vars(N), vars(M)) + ylim(ymin,ymax)
#   g3 <- ggplot(dd2[which(dd2$O == 1),])+geom_line(aes(event_prob,log_density))+ facet_grid(vars(N), vars(M)) + ylim(ymin,ymax)
#   g4 <- ggplot(dd2[which(dd2$O == 2),])+geom_line(aes(event_prob,log_density))+ facet_grid(vars(N), vars(M)) + ylim(ymin,ymax)
#   return(plot_grid(g1,g2,g3,g4))
# }
# 
# get_posterior_marginal_timescape <- function(dd){
#   marginal = ddply(dd, c("N", "M","O"),numcolwise(sum))
#   marginal$density = marginal$density/sum(marginal$density)
#   return(marginal[which(marginal$density > 0),])
# }
# 
# get_posterior_marginal_timescape <- function(dd){
#   marginal = ddply(dd, c("N", "M","O"),numcolwise(sum))
#   marginal$density = marginal$density/sum(marginal$density)
#   marginal = marginal[which(marginal$density > 0),]
#   marginal = marginal[order(marginal$density),]
#   return(marginal)
# }
# 
# get_posterior_marginal_timescape_N <- function(dd){
#   marginal = ddply(dd, c("N", "M","O"),numcolwise(sum))
#   marginal$density = marginal$density/sum(marginal$density)
#   marginal = marginal[which(marginal$density > 0),]
#   marginal = marginal[order(marginal$density),]
#   return(marginal)
# }
# 
# get_posterior_marginal_timescape_plot <- function(dd){
#   marginal = get_posterior_marginal_timescape(dd)
#   ylimit = max(marginal$density)*1.01
#   library(cowplot)
#   return(ggplot(marginal)+geom_bar(aes(x=M,y=density),stat='identity')+ facet_grid(vars(N),vars(O)) + ylim(0,ylimit))
# }
# 
# get_posterior_marginal_GD_prob <- function(dd){
#   GD = sum(dd[which(dd$M>-1),"density"])
#   GD2 = sum(dd[which(dd$M>-1 & dd$O>-1),"density"])
#   noGD = sum(dd[which(dd$M==-1),"density"])
#   probability_GD = GD/(GD+noGD)
#   probs012 = c(noGD,GD-GD2,GD2)/(noGD+GD)
#   return(probs012)
# }
# get_CN_counts <- function(CNS,threshold= 0.01){
#   major = unlist(CNS[,1])
#   minor = unlist(CNS[,2])
#   count = unlist(CNS[,3])
#   
#   major[which(CNS[,2] > CNS[,1])] = CNS[which(CNS[,2] > CNS[,1]),2]
#   major[which(CNS[,2] < CNS[,1])] = CNS[which(CNS[,2] < CNS[,1]),1]
#   minor[which(CNS[,2] < CNS[,1])] = CNS[which(CNS[,2] < CNS[,1]),2]
#   minor[which(CNS[,2] > CNS[,1])] = CNS[which(CNS[,2] > CNS[,1]),1]
#   
#   FAs = unlist(major[which(minor == 0)])
#   FACounts = unlist(count[which(minor == 0)])
#   
#   FA = t(sapply(unique(FAs), function(x) c(x,sum(FACounts[which(FAs == x)]))))
#   
#   nonFAs = unlist(c(minor,major[which(minor != 0)]))
#   nonFACounts = unlist(c(count,count[which(minor != 0)]))
#   
#   nonFA = t(sapply(unique(nonFAs), function(x) c(x,sum(nonFACounts[which(nonFAs == x)]))))
#   
#   if(ncol(FA) > 0 & ncol(nonFA) > 0){
#     mat = rbind(
#       cbind(FA,rep(1,nrow(FA))),
#       cbind(nonFA,rep(0,nrow(nonFA)))
#     )    
#   } else if(ncol(FA) > 0){
#     mat = cbind(FA,rep(1,nrow(FA)))
#   } else{
#     mat = cbind(nonFA,rep(1,nrow(nonFA)))
#   }
#   
#   
#   mat = mat[order(mat[,2]),]
#   mat = mat[which(mat[,2] > threshold),]
#   
#   mat = as.data.frame(mat)
#   colnames(mat) <- c("copy","freq","FA")
#   return(mat)
# }
# 
# 
# get_posterior_GD_timing_distribution <- function(dd){
#   domain = 1:1000
#   output = rep(0,length(domain))
#   marginal = ddply(dd, c("N", "M","O"),numcolwise(sum))
#   alpha = 1
#   for( row in 1:nrow(Kdomain)){
#     N=Kdomain[row,"N"]
#     M=Kdomain[row,"M"]
#     O=Kdomain[row,"O"]
#     if(M>-1){
#       this_marginal = marginal[which(marginal$N==N &marginal$M==M &marginal$O==O),]
#       from = round((N) /(N+M+alpha) * (length(domain) -1) + 1)
#       to =  round((N+alpha) / (N+M+alpha) * (length(domain) -1) + 1)
#       if(nrow(this_marginal) == 0){
#         next
#       } # FIX THIS, it hasn't been computed because stack depth has been exceeded.
#       output[from:to] = output[from:to] + rep(this_marginal$density,to-from+1)
#     }
#   }
#   return(output)
# }
# 
# get_posterior_GD_timing_distribution_smoothed <- function(dd){
#   output = get_posterior_GD_timing_distribution(dd)
#   domain=1:length(output)/1000
#   output2 = density(domain,weights=output/sum(output),from=min(domain),to=max(domain))
#   ol = length(output2$y)
#   output2$y = output2$y/(sum(output2$y[2:(ol-1)]) + 0.5*output2$y[1] + 0.5*output2$y[ol])*(ol-1)
#   return(output2)
# }
# 
# 
# get_posterior_GD_timing_distribution_plot <- function(dd){
#   density <- get_posterior_GD_timing_distribution(dd)
#   time = 1:length(density)/length(density)
#   return(ggplot() + geom_point(aes(time,density)))
# }
# 
# get_posterior_GD_timing_distribution_smoothed_plot <- function(dd){
#   density <- get_posterior_GD_timing_distribution_smoothed(dd)$y
#   time = 1:length(density)/length(density)
#   return(ggplot() + geom_point(aes(time,density)))
# }
# 
# get_posterior_GD_timing_ML_estimate <- function(dd){
#   output2 <- get_posterior_GD_timing_distribution_smoothed(dd)
#   output2$x[which(output2$y == max(output2$y))]
# }
# 
# get_posterior_GD_timing_smallest_CI_around_ML <- function(dd,CI){
#   # the aim here is to find the smallest area under which 95% of the area lays
#   output2 <- get_posterior_GD_timing_distribution_smoothed(dd)
#   sum(output2$y)
#   # we can assume that the output is unimodal and this makes things much easier:
#   xsys <- cbind(output2$x[order(output2$y)],output2$y[order(output2$y)])
#   xsys <- cbind(xsys,cumsum(xsys[,2]))
#   xsys <- xsys[max(which(xsys[,3]/nrow(xsys) < 1-CI)):nrow(xsys),]
#   return(c(min(xsys[,1]),max(xsys[,1])))
# }
# 
# 
# 
# 
# get_posterior_plot_2d <- function(dd){
#   library(cowplot)
#   cols <- rev(rainbow(7)[-7])
#   dd$prob_up <- dd$c
#   dd$prob_down <- dd$a
#   low = min(dd$density)
#   high = max(dd$density)
#   g1 <- ggplot(dd[which(dd$O == -1),], aes(prob_down, prob_up, fill= density)) + 
#     geom_tile() +
#     scale_fill_gradientn(colours = cols,limits=c(low,high))+ facet_grid(vars(N), vars(M))
#   
#   # need to find a way to make this uniformally blue.
#   g2 <- ggplot(dd[which(dd$O == 0),], aes(prob_down, prob_up, fill= density)) + 
#     geom_tile() +
#     scale_fill_gradientn(colours = cols,limits=c(low,high))+ facet_grid(vars(N), vars(M))
#   
#   g3 <- ggplot(dd[which(dd$O == 1),], aes(prob_down, prob_up, fill= density)) + 
#     geom_tile() +
#     scale_fill_gradientn(colours = cols,limits=c(low,high))+ facet_grid(vars(N), vars(M))
#   
#   g4 <- ggplot(dd[which(dd$O == 2),], aes(prob_down, prob_up, fill= density)) + 
#     geom_tile() +
#     scale_fill_gradientn(colours = cols,limits=c(low,high))+ facet_grid(vars(N), vars(M))
#   g5 <- ggplot(dd[which(dd$O == 3),], aes(prob_down, prob_up, fill= density)) + 
#     geom_tile() +
#     scale_fill_gradientn(colours = cols,limits=c(low,high))+ facet_grid(vars(N), vars(M))
#   g6 <- ggplot(dd[which(dd$O == 4),], aes(prob_down, prob_up, fill= density)) + 
#     geom_tile() +
#     scale_fill_gradientn(colours = cols,limits=c(low,high))+ facet_grid(vars(N), vars(M))
#   return(plot_grid(g1,g2,g3,g4,g5,g6))
# }
# 
# get_log_posterior_plot_2d <- function(dd){
#   library(cowplot)
#   cols <- rev(rainbow(7)[-7])
#   dd$prob_up <- dd$c
#   dd$prob_down <- dd$a
#   dd$log_density <- dd$res
#   high = max(dd$log_density)
#   low=2*high
#   g1 <- ggplot(dd[which(dd$O == -1),], aes(prob_down, prob_up, fill= log_density)) + 
#     geom_tile() +
#     scale_fill_gradientn(colours = cols,limits=c(low,high))+ facet_grid(vars(N), vars(M))
#   
#   # need to find a way to make this uniformally blue.
#   g2 <- ggplot(dd[which(dd$O == 0),], aes(prob_down, prob_up, fill= log_density)) + 
#     geom_tile() +
#     scale_fill_gradientn(colours = cols,limits=c(low,high))+ facet_grid(vars(N), vars(M))
#   
#   g3 <- ggplot(dd[which(dd$O == 1),], aes(prob_down, prob_up, fill= log_density)) + 
#     geom_tile() +
#     scale_fill_gradientn(colours = cols,limits=c(low,high))+ facet_grid(vars(N), vars(M))
#   
#   g4 <- ggplot(dd[which(dd$O == 2),], aes(prob_down, prob_up, fill= log_density)) + 
#     geom_tile() +
#     scale_fill_gradientn(colours = cols,limits=c(low,high))+ facet_grid(vars(N), vars(M))
#   return(plot_grid(g1,g2,g3,g4))
# }
# 
# 
