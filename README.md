Chen Wang
2022-03-22

RATES
=====

**R**eplicability **A**ssessment in **T**rans-Ancestry **G**enetic **S**tudies

RATES is designed to assess repliabilty of association signals from trans-ancestry GWAS meta-analysis by calculating a posterior probability of replicability for each analyzed variant.

Installation
------------

    library(devtools)
    install_github("wangc29/RATES/rates")

Detailed example
--------------

We will use an example dataset to illurstrate the procedures for using RATES to assess the replicality of association signals identified from trans-ancestry GWAMA (Genome-wide association meta-analysis). The example dataset consists of three pieces of summary statistics, i.e. effect size estimate matrix, standard error estimate matrix and allele frequency matrix. For each matrix, each row represents a variant and each column represents a participating study in meta-analysis. These studies were conducted in populations of five continental ancestries including European, African American, East Asian and American ancestries.

**Step 1: Initial trans-ancestry GWAMA**

By utilizing the summary statistics data mentioned above, we first conduct a trans-ancestry GWAMA using the METASOFT software to marginally associated each variant to the trait of interest. Then we load the mariginal assocation results in R and retain the result from RE2 (an improved random effects model implemented in METASOFT) method for downstream analysis.

```{r}
meta.res.dt<-read.table("./meta_demo_res_cleaned.txt",header = T,sep = "\t",stringsAsFactors = F)
meta_re2_dt<-meta.res.dt[,c("Variant","CHR","BP","PVALUE_RE2")]
###Remove NA
meta_re2_dt<-meta_re2_dt[!is.na(meta_re2_dt[,"PVALUE_RE2"]),]
head(meta_re2_dt)
#            Variant  CHR       BP PVALUE_RE2
#1 chr2:79308790_T/C chr2 79308790   0.973941
#2 chr2:79308954_A/G chr2 79308954   0.366667
#3 chr2:79308981_T/G chr2 79308981   0.375199
#4 chr2:79309040_A/G chr2 79309040   0.931090
#5 chr2:79309067_A/G chr2 79309067   0.436295
#6 chr2:79309082_T/A chr2 79309082   0.809719
```
We set significant threshold for association as 1e-5 and split the results into significantly associated variants and null variants

```{r}
meta_sig_re2_dt<-meta_re2_dt[meta_re2_dt[,"PVALUE_RE2"]<=sig_thres,]
meta_null_re2_dt<-meta_re2_dt[meta_re2_dt[,"PVALUE_RE2"]>sig_thres,]
```

**Step 2: Identify sentinel variants and genetic loci through clumping**

Then we conduct a distance based clumping procedure to identify sentinel variants and corresponding genetic loci based on those significantly associated variants. Briefly, the most significantly associated variant in each genetic locus is chosen as sentinel. The corresponing genetic locus will span a predefined range on each side of the sentinel variant. We implemented this procedure into the following fuction, which can be utilized to identify sentinel variants by chromosomes

```{r}
PriorityQueue <- function() {
  keys <- values <- NULL
  insert <- function(key, value) {
    ord <- findInterval(key, keys)
    keys <<- append(keys, key, ord)
    values <<- append(values, value, ord)
  }
  pop <- function() {
    head <- list(key=keys[1],value=values[[1]])
    values <<- values[-1]
    keys <<- keys[-1]
    return(head)
  }
  empty <- function() length(keys) == 0
  environment()
}
is.merged<-function(cand.bin,interval.list){
  key.vec<-c()
  old_bin.vec<-c()
  for (key in names(interval.list)) {
    old_bin<-interval.list[[key]]
    if ((cand.bin[1]<old_bin[2]) & (cand.bin[2]>old_bin[1])) {
      key.vec<-c(key.vec,key)
      old_bin.vec<-c(old_bin.vec,paste0(old_bin[1],"-",old_bin[2]))
    }
  }
  merged_res<-as.data.frame(list(old_senti=key.vec,old_bin=old_bin.vec),stringsAsFactors =F)
  if(dim(merged_res)[1]!=0){
    return(merged_res)
  }
  return(NULL)
}

get.sentinelbyChrom<-function(sig.hit.df,bin_size=5e5){
  ##prepare the dataframe
  colnames(sig.hit.df)<-c('UNIQ_ID',"CHR","BP",'PVALUE')
  ##
  priori<-rank(sig.hit.df[,"PVALUE"],ties.method = "first")
  sig.hit.df[,"PRIORITY"]<-priori
  ##create a pq of variants
  pq<-PriorityQueue()
  for (p in 1:length(priori)) {
    pq$insert(p,sig.hit.df[sig.hit.df$PRIORITY==p,]$UNIQ_ID)
  }
  #create a list to hold sentinels
  S<-vector("list")
  #loop to insert proper value to S.
  while (!pq$empty()) {
    senti<-pq$pop()
    s<-senti$value
    pos<-as.numeric(sig.hit.df[sig.hit.df[,"UNIQ_ID"]==s,"BP"])
    bin<-c(-bin_size,bin_size)+pos
    if (bin[1]<0) {
      bin[1]=0
    }
    merged.res<-is.merged(bin,S)
    if (is.null(merged.res)) {
      S[[s]]<-list(start=bin[1],end=bin[2])
      next
    }
    old_s<-merged.res[,1]
    s.all<-c(old_s,s)
    old_start<-as.numeric(as.character(str_split(merged.res[,2],"-",simplify = T)[,1]))
    old_end<-as.numeric(as.character(str_split(merged.res[,2],"-",simplify = T)[,2]))
    start.all<-c(old_start,bin[1])
    end.all<-c(old_end,bin[2])
    start_updated<-min(start.all)
    end_updated<-max(end.all)
    for(i in 1:length(s.all)){
      S[[s.all[i]]]<-NULL
    }
    pvalue.all<-as.numeric(sig.hit.df$PVALUE[sig.hit.df$UNIQ_ID %in% s.all])
    pvalue.all<-sig.hit.df$PVALUE[match(s.all,sig.hit.df$UNIQ_ID)]
    ix<-which.min(pvalue.all)
    new_senti<-s.all[ix]
    S[[new_senti]]<-list(start=start_updated,end=end_updated)
  }
  return(S)
}
```
We perform clumping with those significantly associated variants. The width of genetic locus is defined as 5e4 on each side of sentinel variants.

```{r}
senti_list= get.sentinelbyChrom(sig.hit.df=meta_sig_re2_dt,bin_size = 5e4)
senti.df<-data.frame(list(start=NULL,end=NULL),stringsAsFactors = F)
for (s in names(senti_list)) {
        senti.df<-rbind(senti.df,data.frame(senti_list[[s]]))
}
senti.df=cbind(Variant=names(senti_list),senti.df)
senti.df[,"CHR"]=str_split(senti.df[,1],pattern = ":",simplify = TRUE)[,1]
senti.df<-senti.df[,c("Variant","CHR","start","end")]
senti.df
#             Variant  CHR    start      end
# 1 chr2:79465508_A/G chr2 79415508 79515508
# 2 chr2:80786611_C/G chr2 80677951 81002790
# 3 chr2:80465766_C/T chr2 80382319 80518505
# 4 chr2:81658818_A/G chr2 81464626 81776824
```
For the example data, four sentinel variants, i.e. four genetic loci were identified.

**Step 3: Pruning of null variants**

In this step we select a number of null variants that are mutually independent and are independent from identified sentinels. These selected null variants will be jointly analyzed by RATES for replicability assessment.

We first obtained a series of variants that are mutally indenpendent by pruning against the HRC reference panel using the "“indep-pairwise" function in PLINK, and load these variants in R. Since the example data is from chromosome 2, we will only retain those pruned variants on chromosome 2.

```{r}
pruned_var_file<-"./hg38_HRC_pruned_lifted.bed"
pruned_var_dt<-read.table(pruned_var_file,stringsAsFactors = F,header = F)
pruned_var_selected<-pruned_var_dt[pruned_var_dt[,1]=="chr2",]
head(pruned_var_selected)
#          V1    V2    V3
# 129290 chr2 12821 12821
# 129291 chr2 16078 16078
# 129292 chr2 30703 30703
# 129293 chr2 38021 38021
# 129294 chr2 76530 76530
# 129295 chr2 91076 91076
```
Columns are chromosome, start, and end position of those pruned variants. Next, we selected those pruned null variants that failed to reach specified significance threshold in GWAMA.

```{r}
pruned_null_var= pruned_var_selected[which(pruned_var_selected[,2] %in% meta_null_re2_dt$BP),]
nrow(pruned_null_var)
#218
```
In total 218 pruned variants failed to achieve specified significance threshold in GWAMA. Then we remove those null variants that fall in identified genetic loci to ensure that these null variants are also indenpendent from sentinel variants.

```{r}
library(seqminer)
library(pracma)
pruned_null_pos<-paste0(pruned_null_var[,1],":",pruned_null_var[,2])
indpt_loci<-paste0(senti.df$CHR,":",senti.df$start,"-",senti.df$end)
indpt_null_pos<-pruned_null_pos[!isInRange(pruned_null_pos,paste(indpt_loci,collapse = ","))]
indpt_null_pos
senti_pos<-str_split(senti.df$Variant,pattern = "_",simplify = TRUE)[,1]
###
jointly_analyzed_var= as.data.frame(list(Variant_pos=c(senti_pos,indpt_null_pos),is.sentinel=c(rep(1,length(senti_pos)),rep(0,length(indpt_null_pos)))))
head(jointly_analyzed_var)
#     Variant_pos is.sentinel
# 1 chr2:79465508           1
# 2 chr2:80786611           1
# 3 chr2:80465766           1
# 4 chr2:81658818           1
# 5 chr2:79317210           0
# 6 chr2:79318025           0
nrow(jointly_analyzed_var)
#184
```
These 184 variants (4 sentinels and 180 null variants) will be jointly analyzed by RATES in following procedures.

**Step 4: Replicability assessment by RATES**

First we gather the summary statistics for those variants that were selected to be analyzed by RATES.

```{r}
beta.mt<-read.table("./beta_demo_mt.txt",stringsAsFactors = F,header = TRUE)
sd.mt<-read.table("./sd_demo_mt.txt",stringsAsFactors = F,header = TRUE)
af.mt<-read.table("./af_demo_mt.txt",stringsAsFactors = F,header = TRUE)
###
match.ix<-match(jointly_analyzed_var$Variant_pos,str_split(beta.mt[,1],pattern = "_",simplify = TRUE)[,1])
beta.retained<-beta.mt[match.ix,]
sd.retained<-sd.mt[match.ix,]
```

Next we obtain the genome-wide allele frequency PCs that are needed by RATES to capture variation of genetic effects across populations due to ancestry.

```{r}
dist_mt=dist(t(af.mt[,-1]))
raw_PC = cmdscale(dist_mt,k = 3)
PC<-cbind(1,raw_PC)
```

Finally, lets run RATES to for replicability assessment for those signals identified from trans-ancestry GWAMA.
```{r}
library(rates)
res<-rates.fit(betajk = beta.retained[,-1],sjk2 = (sd.retained[,-1])^2,PC = PC1,SNP = beta.retained[,1])
```
After fitting the model, we can obtain the Posterior Probability of Replicability estimate for all jointly analyzed variants. 

```{r}
res.dt= as.data.frame(list(SNP=res$SNP,PPR=res$ppr),stringsAsFactors = F)
```
Select those variants of interest, i.e. sentinel variants,

```{r}
res.dt[,"is.sentinel"]=0
res.dt$is.sentinel[match(senti.df$Variant,res.dt$SNP,)]<-1
res.dt[res.dt$is.sentinel==1,]
#                 SNP       PPR is.sentinel
# 1 chr2:79465508_A/G 0.3239346           1
# 2 chr2:80786611_C/G 0.9999991           1
# 3 chr2:80465766_C/T 0.9997391           1
# 4 chr2:81658818_A/G 0.9998328           1
```

Out of four identified sentinel variants, three sentinels have PPR>0.99. According to RATES, these three variants are more likely to be replicated in a well powered replication dataset.