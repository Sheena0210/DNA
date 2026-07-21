#version1
library(data.table)
input<-paste0("/Users/sheena/Desktop/summerintern/","final.txt")
genotype_data<-fread(file=input,sep = "\t", header = TRUE, quote = "", na.strings = c("---","" ,"NA"),showProgress = TRUE)
dim(genotype_data)
head(genotype_data)

# 因為後續看chr start時發現三個Id是missimg 因此先進行補值 且position_source：hg_19----
target_ids <- c(
  "AX-123334869",
  "AX-158663204",
  "AX-123355693"
)
genotype_data[
  probeset_id %in% target_ids &
    (is.na(Chr_id) | trimws(as.character(Chr_id)) == "") &
    !is.na(hg19_chromosome) &
    trimws(as.character(hg19_chromosome)) != "",
  Chr_id := as.character(hg19_chromosome)
]

genotype_data[
  probeset_id %in% target_ids &
    is.na(Start) &
    !is.na(hg19_position),
  Start := as.numeric(hg19_position)
]

genotype_data[
  probeset_id %in% target_ids &
    !is.na(Chr_id) &
    trimws(as.character(Chr_id)) != "" &
    !is.na(Start),
  position_source := "hg_19"
]
#step1:見一個資料表,必要的欄位----
required_annotation<-c("probeset_id",
                       "dbSNP_RS_ID",
                       "Chr_id",
                       "Start",
                       "Allele_A",
                       "Allele_B",
                       "specialSNP_chr",
                       "hemizygous")
#看有沒有缺哨的欄位
data.frame(column=required_annotation,exists=required_annotation %in% names(genotype_data))

#找出10個樣本的genotype欄位
grep("DM|EOAD",names(genotype_data),value = TRUE)
sample_cols<-c("DM-002_(Axiom_TPM)_C03.CEL_call_code" ,"EOAD_P10_(Axiom_TPM)_G08.CEL_call_code",
              "EOAD_P11_(Axiom_TPM)_H08.CEL_call_code" ,"EOAD_P12_(Axiom_TPM)_A09.CEL_call_code",
              "EOAD_P13_(Axiom_TPM)_B09.CEL_call_code" ,"EOAD_P14_(Axiom_TPM)_C09.CEL_call_code",
              "EOAD_P15_(Axiom_TPM)_D09.CEL_call_code" ,"EOAD_P16_(Axiom_TPM)_E09.CEL_call_code",
              "EOAD_P17_(Axiom_TPM)_F09.CEL_call_code", "EOAD_P18_(Axiom_TPM)_G09.CEL_call_code")

#check genotype
genotype_data[, (sample_cols):= lapply(
              .SD,
              function(x) trimws(as.character(x))
              ),
              .SDcols= sample_cols]
allowed_genotypes <- c(
  "AA",
  "AB",
  "BB",
  "NoCall"
)
lapply(genotype_data[, ..sample_cols],function(x) sort(unique(x)))

#排序之前先新增列號
genotype_data[,original_row_id := .I]
genotype_data[1:6, .(original_row_id ,probeset_id,Chr_id,Start)]

#需要的annotation 
annotation_cols <- c(
  "original_row_id",
  "probeset_id",
  "dbSNP_RS_ID",
  "Chr_id",
  "Start",
  "Allele_A",
  "Allele_B",
  "specialSNP_chr",
  "hemizygous",
  "position_source"
)

#把annotation 跟10個genotye合併
master_cols<-c(annotation_cols,sample_cols)

#version1 table(680865*20)
snp_master_v1 <- genotype_data[, ..master_cols]
dim(snp_master_v1)
head(snp_master_v1)

#check prodeset_id
sum(is.na(snp_master_v1$probeset_id) | snp_master_v1$probeset_id=="")
sum(duplicated(snp_master_v1$probeset_id))

#check  dpsnp_rs_id
#missing:62515
sum(is.na(snp_master_v1$dbSNP_RS_ID) | snp_master_v1$dbSNP_RS_ID=="")
snp_master_v1[is.na(dbSNP_RS_ID)|dbSNP_RS_ID=="",.(original_row_id,probeset_id,dbSNP_RS_ID,Chr_id,Start)][1:10]

#check chr_id
#看染色體分佈
table(snp_master_v1$Chr_id, useNA = "ifany")
#missing=3 (AX-123355693, AX-123334869, AX-158663204)
sum(is.na(snp_master_v1$Chr_id)| snp_master_v1$Chr_id=="")
missing_chr <-snp_master_v1[is.na(Chr_id)| Chr_id==""]

#check start
#missing=3 (AX-123355693 AX-123334869 AX-158663204)
sum(is.na(snp_master_v1$Start))
missing_start <-snp_master_v1[is.na(Start)| Start==""]


#確認chr_id start 在genotype_data這四個有沒有hg19
sum(is.na(genotype_data$hg19_chromosome) | genotype_data$hg19_chromosome=="")
sum(is.na(genotype_data$hg19_position) | genotype_data$hg19_position=="")

target<-c("AX-123334869", "AX-158663204", "AX-123355693")
check<-genotype_data[probeset_id %in% target,
                     .(probeset_id,
                       Chr_id,
                       Start,
                       hg19_chromosome,
                       hg19_position,
                       
                       hg19_chr_missing =
                         is.na(hg19_chromosome) |
                         trimws(as.character(hg19_chromosome)) == "",
                       
                       hg19_position_missing =
                         is.na(hg19_position) |
                         trimws(as.character(hg19_position)) == ""
                     )
]
check[,
  both_hg19_missing :=
    hg19_chr_missing & hg19_position_missing
]

setdiff(
  target,
  genotype_data$probeset_id
)

#check allele_a allele_b
sum(is.na(snp_master_v1$Allele_A)|snp_master_v1$Allele_A=="")
sum(is.na(snp_master_v1$Allele_B)|snp_master_v1$Allele_B=="")

#check specialSNP_chr
sum(is.na(snp_master_v1$specialSNP_chr)|snp_master_v1$specialSNP_chr=="")

#check hemizygous
sum(is.na(snp_master_v1$hemizygous)|snp_master_v1$hemizygous=="")

genotype_count_v1 <- rbindlist(
  lapply(sample_cols, function(sample_name) {
    
    counts <- table(
      factor(
        snp_master_v1[[sample_name]],
        levels = allowed_genotypes
      )
    )
    
    data.table(
      sample = sample_name,
      AA = as.integer(counts["AA"]),
      AB = as.integer(counts["AB"]),
      BB = as.integer(counts["BB"]),
      NoCall = as.integer(counts["NoCall"]),
      R_NA = sum(is.na(snp_master_v1[[sample_name]])),
      total = nrow(snp_master_v1)
    )
  })
)
#建立摘要表
library(data.table)

step1_summary <- data.table(
  item = c(
    "Number of SNP rows",
    "Number of master-table columns",
    "Missing probeset_id",
    "Duplicated probeset_id",
    "Missing dbSNP_RS_ID",
    "Missing Chr_id",
    "Missing Start",
    "Start <= 0"
  ),
  
  value = c(
    nrow(snp_master_v1),
    ncol(snp_master_v1),
    
    sum(
      is.na(snp_master_v1$probeset_id) |
        trimws(as.character(snp_master_v1$probeset_id)) == ""
    ),
    
    sum(
      duplicated(snp_master_v1$probeset_id)
    ),
    
    sum(
      is.na(snp_master_v1$dbSNP_RS_ID) |
        trimws(as.character(snp_master_v1$dbSNP_RS_ID)) == ""
    ),
    
    sum(
      is.na(snp_master_v1$Chr_id) |
        trimws(as.character(snp_master_v1$Chr_id)) == ""
    ),
    
    sum(
      is.na(snp_master_v1$Start)
    ),
    
    sum(
      !is.na(snp_master_v1$Start) &
        snp_master_v1$Start <= 0
    )
  )
)

step1_summary
#見一個資料夾 保存parallel 的資料 
output_dir <- "/Users/sheena/Desktop/summerintern/parallel2"
dir.create(output_dir,showWarnings = FALSE,recursive = TRUE)

#將摘要表存成txt
fwrite(snp_master_v1,file = file.path(output_dir,"snp_master_v1_step1.txt"),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)
  

#step2:排序染色體-----
#避免覆蓋 所以先複製檔案
snp_master_v2 <-copy(snp_master_v1)

#看目前有哪些染色體
table(snp_master_v2$Chr_id,useNA = "ifany")
#定義染色體順序
valid_chr<-c(as.character(1:22),"X","Y","MT")
unexpected_chr <- setdiff(
  unique(snp_master_v2$Chr_id),
  valid_chr
)

unexpected_chr

snp_master_v2[!Chr_id %in% valid_chr,
  .(
    original_row_id,
    probeset_id,
    dbSNP_RS_ID,
    Chr_id,
    Start
  )
]

#檢查"7_KI270803v1_alt"     "8_KI270821v1_alt"     "22_KI270879v1_alt"    "1_KI270706v1_random"  "1_KI270766v1_alt"     "19_KI270938v1_alt"    "14_GL000009v2_random" "4_GL000008v2_random" 
#在hg19_chromosome、hg19_position 是否有主染色體位置 才能套用在Chr_id中
#利用genotype_data
noncanonical_chr <- c(
  "7_KI270803v1_alt",
  "8_KI270821v1_alt",
  "22_KI270879v1_alt",
  "1_KI270706v1_random",
  "1_KI270766v1_alt",
  "19_KI270938v1_alt",
  "14_GL000009v2_random",
  "4_GL000008v2_random"
)

#找出probeset_id
noncanonical_snps<- snp_master_v2[Chr_id %in% noncanonical_chr,
  .(original_row_id,
    probeset_id,
    dbSNP_RS_ID,
    Chr_id,
    Start
  )
]
#共有21個snp 的chr_id start無法對回hg19 因此先排除 等確認後在加入後續----
nrow(noncanonical_snps)
#去查位置
check_hg19<-genotype_data[probeset_id %in% noncanonical_snps$probeset_id,.(probeset_id,Chr_id,Start,hg19_chromosome,hg19_position,position_source)]
setorder(check_hg19,Chr_id,Start)
check_hg19

#排除21個snp的snp_master_v2->snp_master_v2_n (snp:680844)----
#另外保存這21個snp
excluded_noncanonical <- snp_master_v2[Chr_id %in% noncanonical_chr]
#排除
snp_master_v2_n <- copy(snp_master_v2[!Chr_id %in% noncanonical_chr])
dim(snp_master_v2_n)




#排序染色體 原本chr_id是文字型 建立一個chr_level 使"X","Y","MT"＝23 24 25
chr_level<-c(as.character(1:22),"X","Y","MT")
#依據chr1-22 x y mt排序->chr_order
snp_master_v2_n[,chr_order := match(Chr_id, chr_level)]
sum(is.na(snp_master_v2_n$chr_order))


#先按照染色體順序排列 同一條染色體維持原始順序
setorder(snp_master_v2_n,chr_order,original_row_id)
head(snp_master_v2_n)








#step3：排序start----
#chr1-22/x/y/mt 且start從小到大
#先將start轉乘數值
snp_master_v2_n[,Start :=as.numeric(trimws(Start))]
class(snp_master_v2_n$Start)
setorder(snp_master_v2_n,chr_order,Start,probeset_id)
head(snp_master_v2_n)
#建立染色體跟start排序完後 依照的順序>作為後續影像順序 image_order(從1-680844)----
snp_master_v2_n[,image_order:=.I]

#新增每條染色體內的SNP順序  SNP 在自己所屬染色體內的位置順序 方便之後查看某個 SNP 是該染色體內第幾個
snp_master_v2_n[,snp_order_within_chr:=seq_len(.N),by=.(chr_order,Chr_id)]

#計算相鄰snp距離->distance----
#同一條染色體內，依 Start 排序後，前後相接的兩個 array SNP(距離單位為 bp)
#如果同一個染色體沒有前一個snp/換到下一個染色體  則會是na 
snp_master_v2_n[,distance:=Start-shift(Start), by=Chr_id]
#看NA的情況:25個
snp_master_v2_n[,.(snp_count=.N, distance=sum(is.na(distance)))
                ,by=.(chr_order, Chr_id)][order(chr_order)]

#考慮不同snp 但是在停一條染色體上 且 同一個基因組位置->distance:=Start-shift(Start)=0----
sum(snp_master_v2_n$distance==0,na.rm=TRUE)
snp_master_v2_n[distance==0,.(image_order,probeset_id,dbSNP_RS_ID,Chr_id,Start,distance)]

#利用Chr_id, Start, probeset_id, dbSNP_RS_ID, affy_snp_id, BestProbeset
#有3444個位置重複 7502個probeset->皆保留 且 新增欄位紀錄排列 不是排序----
same_position_key<-snp_master_v2_n[,.N,by=.(Chr_id,Start)][N>1,.(Chr_id,Start)]
same_position_ids<-snp_master_v2_n[same_position_key,on=.(Chr_id,Start),probeset_id]

same_position_check<-genotype_data[probeset_id %in% same_position_ids,.(Chr_id,Start,probeset_id,dbSNP_RS_ID,affy_snp_id,BestProbeset)]
setorder(same_position_check,Chr_id,Start,affy_snp_id,probeset_id)
same_position_summary<-same_position_check[,.(n_probeset=.N,
     n_affy_snp_id=uniqueN(affy_snp_id),
     n_rsid=uniqueN(dbSNP_RS_ID,na.rm=TRUE),
     BestProbeset_values=paste(unique(BestProbeset),collapse=",")),
  by=.(Chr_id,Start)]

same_position_summary


#同一位置共有幾個marker
snp_master_v2_n[,same_position_n:=.N,by=.(Chr_id,Start)]

#同一位置內的排列順序
snp_master_v2_n[,same_position_rank:=seq_len(.N),by=.(Chr_id,Start)]

snp_master_v2_n[same_position_n>1,
                .(Chr_id,Start,probeset_id,dbSNP_RS_ID,same_position_n,same_position_rank,image_order)]


same_position_distribution<-snp_master_v2_n[
  ,.(same_position_n=.N),by=.(Chr_id,Start)
][same_position_n>1,
  .N,by=same_position_n
][order(same_position_n)]

same_position_distribution









#step4:確認13404個nocall----

#每個SNP在10位樣本中有幾個NoCall
snp_master_v2_n[,missing_count_10:=rowSums(.SD=="NoCall"),.SDcols=sample_cols]

table(snp_master_v2_n$missing_count_10)

sum(snp_master_v2_n$missing_count_10==6)
six_nocall_chr<-snp_master_v2_n[missing_count_10==6,.N,by=.(chr_order,Chr_id)][order(chr_order)]
six_nocall_chr

snp_master_v2_n[missing_count_10==6 & Chr_id=="Y",.N]
six_nocall<-snp_master_v2_n[missing_count_10==6]
six_nocall[,missing_pattern:=apply(.SD,1,function(x) paste(sample_cols[x=="NoCall"],collapse=" | ")),.SDcols=sample_cols]

six_nocall_pattern<-six_nocall[,.N,by=missing_pattern][order(-N)]

six_nocall_pattern

nrow(six_nocall_pattern)


y_snp_summary<-snp_master_v2_n[Chr_id=="Y",.(total_Y_SNP=.N,
  six_sample_NoCall=sum(missing_count_10==6),
  other_missing_pattern=sum(missing_count_10!=6)
)]

y_snp_summary

snp_master_v2_n[Chr_id=="Y" & missing_count_10!=6,
                .(image_order,probeset_id,dbSNP_RS_ID,Chr_id,Start,missing_count_10)]

#發現14個snp 來自7人 nocall 也是y染色體的部分
# 篩選出位於Y染色體，而且10位樣本中有7位為NoCall的14個SNP
seven_nocall<-snp_master_v2_n[Chr_id=="Y"&missing_count_10==7]

# 對每一個SNP逐列檢查10位樣本
# 找出genotype為NoCall的樣本名稱
# 再用「 | 」將7位NoCall樣本名稱合併成一個字串
seven_nocall[,missing_pattern:=apply(
  .SD,
  1,
  function(x){
    paste(sample_cols[x=="NoCall"],collapse=" | ")
  }
),.SDcols=sample_cols]

# 統計每一種NoCall樣本組合出現幾個SNP
# N代表該missing pattern出現的SNP數量
seven_nocall_pattern<-seven_nocall[,.N,by=missing_pattern][order(-N)]

# 顯示結果
seven_nocall_pattern

#找出10位樣本的genotype call欄位
sample_cols<-grep("_call_code$",names(genotype_data),value=TRUE)

#計算每位樣本在male類別SNP中的genotype成功率
#如果某人在 male 類別 SNP 幾乎都有 genotype，通常表示具有 Y 染色體
#如果某人在這些 SNP 幾乎全部是 NoCall，通常表示沒有 Y 染色體
male_marker_qc<-rbindlist(lapply(sample_cols,function(x){
  
  #取出目前樣本在male類別SNP中的genotype
  calls<-genotype_data[gender_metrics=="male",get(x)]
  
  data.table(
    sample=x,
    male_marker_total=length(calls),
    male_marker_called=sum(!is.na(calls)&calls!="NoCall"),
    male_marker_NoCall=sum(is.na(calls)|calls=="NoCall"),
    male_marker_call_rate=round(
      mean(!is.na(calls)&calls!="NoCall")*100,
      3
    )
  )
}))

male_marker_qc


