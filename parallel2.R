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


#
#根據Y染色體marker的call rate推定樣本性別型態
#>=90%：具有明顯Y染色體訊號
#<=10%：沒有Y染色體訊號
#介於兩者之間：暫時標記為Uncertain
male_marker_qc[,inferred_sex:=fifelse(male_marker_call_rate>=90,"Male",fifelse(male_marker_call_rate<=10,"Female","Uncertain"))]

#顯示每位樣本的Y marker結果與推定性別
male_marker_qc[,.(sample,male_marker_total,male_marker_called,male_marker_NoCall,male_marker_call_rate,inferred_sex)]

#沒有Y染色體訊號的6位樣本
female_samples<-male_marker_qc[inferred_sex=="Female",sample]

#具有Y染色體訊號的4位樣本
male_samples<-male_marker_qc[inferred_sex=="Male",sample]

#確認分組人數
length(female_samples)
length(male_samples)

#查看實際樣本名稱
female_samples
male_samples


#影像使用的 genotype 分類表 把「沒有 Y 染色體造成的 NoCall」和「真正判讀失敗的 NoCall」分開----
#處理欄位
#逐一處理6位沒有Y染色體訊號的樣本
for(x in female_samples){
  
  #只有在Y染色體且原本為NoCall時
  #才改成Y_not_applicable
  snp_master_v2_n[
    Chr_id=="Y"&get(x)=="NoCall",
    (x):="Y_not_applicable"
  ]
}

#合併10位樣本的分類值
all_image_classes<-unlist(snp_master_v2_n[,..sample_cols],use.names=FALSE)
#查看全部AA、AB、BB、NoCall及Y_not_applicable數量
table(all_image_classes,useNA="ifany")

#合併10位樣本的分類值
all_image_classes<-unlist(snp_master_v2_n[,..sample_cols],use.names=FALSE)

#查看全部AA、AB、BB、NoCall及Y_not_applicable數量
table(all_image_classes,useNA="ifany")

head(snp_master_v2_n)
#統計每位樣本在Y染色體上的分類數量
y_image_class_by_sample<-rbindlist(
  lapply(sample_cols,function(x){
    
    snp_master_v2_n[
      Chr_id=="Y",
      .(
        sample=x,
        AA=sum(get(x)=="AA"),
        AB=sum(get(x)=="AB"),
        BB=sum(get(x)=="BB"),
        NoCall=sum(get(x)=="NoCall"),
        Y_not_applicable=sum(get(x)=="Y_not_applicable")
      )
    ]
  })
)

y_image_class_by_sample

#建立座標









#step5 建立圖像顏色 兩種版本----
#考慮是否要把這兩種分兩個顏色 因為看起來 如果兩種nocall上兩個顏色感覺是提供了性別的資訊 因此可能模型訓練會overfitting
#因此做兩種
#v1兩種nocall同色
#v2兩種nocall不同色

#Version A：兩種 NoCall 使用同顏色(灰色)
#"Padding","AA","AB","BB","NoCall","Y_not_applicable" : black, red, green, blue, gray
color_encoding_A<-data.table(
  genotype_class=c("Padding","AA","AB","BB","NoCall","Y_not_applicable"),code=c(0,1,2,3,4,4), 
  R=c(0,255,0,0,128,128),
  G=c(0,0,255,0,128,128),
  B=c(0,0,0,255,128,128))
color_encoding_A


#Version B：一般NoCall與Y_not_applicable分開編碼
#"Padding","AA","AB","BB","NoCall","Y_not_applicable" : black, red, green, blue, gray, yellow
color_encoding_B<-data.table(
  genotype_class=c("Padding","AA","AB","BB","NoCall","Y_not_applicable"),code=c(0,1,2,3,4,5),
  R=c(0,255,0,0,128,255),
  G=c(0,0,255,0,128,255),
  B=c(0,0,0,255,128,0))
color_encoding_B


#step6 將genotype轉數值(AA/AB/BB->1/2/3 最初版)----
#genotype_code_A：兩種 NoCall 合併
#genotype_code_B：兩種 NoCall 分開
#先排序,每位樣本都必須按照相同的 image_order 轉換，才能確保所有人的第 1 個數值都對應同一個 SNP
head(snp_master_v2_n)
#依照先前建立的image_order排序，確保兩個版本使用完全相同的SNP順序
#確認目前資料列順序與image_order完全一致
stopifnot(all(snp_master_v2_n$image_order==seq_len(nrow(snp_master_v2_n))))
#確認目前共有多少個SNP
nrow(snp_master_v2_n)

#確認 genotype 類別沒有異常值
#合併10位樣本的genotype欄位，找出實際出現的所有類別
observed_classes<-sort(unique(unlist(snp_master_v2_n[,..sample_cols],use.names=FALSE)))
observed_classes
#檢查是否有不在顏色編碼表中的異常類別
setdiff(observed_classes,color_encoding_B$genotype_class)

#建立數值轉換規則----
#version1(color_encoding_A/genotype_code_A)----
#AA->1
#AB->2
#BB->3
#NoCall->4
#Y_not_applicable->4
#NoCall,Y_not_applicable轉成4
code_map_A<-setNames(color_encoding_A$code,color_encoding_A$genotype_class)
#文字轉成數值矩陣
#按照image_order把文字類別轉成Version A數值
#結果為：row=SNP，column=sample
genotype_code_A<-sapply(sample_cols,function(x){unname(code_map_A[snp_master_v2_n[[x]]])})

#將矩陣資料型態固定為整數
storage.mode(genotype_code_A)<-"integer"


#Version 2(color_encoding_B/genotype_code_B)----
#AA->1
#AB->2
#BB->3
#NoCall->4
#Y_not_applicable->5
#NoCall轉成4, Y_not_applicable轉成5
code_map_B<-setNames(color_encoding_B$code,color_encoding_B$genotype_class)
#轉序列
#按照image_order把文字類別轉成Version B數值
#結果為：row=SNP，column=sample
genotype_code_B<-sapply(sample_cols,function(x){unname(code_map_B[snp_master_v2_n[[x]]])})
storage.mode(genotype_code_B)<-"integer"


#Step 7：轉成二維矩陣->把原本的一長串 SNP，整理成 CNN 能讀取的影像格式----
#計算影像大小:826*826----
#Padding:1432----
#取得目前每位樣本的SNP數量
n_snp<-nrow(genotype_code_A)
#計算能放下全部SNP的最小正方形邊長
image_side<-ceiling(sqrt(n_snp))
#計算補成完整正方形需要增加多少個Padding
padding_n<-image_side^2-n_snp
#確認Version A與Version B的矩陣大小完全相同
stopifnot(identical(dim(genotype_code_A),dim(genotype_code_B)))











#versionA轉影像矩陣----
#逐一處理Version A中的10位樣本
image_matrix_A<-lapply(seq_len(ncol(genotype_code_A)),function(i){
  #取出第i位樣本依image_order排列的680844個SNP數值
  snp_vector<-genotype_code_A[,i]
  #在序列最後補上1432個Padding，Padding代碼為0
  padded_vector<-c(snp_vector,rep(0L,padding_n))
  #按照由左到右、由上到下的順序排成826×826矩陣
  matrix(padded_vector,nrow=image_side,ncol=image_side,byrow=TRUE) #byrow=TRUE 希望排序方式第1個SNP → 第1列第1格 第2個SNP → 第1列第2格 而不是default從上到下
})

#使用原始樣本欄位名稱命名每個影像矩陣
names(image_matrix_A)<-colnames(genotype_code_A)

#versionB轉影像矩陣----
#逐一處理Version B中的10位樣本
image_matrix_B<-lapply(seq_len(ncol(genotype_code_B)),function(i){
  #取出第i位樣本的SNP數值序列
  snp_vector<-genotype_code_B[,i]
  #在序列最後補上Padding代碼0
  padded_vector<-c(snp_vector,rep(0L,padding_n))
  #byrow=TRUE:按照由左到右、由上到下排列成影像矩陣
  matrix(padded_vector,nrow=image_side,ncol=image_side,byrow=TRUE)
})

#使用原始樣本欄位名稱命名每個影像矩陣
names(image_matrix_B)<-colnames(genotype_code_B)

#確認兩個版本都包含10位樣本
length(image_matrix_A)
length(image_matrix_B)

#確認第一位樣本的影像尺寸
dim(image_matrix_A[[1]])
dim(image_matrix_B[[1]])

#check snp排序正確
#比較原始數值序列前10個SNP與影像第一列前10格-> is the same
genotype_code_A[1:10,1]
image_matrix_A[[1]][1,1:10]

#確認第一列與第二列的銜接
#第826個SNP應位於第一列最後一格
genotype_code_A[826,1]
image_matrix_A[[1]][1,826]

#第827個SNP應位於第二列第一格
genotype_code_A[827,1]
image_matrix_A[[1]][2,1]

#確認最後一個SNP與Padding
#located at 第825列、第220欄
genotype_code_A[n_snp,1]
image_matrix_A[[1]][825,220]
#最後一個SNP的下一格開始 Padding=0
image_matrix_A[[1]][825,221]









#Step8:數值矩陣轉成RGB三通道矩陣----
#code_matrix：一位樣本的826×826數值矩陣
#color_table：顏色編碼表(color_encoding_A, color_encoding_B)

code_to_rgb<-function(code_matrix, color_table){
  #取出code R G B四個欄位
  #unique(...,by="code")：每個code只保留一組RGB 並且照code小到大排列
  lookup<-unique(color_table[,.(code,R,G,B)],by="code")[order(code)]
  #將code統一轉成整數，避免numeric與integer型態不同
  lookup[,code:=as.integer(code)]
  #建立預期的連續code，例如0、1、2、3、4
  expected_codes<-seq.int(0L,max(lookup$code))
  #確認顏色表沒有缺少任何code
  stopifnot(identical(lookup$code,expected_codes))
  #三維矩陣, 初始值=0  前兩維度跟原本數值矩陣相同
  #第三維度 放R/G/B
  rgb_array<-array(0, dim=c(nrow(code_matrix),ncol(code_matrix),3))
  
  #轉成對應red 
  #code_matrix+1L：因為R索引從1開始，但code從0開始
  #例如code 0要對應lookup第1列，因此要加1
  #除以255：將0–255轉成影像常用的0–1
  rgb_array[,,1]<-matrix(lookup$R[code_matrix+1L],nrow=nrow(code_matrix))/255 #1L因為 R 的索引從 1 開始，但你的 code 從 0 開始
  #將每一格code轉成對應的綠色數值
  rgb_array[,,2]<-matrix(lookup$G[code_matrix+1L],nrow=nrow(code_matrix))/255 #1L 中的 L 代表 integer，亦即整數 1。
  #將每一格code轉成對應的藍色數值
  rgb_array[,,3]<-matrix(lookup$B[code_matrix+1L],nrow=nrow(code_matrix))/255
  
  #第三維加上名稱 第一層叫R、第二層叫G、第三層叫B
  dimnames(rgb_array)<-list(NULL,NULL,c("R","G","B"))
  rgb_array}

#轉換 Version A
#逐一將10位樣本的Version A數值矩陣轉成RGB
#NoCall與Y_not_applicable → 灰色
rgb_image_A<-lapply(image_matrix_A,function(x){code_to_rgb(x,color_encoding_A)})

#轉換Version B
#逐一將10位樣本的Version B數值矩陣轉成RGB
#NoCall → 灰色
#Y_not_applicable → 黃色
rgb_image_B<-lapply(image_matrix_B,function(x){code_to_rgb(x,color_encoding_B)})

#check
#--1.選擇第一位樣本
check_sample<-sample_cols[1]
#查看第一格原本的genotype code
image_matrix_A[[check_sample]][1,1]

#--2.查看第一格轉換後的RGB，乘回255方便判讀顏色
round(rgb_image_A[[check_sample]][1,1,]*255)
#影像最後一格沒有對應SNP，因此應為Padding code 0
image_matrix_A[[check_sample]][826,826]
#Padding的RGB應為黑色：0、0、0
round(rgb_image_A[[check_sample]][826,826,]*255)

#--3.女性型 Y 位點的差異
first_Y_order<-snp_master_v2_n[Chr_id=="Y",min(image_order)]
#將image_order換算成826×826矩陣的位置
first_Y_row<-((first_Y_order-1L)%/%image_side)+1L
first_Y_col<-((first_Y_order-1L)%%image_side)+1L
#選擇第一位沒有Y染色體訊號的樣本
check_female<-female_samples[1]
#比較同一個Y位點在兩個版本中的code與RGB
Y_color_check<-data.table(
  version=c("Version_A","Version_B"),
  code=c(
    image_matrix_A[[check_female]][first_Y_row,first_Y_col],
    image_matrix_B[[check_female]][first_Y_row,first_Y_col]
  ),
  RGB=c(
    paste(round(rgb_image_A[[check_female]][first_Y_row,first_Y_col,]*255),collapse=","),
    paste(round(rgb_image_B[[check_female]][first_Y_row,first_Y_col,]*255),collapse=",")
  )
)

Y_color_check

#step9把RGB矩陣輸出成PNG圖片檔----
library(png)
#建立Version A輸出資料夾；showWarnings=FALSE表示資料夾已存在時不要警告
dir.create("image_version_A_shared_NoCall",showWarnings=FALSE)
#建立Version B輸出資料夾
dir.create("image_version_B_split_NoCall",showWarnings=FALSE)

#取得RGB影像list中的完整樣本名稱
full_sample_names<-names(rgb_image_A)
#刪除_(Axiom_TPM)之後的文字，留下DM-002或EOAD_P10等名稱
short_sample_names<-sub("_\\(Axiom_TPM\\).*","",full_sample_names)
#確認完整名稱與簡短名稱的對應關係
data.table(full_sample_names,short_sample_names)
#建立第一位樣本Version A的檔案路徑
test_file_A<-file.path("image_version_A_shared_NoCall",paste0(short_sample_names[1],"_A.png"))
#輸出第一位樣本的RGB影像
writePNG(rgb_image_A[[1]],target=test_file_A)
#確認檔案是否成功建立
file.exists(test_file_A)

#建立第一位樣本Version B的檔案路徑
test_file_B<-file.path("image_version_B_split_NoCall",paste0(short_sample_names[1],"_B.png"))
#輸出第一位樣本的Version B影像
writePNG(rgb_image_B[[1]],target=test_file_B)
#確認檔案是否成功建立
file.exists(test_file_B)
#重新讀取剛輸出的Version A圖片
test_read_A<-readPNG(test_file_A)
#確認圖片尺寸為826×826×3
dim(test_read_A)
#確認讀回的RGB值範圍仍為0到1
range(test_read_A)
#比較原始RGB矩陣和輸出後重新讀取的圖片
all.equal(test_read_A,rgb_image_A[[1]],tolerance=1/255,check.attributes=FALSE)

#輸出10張Version A影像
for(i in seq_along(rgb_image_A)){
  #建立目前樣本的輸出檔名
  output_file<-file.path("image_version_A_shared_NoCall",paste0(short_sample_names[i],"_A.png"))
  #寫出PNG影像
  writePNG(rgb_image_A[[i]],target=output_file)
}

#輸出10張Version B影像
for(i in seq_along(rgb_image_B)){
  #建立目前樣本的輸出檔名
  output_file<-file.path("image_version_B_split_NoCall",paste0(short_sample_names[i],"_B.png"))
  #寫出PNG影像
  writePNG(rgb_image_B[[i]],target=output_file)
}
#列出Version A的PNG檔案
files_A<-list.files("image_version_A_shared_NoCall",pattern="\\.png$",full.names=TRUE)
#列出Version B的PNG檔案
files_B<-list.files("image_version_B_split_NoCall",pattern="\\.png$",full.names=TRUE)
#確認兩個版本的圖片數量
length(files_A)
length(files_B)


#Step10 建立並驗證 pixel 與 SNP 的固定對照表----
#建立 pixel–SNP mapping----
#在 snp_master_v2_n 新增 pixel 座標

#將每個SNP的image_order換算成影像列座標
#第1到826個SNP位於第1列，第827個SNP從第2列開始
snp_master_v2_n[,pixel_row:=((image_order-1L)%/%image_side)+1L]

#將每個SNP的image_order換算成影像欄座標
#每列的欄位由1到826循環
snp_master_v2_n[,pixel_col:=((image_order-1L)%%image_side)+1L]

#建立pixel與SNP annotation的對照表
pixel_snp_map<-snp_master_v2_n[,.(image_order,pixel_row,pixel_col,probeset_id,dbSNP_RS_ID,Chr_id,Start,Allele_A,Allele_B,distance)]

#依image_order確認對照表順序
setorder(pixel_snp_map,image_order)

#查看前10個pixel對應的SNP
pixel_snp_map[1:10]

#確定image_order與pixel座標能互換
#根據pixel_row與pixel_col重新計算image_order
pixel_snp_map[,recovered_image_order:=(pixel_row-1L)*image_side+pixel_col]

#確認重新計算的順序與原始image_order完全一致
all(pixel_snp_map$image_order==pixel_snp_map$recovered_image_order)

#check 位置都沒問題
#固定隨機種子，讓每次抽到相同的100個SNP
set.seed(123)
#隨機抽取100個image_order
check_orders<-sample(pixel_snp_map$image_order,100)
#逐一確認這100個SNP在一維序列與二維矩陣中的code完全相同
mapping_check<-sapply(check_orders,function(i){
  r<-pixel_snp_map[image_order==i,pixel_row]
  c<-pixel_snp_map[image_order==i,pixel_col]
  genotype_code_A[i,check_sample]==image_matrix_A[[check_sample]][r,c]
})
all(mapping_check)

#step6看影像
# radomly select ane male/female
#選擇第一位沒有Y染色體訊號的樣本"DM-002_(Axiom_TPM)_C03.CEL_call_code"
female_view<-female_samples[1]
#選擇第一位具有Y染色體訊號的樣本 "EOAD_P13_(Axiom_TPM)_B09.CEL_call_code"
male_view<-male_samples[1]


#將Version A的RGB矩陣直接顯示成影像
#axes=FALSE表示不顯示一般座標軸
plot(as.raster(rgb_image_A[[female_view]]),axes=FALSE,main=paste0(sub("_\\(Axiom_TPM\\).*","",female_view),"－Version A"))
#Version B將Y_not_applicable獨立顯示為黃色
plot(as.raster(rgb_image_B[[female_view]]),axes=FALSE,main=paste0(sub("_\\(Axiom_TPM\\).*","",female_view),"－Version B"))


#female:"DM-002_(Axiom_TPM)_C03.CEL_call_code"-----


#建立左右兩個繪圖區，並縮小邊界
par(mfrow=c(1,2),mar=c(1,1,3,1))
#左邊顯示Version A
plot.new()
plot.window(xlim=c(0,1),ylim=c(0,1),asp=1)
rasterImage(as.raster(rgb_image_A[[female_view]]),0,0,1,1,interpolate=FALSE)
title(main="Female:Version A")

#右邊顯示Version B
plot.new()
plot.window(xlim=c(0,1),ylim=c(0,1),asp=1)
rasterImage(as.raster(rgb_image_B[[female_view]]),0,0,1,1,interpolate=FALSE)
title(main="Female:Version B")

#恢復單張圖設定
par(mfrow=c(1,1))



#male:"EOAD_P13_(Axiom_TPM)_B09.CEL_call_code"----
#建立左右兩個繪圖區，並縮小邊界
par(mfrow=c(1,2),mar=c(1,1,3,1))
#左邊顯示Version A
plot.new()
plot.window(xlim=c(0,1),ylim=c(0,1),asp=1)
rasterImage(as.raster(rgb_image_A[[male_view]]),0,0,1,1,interpolate=FALSE)
title(main="Male:Version A")

#右邊顯示Version B
plot.new()
plot.window(xlim=c(0,1),ylim=c(0,1),asp=1)
rasterImage(as.raster(rgb_image_B[[male_view]]),0,0,1,1,interpolate=FALSE)
title(main="Male:Version B")

#恢復單張圖設定
par(mfrow=c(1,1))

#10個人一起看----
#version a:
#將10位樣本排成2列×5欄，不另外標示性別
par(mfrow=c(2,5),mar=c(1,1,2,1))
for(x in sample_cols){
  #建立空白正方形繪圖區
  plot.new()
  plot.window(xlim=c(0,1),ylim=c(0,1),asp=1)
  #顯示目前樣本的Version A影像
  rasterImage(as.raster(rgb_image_A[[x]]),0,0,1,1,interpolate=FALSE)
  #標題只顯示簡短樣本名稱
  title(main=sub("_\\(Axiom_TPM\\).*","",x),cex.main=0.8)
}
#恢復單張圖設定
par(mfrow=c(1,1))

#version b:
#顯示Version B的10位樣本
par(mfrow=c(2,5),mar=c(1,1,2,1))
for(x in sample_cols){
  plot.new()
  plot.window(xlim=c(0,1),ylim=c(0,1),asp=1)
  rasterImage(as.raster(rgb_image_B[[x]]),0,0,1,1,interpolate=FALSE)
  title(main=sub("_\\(Axiom_TPM\\).*","",x),cex.main=0.8)
}
par(mfrow=c(1,1))
