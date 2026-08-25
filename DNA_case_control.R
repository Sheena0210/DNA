#150人的case control(1:2)


#--------------------------------------------#



install.packages("pgenlibr")
library("pgenlibr")
library("data.table")
#ReadIntList() 回傳的矩陣是「row = sample、column = variant」，值為 0/1/2/NA 的 ALT allele hardcall 數

#分每一條染色體讀入


pgen<-file.choose() #選.pgen
plink_prefix<-sub("\\.pgen$","",pgen) #抓出路徑

#找出三個檔案的路徑
pgen_path<-paste0(plink_prefix,".pgen")
pvar_path<-paste0(plink_prefix,".pvar")
psam_path<-paste0(plink_prefix,".psam")

#.psam 是受試者資料
psam_dt<-fread(psam_path)
#.pvar 是 SNP/variant annotation
pvar_dt<-fread(pvar_path)

head(psam)
head(pvar)
dim(pvar_dt)

#改欄位名稱
##CHROM -> Chr_id
if("#CHROM"%in%names(pvar_dt)){setnames(pvar_dt,"#CHROM","Chr_id")}
setnames(pvar_dt,c("POS","ID"),c("Start","variant_id"))
#非常重要：記錄每顆variant在pgen中的原始index
pvar_dt[,pgen_index:=.I]

#pgen
#把 talent_150.pvar 這個 variant 資訊檔，用 pgenlibr 打開，建立一個 R 可以拿來讀取的物件
pvar_reader<-pgenlibr::NewPvar(pvar_path)
#open pgen
pgen_reader<-pgenlibr::NewPgen(pgen_path,pvar = pvar_reader)

pgenlibr::GetRawSampleCt(pgen_reader) #n=150
pgenlibr::GetVariantCt(pgen_reader) #snp=684603

#chr1----
chr1<-pvar_dt[Chr_id=="1",pgen_index]
head(pvar)
length(chr1)
head(chr1)

geno_chr1<-pgenlibr::ReadIntList(pgen_reader,chr1)
dim(geno_chr1)
table(geno_chr1,useNA = "ifany")
#加入psam的資料
rownames(geno_chr1)<-psam$IID
geno_chr1[1:5,1:10]

#pvar QC:
#[1] "Chr_id""Start"   "variant_id" "REF"   "ALT"    "pgen_index"
head(pvar_dt)
names(pvar_dt)
#"Start"轉數值
pvar_dt[,Start:=as.numeric(Start)]
#"Chr_id排除y染色體
chr_keep<-c(as.character(1:22),"X","MT")
#先看原始chromosome分布
pvar_dt[,.(variant_n=.N),by=Chr_id][
  order(match(Chr_id,c(chr_keep,"Y")))
]
#排除y
pvar_dt[,standard_chr:=Chr_id%in%chr_keep]

#biallelic：ALT不能有多個alleles
pvar_dt[,biallelic:=!grepl(",",ALT)]

#SNV：REF與ALT都只有1個base
pvar_dt[,is_snv:=nchar(REF)==1&nchar(ALT)==1]

#只接受A/C/G/T
pvar_dt[,is_acgt:=REF%in%c("A","C","G","T")&ALT%in%c("A","C","G","T")]

pvar_qc_summary<-data.table(
  total_variant=nrow(pvar_dt),
  standard_chr=sum(pvar_dt$standard_chr),
  biallelic=sum(pvar_dt$biallelic),
  SNV=sum(pvar_dt$is_snv),
  ACGT=sum(pvar_dt$is_acgt),
  final_keep=sum(
    pvar_dt$standard_chr&
      pvar_dt$biallelic&
      pvar_dt$is_snv&
      pvar_dt$is_acgt
  )
)

pvar_qc_summary

#真正要進影像的 SNP annotation
pvar_snv<-pvar_dt[
  standard_chr==TRUE&
    biallelic==TRUE&
    is_snv==TRUE&
    is_acgt==TRUE
]

#依染色體與位置排序
pvar_snv[,chr_order:=match(Chr_id,chr_keep)]
setorder(pvar_snv,chr_order,Start,pgen_index)
#確認
dim(pvar_snv)
head(pvar_snv)

pvar_snv[
  ,.(SNP_n=.N),
  by=Chr_id
][order(match(Chr_id,chr_keep))]

#看是否有同一 chromosome + position 多個 variants
duplicate_position<-pvar_snv[
  ,N:=.N,
  by=.(Chr_id,Start)
][N>1]

nrow(duplicate_position)

head(
  duplicate_position[
    order(Chr_id,Start)
  ]
)

anyDuplicated(pvar_snv$pgen_index)

range(pvar_snv$pgen_index)

#抓Chr1 SNP index
pvar_chr1<-copy(pvar_snv[Chr_id=="1"])
chr1_index<-pvar_chr1$pgen_index

length(chr1_index)
geno_chr1<-pgenlibr::ReadIntList(
  pgen_reader,
  chr1_index
)

dim(geno_chr1)

table(
  geno_chr1,
  useNA="ifany"
)

if("#IID"%in%names(psam_dt)){
  setnames(psam_dt,"#IID","IID")
}

rownames(geno_chr1)<-psam_dt$IID


#0/1/2 分別是 ALT allele copies，NA 是 missing
#encoding:
#1 = AA / TT
#2 = CC / GG
#3 = AG / CT
#4 = AC / GT
#5 = AT
#6 = CG
#7 = NoCall

#先建立 Chr1 的 REF/ALT → code 對照
pvar_chr1<-copy(pvar_snv[Chr_id=="1"])
#Chr1在pgen中的原始index
chr1_index<-pvar_chr1$pgen_index

#讀取150人的Chr1 genotype
geno_chr1<-pgenlibr::ReadIntList(
  pgen_reader,
  chr1_index
)

dim(geno_chr1)

#REF與ALT
ref_chr1<-pvar_chr1$REF
alt_chr1<-pvar_chr1$ALT

#heterozygous genotype統一成不分順序的pair
pair_chr1<-ifelse(
  ref_chr1<alt_chr1,
  paste0(ref_chr1,alt_chr1),
  paste0(alt_chr1,ref_chr1)
)

#pgen=0：REF/REF
code_0_chr1<-fifelse(
  ref_chr1%in%c("A","T"),
  1L,
  2L
)

#pgen=1：REF/ALT
code_1_chr1<-fcase(
  pair_chr1%in%c("AG","CT"),3L,
  pair_chr1%in%c("AC","GT"),4L,
  pair_chr1=="AT",5L,
  pair_chr1=="CG",6L,
  default=NA_integer_
)

#pgen=2：ALT/ALT
code_2_chr1<-fifelse(
  alt_chr1%in%c("A","T"),
  1L,
  2L
)

code_lookup_chr1<-rbind(
  code_0_chr1,
  code_1_chr1,
  code_2_chr1
)

#先全部設為7 = NoCall
geno_code_chr1<-matrix(
  7L,
  nrow=nrow(geno_chr1),
  ncol=ncol(geno_chr1)
)

#找非missing genotype位置
valid_idx_chr1<-which(
  !is.na(geno_chr1),
  arr.ind=TRUE
)

#依0/1/2轉成Version 7的1–6 code
geno_code_chr1[valid_idx_chr1]<-
  code_lookup_chr1[
    cbind(
      geno_chr1[valid_idx_chr1]+1L,
      valid_idx_chr1[,2]
    )
  ]
table(geno_code_chr1,useNA="ifany")

names(psam_dt)
rownames(geno_code_chr1)<-psam_dt$IID
geno_code_chr1[1:5,1:10]

#記錄Chr1 variant在geno_code_chr1中的column
pvar_chr1[,variant_col_chr1:=.I]
#20% overlapping sliding window
window_size_bp_v7<-1000000L
overlap_ratio_v7<-0.2
overlap_bp_v7<-as.integer(window_size_bp_v7*overlap_ratio_v7) #200,000
stride_bp_v7<-window_size_bp_v7-overlap_bp_v7   


#Chr1最大位置
chr1_max_start<-max(pvar_chr1$Start,na.rm=TRUE)

#每800 kb建立一個1 Mb candidate window
window_start_chr1<-seq(
  1,
  chr1_max_start,
  by=stride_bp_v7
)

window_map_chr1_v7<-data.table(
  Chr_id="1",
  window_id_v7=seq_along(window_start_chr1),
  window_start_v7=window_start_chr1,
  window_end_v7=window_start_chr1+window_size_bp_v7-1L
)

#window center
window_map_chr1_v7[,window_center_v7:=
                     window_start_v7+(window_size_bp_v7-1)/2]

head(window_map_chr1_v7)

#SNP位置轉成interval
snp_interval_chr1_v7<-pvar_chr1[,.
                                (
                                  variant_col_chr1,
                                  pgen_index,
                                  variant_id,
                                  Chr_id,
                                  Start,
                                  REF,
                                  ALT,
                                  snp_start_v7=Start,
                                  snp_end_v7=Start
                                )
]

#foverlaps
setkey(
  window_map_chr1_v7,
  Chr_id,
  window_start_v7,
  window_end_v7
)

snp_window_chr1_v7<-foverlaps(
  snp_interval_chr1_v7,
  window_map_chr1_v7,
  by.x=c("Chr_id","snp_start_v7","snp_end_v7"),
  by.y=c("Chr_id","window_start_v7","window_end_v7"),
  type="within",
  nomatch=NULL
)

#看每顆 SNP 最多落進幾個 windows
membership_chr1_v7<-snp_window_chr1_v7[
  ,.(window_membership_n=.N),
  by=variant_col_chr1
]

table(membership_chr1_v7$window_membership_n)
max(membership_chr1_v7$window_membership_n)
#計算SNP與candidate window center的距離
snp_window_chr1_v7[,distance_to_center_v7:=
                     abs(snp_start_v7-window_center_v7)]

#距離最近者優先；tie時保留前一個window
setorder(
  snp_window_chr1_v7,
  variant_col_chr1,
  distance_to_center_v7,
  window_id_v7
)

#每顆SNP只保留一個owner window
snp_owner_chr1_v7<-snp_window_chr1_v7[
  ,.SD[1L],
  by=variant_col_chr1
]

#QC：每顆Chr1 SNP應只出現一次
nrow(snp_owner_chr1_v7)
nrow(pvar_chr1)

nrow(snp_owner_chr1_v7)==nrow(pvar_chr1)

#排進 64×56 serpentine patch

patch_width_v7<-64L
patch_height_v7<-56L
patch_capacity_v7<-patch_width_v7*patch_height_v7 #3584

#先看unique ownership後，每個window最多有幾顆SNP
owner_count_chr1_v7<-snp_owner_chr1_v7[
  ,.(SNP_n=.N),
  by=window_id_v7
]

summary(owner_count_chr1_v7$SNP_n)
max(owner_count_chr1_v7$SNP_n)



#每個owner window內依Start排序
setorder(
  snp_owner_chr1_v7,
  window_id_v7,
  Start,
  pgen_index
)

#window內SNP順序
snp_owner_chr1_v7[,snp_order_v7:=
                    seq_len(.N),
                  by=window_id_v7]

#row
snp_owner_chr1_v7[,local_row_v7:=
                    ((snp_order_v7-1L)%/%patch_width_v7)+1L]

#row內位置
snp_owner_chr1_v7[,position_in_row_v7:=
                    ((snp_order_v7-1L)%%patch_width_v7)+1L]

#蛇行column
snp_owner_chr1_v7[,local_col_v7:=
                    fifelse(
                      local_row_v7%%2L==1L,
                      position_in_row_v7,
                      patch_width_v7-position_in_row_v7+1L
                    )]


#Step 15：Chr1建立完整rectangular strip座標----------------

#Chr1總window數
chr1_window_n_v7<-nrow(window_map_chr1_v7)

#整條Chr1 strip寬度
chr1_strip_width_v7<-chr1_window_n_v7*patch_width_v7


#window內座標 → chromosome-level座標
snp_owner_chr1_v7[,global_row_v7:=local_row_v7]

snp_owner_chr1_v7[,global_col_v7:=
                    (window_id_v7-1L)*patch_width_v7+
                    local_col_v7]

#QC
range(snp_owner_chr1_v7$global_row_v7)
range(snp_owner_chr1_v7$global_col_v7)

chr1_strip_width_v7



sample_index_test<-1L
sample_id_test<-rownames(geno_code_chr1)[sample_index_test]

sample_id_test

#0 = padding
chr1_strip_test<-matrix(
  0L,
  nrow=patch_height_v7,
  ncol=chr1_strip_width_v7
)
#依variant_col_chr1抓對應genotype
sample_code_chr1<-
  geno_code_chr1[
    sample_index_test,
    snp_owner_chr1_v7$variant_col_chr1
  ]
chr1_strip_test[
  cbind(
    snp_owner_chr1_v7$global_row_v7,
    snp_owner_chr1_v7$global_col_v7
  )
]<-sample_code_chr1

#Version 7 獨立顏色表：名稱要和code完全對應
genotype_color_v7<-c(
  "white",     #0 Padding
  "#E41A1C",   #1 AA/TT
  "#377EB8",   #2 CC/GG
  "#4DAF4A",   #3 AG/CT
  "#984EA3",   #4 AC/GT
  "#FF7F00",   #5 AT
  "#A65628",   #6 CG
  "grey40"     #7 NoCall
)


#Step 18：Chr1 QC image----------------

chr1_test_file<-paste0(
  sample_id_test,
  "_Chr1_Version7_test.png"
)

png(
  chr1_test_file,
  width=16000,
  height=1500,
  res=150
)

par(
  mar=c(5,5,3,1),
  xaxs="i",
  yaxs="i"
)

image(
  x=seq_len(ncol(chr1_strip_test)),
  y=1:patch_height_v7,
  z=t(chr1_strip_test),
  col=genotype_color_v7,
  breaks=seq(-0.5,7.5,by=1),
  axes=FALSE,
  xlab="",
  ylab="",
  useRaster=TRUE
)

#每64 columns = 一個sliding window
window_boundary_chr1<-seq(
  patch_width_v7+0.5,
  chr1_strip_width_v7-0.5,
  by=patch_width_v7
)

abline(
  v=window_boundary_chr1,
  col="grey80",
  lwd=0.3
)

#Y軸：row 1在下，56在上
axis(
  side=2,
  at=c(1,14,28,42,56),
  labels=c(1,14,28,42,56),
  las=1
)

mtext(
  "Snake row",
  side=2,
  line=3
)

title(
  main=paste0(
    sample_id_test,
    " - Chr1 Version 7 Genotype Strip"
  )
)

box()

dev.off()

system(
  paste(
    "open",
    shQuote(chr1_test_file)
  )
)

#輸出全部150人的Chr1 Version 7 strips----------------

#若psam的sample ID還沒設成rownames，先設定
if(is.null(rownames(geno_code_chr1))){
  if("#IID"%in%names(psam_dt)){setnames(psam_dt,"#IID","IID")}
  rownames(geno_code_chr1)<-psam_dt$IID
}

#建立輸出資料夾
output_dir_chr1_v7<-file.path(getwd(),"V7_Chr1_150samples")
dir.create(output_dir_chr1_v7,showWarnings=FALSE)

#記錄QC結果
chr1_qc_v7<-vector("list",nrow(geno_code_chr1))

#Chr1共有幾個window
chr1_window_n_v7<-nrow(window_map_chr1_v7)

#每10 Mb標線位置（以window patch邊界估計）
major_window_chr1_v7<-seq(10L,chr1_window_n_v7,by=10L)
major_x_chr1_v7<-major_window_chr1_v7*patch_width_v7+0.5

#每個window邊界
window_boundary_chr1_v7<-seq(patch_width_v7+0.5,chr1_strip_width_v7-0.5,by=patch_width_v7)

for(i in seq_len(nrow(geno_code_chr1))){
  
  #sample ID
  sample_id_now<-rownames(geno_code_chr1)[i]
  
  #建立空白矩陣；0=padding
  chr1_strip_now<-matrix(
    0L,
    nrow=patch_height_v7,
    ncol=chr1_strip_width_v7
  )
  
  #抓這位sample在Chr1、且依owner scaffold排序後的genotype code
  sample_code_now<-geno_code_chr1[
    i,
    snp_owner_chr1_v7$variant_col_chr1
  ]
  
  #填入固定位置
  chr1_strip_now[
    cbind(
      snp_owner_chr1_v7$global_row_v7,
      snp_owner_chr1_v7$global_col_v7
    )
  ]<-sample_code_now
  
  #QC
  nonpadding_n_now<-sum(chr1_strip_now!=0L)
  expected_n_now<-nrow(snp_owner_chr1_v7)
  match_qc_now<-nonpadding_n_now==expected_n_now
  has_na_now<-anyNA(chr1_strip_now)
  
  #輸出檔名
  output_file_now<-file.path(
    output_dir_chr1_v7,
    paste0(sample_id_now,"_Chr1_Version7_strip.png")
  )
  
  #依圖寬調整png寬度
  png_width_now<-max(3000,min(18000,round(ncol(chr1_strip_now)*2)))
  
  png(
    filename=output_file_now,
    width=png_width_now,
    height=1400,
    res=150
  )
  
  par(mar=c(5,5,3,1),xaxs="i",yaxs="i")
  
  image(
    x=seq_len(ncol(chr1_strip_now)),
    y=1:patch_height_v7,
    z=t(chr1_strip_now),
    col=genotype_color_v7,
    breaks=seq(-0.5,7.5,by=1),
    axes=FALSE,
    xlab="",
    ylab="",
    useRaster=TRUE
  )
  
  #每個window邊界
  if(length(window_boundary_chr1_v7)>0){
    abline(v=window_boundary_chr1_v7,col="grey85",lwd=0.3)
  }
  
  #每10 Mb粗線
  if(length(major_x_chr1_v7)>0){
    abline(v=major_x_chr1_v7,col="black",lwd=1)
    axis(side=1,at=major_x_chr1_v7,labels=paste0(major_window_chr1_v7," Mb"),tick=FALSE,cex.axis=0.7)
  }
  
  axis(side=2,at=c(1,14,28,42,56),labels=c(1,14,28,42,56),las=1,cex.axis=0.8)
  mtext("Genomic position",side=1,line=3)
  mtext("Snake row",side=2,line=3)
  
  title(main=paste0(sample_id_now," - Chr1 Version 7 Genotype Strip"))
  
  box()
  dev.off()
  
  #儲存QC
  chr1_qc_v7[[i]]<-data.table(
    sample_id=sample_id_now,
    SNP_n_expected=expected_n_now,
    nonpadding_n=nonpadding_n_now,
    match_qc=match_qc_now,
    has_na=has_na_now,
    output_file=output_file_now
  )
  
  #釋放這一輪矩陣
  rm(chr1_strip_now,sample_code_now)
  gc()
}

#整理QC結果
chr1_qc_v7<-rbindlist(chr1_qc_v7,use.names=TRUE,fill=TRUE)

#查看QC
chr1_qc_v7
all(chr1_qc_v7$match_qc)
any(chr1_qc_v7$has_na)

#打開輸出資料夾
system(paste("open",shQuote(output_dir_chr1_v7)))


#---------------------------------------------#
# data:sample information.xlsx
#先進行欄位轉置
library(data.table)
install.packages("readxl")
library(readxl)
sample_information<-file.choose()
sample_info<-as.data.table(read_excel(sample_information,sheet="sample_list"))

#step1 新增group_id (一組有一個case 兩個cotrol)
#轉置
sample_info[,group_id:=.I]
head(sample_info[,.(group_id,case_id,control_id1,control_id2)])

#step2 分割case
#新增id=case_id
case<-sample_info[,.(id=case_id,group_id,candi,Age_lynn,GENDER,SMOKE,AGE)]
#新增sample_type=case
case[,sample_type:="case"]

#step3 分割control1
control1<-sample_info[,.(id=control_id1,group_id, candi=candi.ctrl1,Age_lynn=Age_lynn.ctrl1,
                         GENDER=GENDER.ctrl1,SMOKE=SMOKE.ctrl1,AGE=AGE.ctrl1)]
#新增sample_type=control1
control1[,sample_type:="control1"]


#step4 分割control2
control2<-sample_info[,.(id=control_id2,group_id, candi=candi.ctrl2,Age_lynn=Age_lynn.ctrl2,
                         GENDER=GENDER.ctrl2,SMOKE=SMOKE.ctrl2,AGE=AGE.ctrl2)]
#新增sample_type=control2
control2[,sample_type:="control2"]

#step5 合併
sample<-rbindlist(list(case, control1, control2),use.names = TRUE)
dim(sample)

#按照group排序
setorder(sample,group_id,-sample_type)
head(sample)

#找出最年輕的case 
sample[,case_control:=
         fifelse(sample_type=="case",
                 "case",
                 "control")]
#min:49 max63
sample[case_control=="case",range(AGE,na.rm = TRUE)] #NA remove
youngest_case<-sample[
  case_control=="case"
][
  AGE==min(AGE,na.rm=TRUE)
]



#找出最老的control 共有14個 但跟case資訊一樣的有12個

sample[case_control=="control",range(AGE,na.rm = TRUE)] #NA remove

oldest_control<-sample[
  case_control=="control"
][
  AGE==max(AGE,na.rm=TRUE)
]



#Table 1：case vs control demographic characteristics----
#表示matching結構是相同的：同一個 group_id 裡面只有一種 AGE，也就是三個人年齡相同
match_check<-sample[
  ,.(age_same=uniqueN(AGE)==1,
     gender_same=uniqueN(GENDER)==1,
     smoke_same=uniqueN(SMOKE)==1),
  by=group_id
]

head(match_check)
#年齡分佈
age_summary<-sample[
  ,.(N=.N,
     Mean=mean(AGE,na.rm=TRUE),
     SD=sd(AGE,na.rm=TRUE),
     Median=median(AGE,na.rm=TRUE),
     Min=min(AGE,na.rm=TRUE),
     Max=max(AGE,na.rm=TRUE)),
  by=case_control
]

age_summary
#Gender
sample[,Gender_label:=
         fifelse(GENDER==1,"Male",
                 fifelse(GENDER==2,"Female",NA_character_))]
gender_summary<-sample[
  ,.(N=.N),
  by=.(case_control,Gender_label)
]

gender_summary[
  ,Percent:=round(N/sum(N)*100,2),
  by=case_control
]

gender_summary
#smoking
sample[,Smoke_label:=
         fifelse(SMOKE==1,"Smoker",
                 fifelse(SMOKE==2,"Non-smoker",NA_character_))]
smoke_summary<-sample[
  ,.(N=.N),
  by=.(case_control,Smoke_label)
]

smoke_summary[
  ,Percent:=round(N/sum(N)*100,2),
  by=case_control
]

smoke_summary

#-------------------------------------------#

#sample_information跟chip_info的Id名稱沒有統一
#建立標準的id:match_id------------------------
#建立matching用ID
sample[,match_id:=sub("^CT-","CT",id)]
head(sample[,.(id,match_id)])

#psam 也建立相同規則的 ID
psam_dt[,match_id:=sub("^CT-","CT",IID)]
head(psam_dt[,.(IID,match_id)])

sample$match_id%in%psam_dt$match_id


#youngest case：CT-06-01508
#oldest control：CT-01-00418
target_id<-c(
  "CT-06-01508",
  "CT-01-00418"
)
#用match_id抓
target_match_id<-sample[
  match(target_id,id),
  match_id
]

target_match_id

target_index<-match(
  target_match_id,
  psam_dt$match_id
)

target_index

target_check<-data.table(
  id=target_id,
  match_id=target_match_id,
  pgen_sample_index=target_index,
  psam_IID=psam_dt$IID[target_index]
)

target_check


#case-control difference image
difference_chr1<-
  case_strip!=control_strip
table(difference_chr1)

#Step 8：建立case-control difference matrix----------------

#相同=0，不同=1
diff_strip_chr1<-(case_strip!=control_strip)*1L

#查看0和1的數量
table(diff_strip_chr1,useNA="ifany")

#只看真實SNP位置中的差異數量，不把padding算進去
real_snp_mask_chr1<-(case_strip!=0L)&(control_strip!=0L)

diff_summary_chr1<-data.table(
  total_real_SNP=sum(real_snp_mask_chr1),
  same_genotype=sum(diff_strip_chr1[real_snp_mask_chr1]==0L),
  different_genotype=sum(diff_strip_chr1[real_snp_mask_chr1]==1L)
)

#加上百分比
diff_summary_chr1[,different_percent:=round(different_genotype/total_real_SNP*100,2)]
diff_summary_chr1[,same_percent:=round(same_genotype/total_real_SNP*100,2)]

diff_summary_chr1

#Step 9：padding / same / different 三類版本----------------

diff3_strip_chr1<-matrix(
  0L,
  nrow=nrow(case_strip),
  ncol=ncol(case_strip)
)

#真實SNP位置且相同
diff3_strip_chr1[real_snp_mask_chr1&(case_strip==control_strip)]<-1L

#真實SNP位置且不同
diff3_strip_chr1[real_snp_mask_chr1&(case_strip!=control_strip)]<-2L

table(diff3_strip_chr1,useNA="ifany")


diff3_file_chr1<-file.path(
  output_dir_target_chr1,
  "YoungestCase_vs_OldestControl_Chr1_difference_with_padding.png"
)

png(
  diff3_file_chr1,
  width=6000,
  height=1200,
  res=150
)

par(
  mar=c(4,5,3,1),
  xaxs="i",
  yaxs="i"
)

image(
  x=seq_len(ncol(diff3_strip_chr1)),
  y=1:nrow(diff3_strip_chr1),
  z=t(diff3_strip_chr1),
  col=c("grey85","white","black"),
  breaks=c(-0.5,0.5,1.5,2.5),
  axes=FALSE,
  xlab="",
  ylab="",
  useRaster=TRUE
)

axis(
  side=2,
  at=c(1,14,28,42,56),
  labels=c(1,14,28,42,56),
  las=1
)

mtext(
  "Snake row",
  side=2,
  line=3
)

mtext(
  "Genomic order",
  side=1,
  line=2
)

title(
  main="Chr1 Difference Map with Padding: Youngest Case vs Oldest Control"
)

box()

legend(
  "topright",
  legend=c("Padding","Same genotype","Different genotype"),
  fill=c("grey85","white","black"),
  bty="n",
  cex=0.9
)

dev.off()

system(
  paste(
    "open",
    shQuote(diff3_file_chr1)
  )
)

#24條染色體----

#為了後續模型學習->Patch 64×48:3072

chr_levels_v9<-c(as.character(1:22),"X","MT")
#1 Mb sliding window 目前還是20%overlap
window_size_bp_v9<-1000000L
overlap_ratio_v9<-0.2
overlap_bp_v9<-as.integer(window_size_bp_v9*overlap_ratio_v9) #200 kb
stride_bp_v9<-window_size_bp_v9-overlap_bp_v9                 #800 kb


patch_width_v9<-64L
patch_height_v9<-48L
patch_capacity_v9<-patch_width_v9*patch_height_v9
patch_capacity_v9

##沿用已經讀好的24條染色體資料
pvar_chr_list_v9<-lapply(
  pvar_chr_list_v7,
  copy
)

geno_code_list_v9<-geno_code_list_v7

names(pvar_chr_list_v9)
names(geno_code_list_v9)

#Step 3：建立 V9 scaffold function

build_chr_scaffold_v9<-function(pvar_chr_now){
  
  #依實體位置排序
  setorder(pvar_chr_now,Start,pgen_index)
  
  #這條染色體
  chr_now<-unique(pvar_chr_now$Chr_id)
  
  #重新記錄這條染色體在genotype matrix中的column
  pvar_chr_now[,variant_col_chr:=.I]
  
  #染色體最後一顆SNP位置
  chr_max_now<-max(
    pvar_chr_now$Start,
    na.rm=TRUE
  )
  
  
  #-------------------------
  #1.建立1 Mb sliding windows
  
  window_start_now<-seq(
    1,
    chr_max_now,
    by=stride_bp_v9
  )
  
  window_map_now<-data.table(
    Chr_id=chr_now,
    window_id_v9=seq_along(window_start_now),
    window_start_v9=window_start_now,
    window_end_v9=window_start_now+
      window_size_bp_v9-1L
  )
  
  #window中心
  window_map_now[
    ,window_center_v9:=
      window_start_v9+
      (window_size_bp_v9-1)/2
  ]
  
  
  #-------------------------
  #2.SNP轉成point interval
  
  snp_interval_now<-pvar_chr_now[,.
                                 (
                                   variant_col_chr,
                                   pgen_index,
                                   variant_id,
                                   Chr_id,
                                   Start,
                                   REF,
                                   ALT,
                                   snp_start_v9=Start,
                                   snp_end_v9=Start
                                 )
  ]
  
  
  #-------------------------
  #3.找每顆SNP落在哪些candidate windows
  
  setkey(
    window_map_now,
    Chr_id,
    window_start_v9,
    window_end_v9
  )
  
  snp_window_now<-foverlaps(
    snp_interval_now,
    window_map_now,
    by.x=c(
      "Chr_id",
      "snp_start_v9",
      "snp_end_v9"
    ),
    by.y=c(
      "Chr_id",
      "window_start_v9",
      "window_end_v9"
    ),
    type="within",
    nomatch=NULL
  )
  
  
  #-------------------------
  #4.nearest-center
  
  snp_window_now[
    ,distance_to_center_v9:=
      abs(Start-window_center_v9)
  ]
  
  #距離最近優先
  #若距離完全相同，選前面的window
  setorder(
    snp_window_now,
    variant_col_chr,
    distance_to_center_v9,
    window_id_v9
  )
  
  #每顆SNP只留一個owner window
  owner_now<-snp_window_now[
    ,.SD[1L],
    by=variant_col_chr
  ]
  
  
  #-------------------------
  #5.window內依實體位置排序
  
  setorder(
    owner_now,
    window_id_v9,
    Start,
    pgen_index
  )
  
  owner_now[
    ,snp_order_v9:=seq_len(.N),
    by=window_id_v9
  ]
  
  
  #-------------------------
  #6.64×48 serpentine coordinates
  
  #row
  owner_now[
    ,local_row_v9:=
      ((snp_order_v9-1L)%/%patch_width_v9)+1L
  ]
  
  #row內的位置
  owner_now[
    ,position_in_row_v9:=
      ((snp_order_v9-1L)%%patch_width_v9)+1L
  ]
  
  #蛇行column
  owner_now[
    ,local_col_v9:=
      fifelse(
        local_row_v9%%2L==1L,
        position_in_row_v9,
        patch_width_v9-position_in_row_v9+1L
      )
  ]
  
  
  #-------------------------
  #7.chromosome-level global coordinates
  
  owner_now[
    ,global_row_v9:=local_row_v9
  ]
  
  owner_now[
    ,global_col_v9:=
      (window_id_v9-1L)*patch_width_v9+
      local_col_v9
  ]
  
  
  #回傳
  return(
    list(
      pvar=pvar_chr_now,
      window_map=window_map_now,
      owner_map=owner_now
    )
  )
}


#建立空list
pvar_chr_list_v9_new<-list()
window_map_list_v9<-list()
owner_map_list_v9<-list()

for(chr_now in chr_levels_v9){
  
  cat("Building Chr",chr_now,"...\n")
  
  scaffold_now<-build_chr_scaffold_v9(
    copy(pvar_chr_list_v9[[chr_now]])
  )
  
  pvar_chr_list_v9_new[[chr_now]]<-
    scaffold_now$pvar
  
  window_map_list_v9[[chr_now]]<-
    scaffold_now$window_map
  
  owner_map_list_v9[[chr_now]]<-
    scaffold_now$owner_map
  
  rm(scaffold_now)
  gc()
}

#正式取代
pvar_chr_list_v9<-pvar_chr_list_v9_new
rm(pvar_chr_list_v9_new)



#Step 5：檢查 64×48 是否全部裝得下
#V9 capacity QC----------------

owner_capacity_qc_v9<-rbindlist(
  lapply(chr_levels_v9,function(chr_now){
    
    owner_now<-owner_map_list_v9[[chr_now]]
    
    window_count_now<-owner_now[
      ,.(SNP_n=.N),
      by=window_id_v9
    ]
    
    data.table(
      Chr_id=chr_now,
      window_n=nrow(window_map_list_v9[[chr_now]]),
      max_SNP_per_window=max(window_count_now$SNP_n),
      capacity=patch_capacity_v9,
      fit=max(window_count_now$SNP_n)<=patch_capacity_v9
    )
  })
)

owner_capacity_qc_v9[
  order(-max_SNP_per_window)
]


#Version 9 genotype colors----------------

genotype_color_v9<-c(
  "white",     #0 Padding
  "#E41A1C",   #1 AA/TT
  "#377EB8",   #2 CC/GG
  "#4DAF4A",   #3 AG/CT
  "#984EA3",   #4 AC/GT
  "#FF7F00",   #5 AT
  "#A65628",   #6 CG
  "grey40"     #7 NoCall
)

#Step 9：指定要比較的兩個人
#target samples----------------

target_id<-c(
  "CT-06-01508",
  "CT-01-00418"
)

target_match_id<-sample[
  match(target_id,id),
  match_id
]

target_index<-match(
  target_match_id,
  psam_dt$match_id
)

target_check<-data.table(
  id=target_id,
  match_id=target_match_id,
  pgen_sample_index=target_index,
  psam_IID=psam_dt$IID[target_index]
)

target_check

#Step 10：建立兩人的 24 條 V9 chromosome strips
#建立target strips----------------

target_strips_v9<-list()

for(chr_now in chr_levels_v9){
  
  cat("Creating Chr",chr_now,"...\n")
  
  geno_now<-geno_code_list_v9[[chr_now]]
  owner_now<-owner_map_list_v9[[chr_now]]
  window_now<-window_map_list_v9[[chr_now]]
  
  #這條染色體總寬度
  strip_width_now<-
    nrow(window_now)*patch_width_v9
  
  #兩位sample
  target_chr_now<-list()
  
  
  for(i in seq_along(target_index)){
    
    sample_id_now<-target_id[i]
    
    #空白matrix
    strip_now<-matrix(
      0L,
      nrow=patch_height_v9,
      ncol=strip_width_now
    )
    
    #這個人的genotype
    sample_code_now<-geno_now[
      target_index[i],
      owner_now$variant_col_chr
    ]
    
    #放進固定V9座標
    strip_now[
      cbind(
        owner_now$global_row_v9,
        owner_now$global_col_v9
      )
    ]<-sample_code_now
    
    target_chr_now[[sample_id_now]]<-strip_now
  }
  
  target_strips_v9[[chr_now]]<-target_chr_now
  
  gc()
}
#Step 12：建立 24 張 Case vs Control 圖
#上 = youngest case
#下 = oldest control

#V9 24 chromosomes output----------------

output_dir_v9<-file.path(
  getwd(),
  "V9_64x48_case_control"
)

dir.create(
  output_dir_v9,
  showWarnings=FALSE
)

for(chr_now in chr_levels_v9){
  
  cat("Plotting Chr",chr_now,"...\n")
  
  case_strip_now<-
    target_strips_v9[[chr_now]][[target_id[1]]]
  
  control_strip_now<-
    target_strips_v9[[chr_now]][[target_id[2]]]
  
  window_n_now<-
    nrow(window_map_list_v9[[chr_now]])
  
  strip_width_now<-
    ncol(case_strip_now)
  
  
  #輸出檔名
  output_file_now<-file.path(
    output_dir_v9,
    paste0(
      "Chr",chr_now,
      "_YoungestCase_vs_OldestControl_V9.png"
    )
  )
  
  
  #依chromosome長度調整圖片寬度
  png_width_now<-max(
    2500,
    min(
      18000,
      round(strip_width_now*1.8)
    )
  )
  
  
  png(
    output_file_now,
    width=png_width_now,
    height=2200,
    res=150
  )
  
  par(
    mfrow=c(2,1),
    mar=c(3,5,3,1),
    xaxs="i",
    yaxs="i"
  )
  
  
  #========================
  #Youngest case
  
  image(
    x=seq_len(ncol(case_strip_now)),
    y=1:patch_height_v9,
    z=t(case_strip_now),
    col=genotype_color_v9,
    breaks=seq(-0.5,7.5,by=1),
    axes=FALSE,
    xlab="",
    ylab="",
    useRaster=TRUE
  )
  
  #window boundaries
  if(window_n_now>=2L){
    
    window_boundary_now<-seq(
      patch_width_v9+0.5,
      strip_width_now-0.5,
      by=patch_width_v9
    )
    
    abline(
      v=window_boundary_now,
      col="grey85",
      lwd=0.25
    )
  }
  
  axis(
    side=2,
    at=c(1,12,24,36,48),
    labels=c(1,12,24,36,48),
    las=1
  )
  
  mtext(
    "Snake row",
    side=2,
    line=3
  )
  
  title(
    main=paste0(
      "Chr",chr_now,
      " - ",
      target_id[1],
      " (Youngest Case, Age 49)"
    )
  )
  
  box()
  
  
  #========================
  #Oldest control
  
  image(
    x=seq_len(ncol(control_strip_now)),
    y=1:patch_height_v9,
    z=t(control_strip_now),
    col=genotype_color_v9,
    breaks=seq(-0.5,7.5,by=1),
    axes=FALSE,
    xlab="",
    ylab="",
    useRaster=TRUE
  )
  
  if(window_n_now>=2L){
    
    abline(
      v=window_boundary_now,
      col="grey85",
      lwd=0.25
    )
  }
  
  axis(
    side=2,
    at=c(1,12,24,36,48),
    labels=c(1,12,24,36,48),
    las=1
  )
  
  mtext(
    "Snake row",
    side=2,
    line=3
  )
  
  mtext(
    "Sliding-window order",
    side=1,
    line=2
  )
  
  title(
    main=paste0(
      "Chr",chr_now,
      " - ",
      target_id[2],
      " (Oldest Control, Age 63)"
    )
  )
  
  box()
  
  dev.off()
}
system(paste("open",shQuote(output_dir_v9)))



#儲存目前完整R Environment
save.image(file="/Users/sheenazhengmei/Desktop/summerintern/Version9casecontrol.RData")
#下次開啟
load("/Users/sheenazhengmei/Desktop/summerintern/Version9casecontrol.RData")






