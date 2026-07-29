#parallel v2
#一樣用parallel2的檔案：snp_master_v2_n----
#將genotype 改成nucleotide----
head(snp_master_v2_n)


#確認Version 2需要的欄位都存在
required_v2_cols<-c("probeset_id","Allele_A","Allele_B")
required_v2_cols%in%names(snp_master_v2_n)

#分類genotype 
allele_check<-snp_master_v2_n[,.(image_order,probeset_id,Chr_id,Start,Allele_A,Allele_B)]

allele_check[,Allele_A:=toupper(trimws(Allele_A))]
allele_check[,Allele_B:=toupper(trimws(Allele_B))]

allele_check[,allele_status:=fcase(
  is.na(Allele_A)|Allele_A==""|
    is.na(Allele_B)|Allele_B=="",
  "Missing allele annotation",
  
  Allele_A%in%c("A","C","G","T")&
    Allele_B%in%c("A","C","G","T")&
    Allele_A!=Allele_B,
  "Canonical biallelic SNV",    #標準雙等位單核苷酸變異:變異只涉及一個 nucleotide 位置
  
  Allele_A%in%c("A","C","G","T")&
    Allele_B%in%c("A","C","G","T")&
    Allele_A==Allele_B,
  "Same nucleotide in A/B",
  
  Allele_A=="-"|Allele_B=="-",
  "Gap-containing allele",
  
  default="Multi-base or other allele"
)]

allele_check[,.N,by=allele_status][order(-N)]

#百分比
allele_status<-allele_check[,.N,by=allele_status][order(-N)]
allele_status[,Percentage:=round(N/sum(N)*100,2)]

#分開檢查
#建立Allele_A與Allele_B的組合
allele_check[,allele_pair:=paste0(Allele_A,"/",Allele_B)]

#統計每種allele pair的數量與比例
allele_pair_summary<-allele_check[,.N,by=.(allele_status,allele_pair)][order(allele_status,-N)]
allele_pair_summary[,Percentage_within_status:=round(N/sum(N)*100,2),by=allele_status]

#SNV:631696----
allele_pair_summary[allele_status=="Canonical biallelic SNV"]
#含 gap 的 allele:看它們是不是insertion/deletion marker -/A G/- :48034----
allele_pair_summary[allele_status=="Gap-containing allele"]
#多鹼基或其他 allele :1114----
allele_pair_summary[allele_status=="Multi-base or other allele"]



#先針對SNV:631696----
#取得631696個標準SNV
canonical_snv_order<-allele_check[allele_status=="Canonical biallelic SNV",image_order]
length(canonical_snv_order)

#保留631,696個標準雙等位SNV，copy避免修改原始資料
snp_v2_snv<-copy(snp_master_v2_n[image_order%in%canonical_snv_order])
#依原本的全基因組順序重新排列
setorder(snp_v2_snv,image_order)
#保留原始影像位置，方便未來反查
snp_v2_snv[,original_image_order:=image_order]
#建立Version 2專用的連續影像順序
snp_v2_snv[,image_order_v2:=.I]

#確認筆數與新順序
nrow(snp_v2_snv)
range(snp_v2_snv$image_order_v2)

#確認兩個allele都只能是A、C、G、T
table(snp_v2_snv$Allele_A,useNA="ifany")
table(snp_v2_snv$Allele_B,useNA="ifany")

#確認沒有相同的Allele_A與Allele_B
sum(snp_v2_snv$Allele_A==snp_v2_snv$Allele_B)

#確認所有marker都符合標準雙等位SNV條件
all(snp_v2_snv$Allele_A%in%c("A","C","G","T")&
    snp_v2_snv$Allele_B%in%c("A","C","G","T")&
    snp_v2_snv$Allele_A!=snp_v2_snv$Allele_B)

#找出所有以_call_code結尾的10位樣本genotype欄位
genotype_cols<-grep("_call_code$",names(snp_v2_snv),value=TRUE)

#取得10位樣本所有原始genotype call的不同值
original_call_values<-sort(unique(unlist( snp_v2_snv[,..genotype_cols],use.names=FALSE)))


#建立每列 SNP 的三種 nucleotide genotype
#AA → Allele_A / Allele_A
#AB → Allele_A / Allele_B
#BB → Allele_B / Allele_B
#建立AA對應的實際homozygous nucleotide genotype
snp_v2_snv[,nt_AA:=paste0(Allele_A,"/",Allele_A)]

#建立BB對應的實際homozygous nucleotide genotype
snp_v2_snv[,nt_BB:=paste0(Allele_B,"/",Allele_B)]

#建立AB對應的heterozygous genotype，並統一字母順序
snp_v2_snv[,nt_AB:=paste0(pmin(Allele_A,Allele_B),"/",pmax(Allele_A,Allele_B))]

##查看六種allele pair分別對應的AA、AB與BB
unique(snp_v2_snv[,.(Allele_A,Allele_B,nt_AA,nt_AB,nt_BB)])[order(Allele_A,Allele_B)]
head((snp_v2_snv))

#轉乘nucleotide版本 在snp_v2_snv新增nucleotide欄位----

#逐一轉換10位樣本的genotype call
for(col in genotype_cols){
  #建立新的nucleotide欄位名稱
  new_col<-sub("_call_code$","_nucleotide",col)
  #依AA、AB、BB轉換成該SNP真正的nucleotide genotype
  snp_v2_snv[,(new_col):=fcase(
    get(col)=="AA",nt_AA,
    get(col)=="AB",nt_AB,
    get(col)=="BB",nt_BB,
    get(col)=="NoCall","NoCall",
    get(col)=="Y_not_applicable","Y_not_applicable",
    default=NA_character_)]
}

#建立10個nucleotide genotype欄位名稱
nucleotide_cols<-sub("_call_code$","_nucleotide",genotype_cols)

#取得所有樣本轉換後出現的nucleotide genotype類別
observed_nucleotide_values<-sort(unique(unlist(snp_v2_snv[,..nucleotide_cols],use.names=FALSE)))

#統計10位樣本轉換後每種nucleotide genotype的數量
nucleotide_distribution<-rbindlist(
  lapply(seq_along(nucleotide_cols),function(i){
    col<-nucleotide_cols[i]
    snp_v2_snv[
      ,.N,
      by=.(genotype_value=get(col))
    ][
      ,sample:=sub("_nucleotide$","",col)
    ]
  }),
  use.names=TRUE
)

#調整欄位順序並顯示
setcolorder(nucleotide_distribution,c("sample","genotype_value","N"))

nucleotide_distribution[,Percentage:=round(N/sum(N)*100,3),by=sample]

nucleotide_distribution[order(sample,-N)]

#轉換確認完成後移除暫時輔助欄位
snp_v2_snv[,c("nt_AA","nt_AB","nt_BB"):=NULL]
head(snp_v2_snv)

#image----
#計算Version 2影像邊長與padding數量

#nrow:631696 
#sqrt:多少的平方會大於等於631696->795
image_side_v2<-ceiling(sqrt(nrow(snp_v2_snv))) 
#padding:795*795-631696=329 
padding_n_v2<-image_side_v2^2-nrow(snp_v2_snv)

#上色----
#homozygous 用較明顯的單色 
#heterozygous 用混合色 
#missing 用灰色 
#padding 用黑色 


#建立Version 2A的顏色與代碼對照表；先合併NoCall與Y_not_applicable
#每一種異型合子的 RGB 值，都是其兩種 nucleotide 基礎顏色的平均值
color_encoding_v2A<-data.table(
  genotype_class=c("Padding","A/A","A/C","A/G","A/T","C/C","C/G","C/T","G/G","G/T","T/T","NoCall","Y_not_applicable"),
  code=c(0L,1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,11L),
  R=c(0,255,128,128,255,0,0,128,0,128,255,128,128),
  G=c(0,0,0,128,82,0,128,82,255,210,165,128,128),
  B=c(0,0,255,0,0,255,128,128,0,0,0,128,128)
)
#A/C改成A與C基礎色的平均值
color_encoding_v2A[genotype_class=="A/C",`:=`(R=128,G=0,B=128)]
#查看顏色表
color_encoding_v2A

#取得所有樣本轉換後出現的nucleotide genotype類別
observed_nucleotide_values<-sort(unique(unlist(snp_v2_snv[,..nucleotide_cols],use.names=FALSE)))
#顯示顏色表中定義的類別
sort(unique(color_encoding_v2A$genotype_class))

#建立genotype文字到數值code的命名向量；同名NoCall與Y_not_applicable都能對應
code_lookup_v2A<-setNames(color_encoding_v2A$code,color_encoding_v2A$genotype_class)

#將 10 位樣本的 nucleotide genotype 轉成數值 code
for(col in nucleotide_cols){
  #建立新的code欄位名稱；例如DM-002_..._nucleotide 轉成 DM-002_..._code_v2A
  new_col<-sub("_nucleotide$","_code_v2A",col)
  #利用lookup向量將文字類別轉成數值代碼
  snp_v2_snv[,(new_col):=unname(code_lookup_v2A[get(col)])]
}

code_cols_v2A<-sub("_nucleotide$","_code_v2A",nucleotide_cols)

#做一張每位樣本的 code 分布表

code_distribution_v2A<-rbindlist(
  lapply(seq_along(code_cols_v2A),function(i){
    col<-code_cols_v2A[i]
    snp_v2_snv[
      ,.N,
      by=.(code_value=get(col))
    ][
      ,sample:=sub("_code_v2A$","",col)
    ]
  }),
  use.names=TRUE
)

#調整欄位順序
setcolorder(code_distribution_v2A,c("sample","code_value","N"))

#計算每位樣本內百分比
code_distribution_v2A[,Percentage:=round(N/sum(N)*100,3),by=sample]
code_distribution_v2A[order(sample,code_value)]

#把 code 對回 genotype 名稱
#建立code對應的主要genotype名稱；因Version 2A將兩種missing合併，因此用NoCall代表code 11
code_label_v2A<-data.table(
  code_value=c(0L,1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L),
  genotype_label=c("Padding","A/A","A/C","A/G","A/T","C/C","C/G","C/T","G/G","G/T","T/T","NoCall")
)

#將code名稱合併回分布表
code_distribution_v2A<-merge(
  code_distribution_v2A,
  code_label_v2A,
  by="code_value",
  all.x=TRUE,
  sort=FALSE
)

#重新調整欄位順序
setcolorder(code_distribution_v2A,c("sample","code_value","genotype_label","N","Percentage"))
code_distribution_v2A[order(sample,code_value)]

#建立每位樣本的 795 × 795 數值矩陣----
#從code欄位名稱取出樣本名稱
sample_names_v2A<-sub("_code_v2A$","",code_cols_v2A)

#逐位樣本將631,696個code補329個0，再依列填入795×795矩陣
image_matrix_v2A<-setNames(
  lapply(code_cols_v2A,function(col){
    code_vector<-as.integer(snp_v2_snv[[col]])      #取出一位樣本的631,696個SNV code
    code_vector<-c(code_vector,rep(0L,padding_n_v2)) #在序列最後補329個Padding code 0
    matrix(
      code_vector,
      nrow=image_side_v2,
      ncol=image_side_v2,
      byrow=TRUE#依原始SNV順序由左至右、由上至下填入
    )
  }),
  sample_names_v2A
)


#建立每位樣本的數值矩陣QC結果
matrix_qc_v2A<-rbindlist(
  lapply(seq_along(image_matrix_v2A),function(i){
    m<-image_matrix_v2A[[i]]
    data.table(
      sample=names(image_matrix_v2A)[i],
      matrix_row=nrow(m),
      matrix_col=ncol(m),
      SNV_pixel=sum(m!=0L),
      padding_pixel=sum(m==0L),
      min_code=min(m),
      max_code=max(m),
      check=nrow(m)==795L&
        ncol(m)==795L&
        sum(m!=0L)==631696L&
        sum(m==0L)==329L
    )
  })
)

#確認相同code是否只對應同一組RGB；所有n_RGB都應為1
color_code_qc_v2A<-color_encoding_v2A[
  ,.(n_RGB=uniqueN(paste(R,G,B,sep="_")),
     genotype_label=paste(genotype_class,collapse=", ")),
  by=code
]
#每個code只保留一組RGB，避免code 11重複
rgb_lookup_v2A<-unique(color_encoding_v2A[,.(code,R,G,B)],by="code")

#依code由0至11排列
setorder(rgb_lookup_v2A,code)

#數值矩陣轉 RGB array
code_to_rgb_v2A<-function(code_matrix,lookup=rgb_lookup_v2A){
  #依每個pixel的code找出對應RGB列
  code_index<-match(code_matrix,lookup$code)
  #若有任何code無法配對，立即停止
  if(anyNA(code_index)){
    stop("Code matrix contains undefined color codes.")
  }
  #建立三通道RGB array，第三維依序為R、G、B
  rgb_array<-array(
    0,
    dim=c(nrow(code_matrix),ncol(code_matrix),3),
    dimnames=list(NULL,NULL,c("R","G","B"))
  )
  #writePNG需要0至1，因此將0至255除以255
  rgb_array[,,1]<-matrix(lookup$R[code_index]/255,nrow=nrow(code_matrix),ncol=ncol(code_matrix))
  rgb_array[,,2]<-matrix(lookup$G[code_index]/255,nrow=nrow(code_matrix),ncol=ncol(code_matrix))
  rgb_array[,,3]<-matrix(lookup$B[code_index]/255,nrow=nrow(code_matrix),ncol=ncol(code_matrix))
  return(rgb_array)
}


#將10張數值矩陣全部轉成795×795×3 RGB array
rgb_image_v2A<-lapply(
  image_matrix_v2A,
  code_to_rgb_v2A
)
#保留樣本名稱
names(rgb_image_v2A)<-names(image_matrix_v2A)

#選擇第一位樣本進行視覺檢查
sample_view_v2A<-names(rgb_image_v2A)[1]

#利用rasterImage顯示，interpolate=FALSE避免像素被平滑
old_par<-par(mar=c(0,0,2,0))
plot.new()
plot.window(xlim=c(0,1),ylim=c(0,1),asp=1)
rasterImage(
  as.raster(rgb_image_v2A[[sample_view_v2A]]),
  0,0,1,1,
  interpolate=FALSE
)
title(main=paste0(sample_view_v2A,"Version 2A nucleotide encoding"))
par(old_par)

#10個人的
#建立Version 2A輸出資料夾
dir.create("Version2A_nucleotide_images",showWarnings=FALSE)

#逐位樣本輸出PNG；檔名中的特殊符號改成底線
for(i in seq_along(rgb_image_v2A)){
  sample_name<-names(rgb_image_v2A)[i]
  safe_name<-gsub("[^A-Za-z0-9_-]","_",sample_name)
  writePNG(
    rgb_image_v2A[[i]],
    target=file.path(
      "Version2A_nucleotide_images",
      paste0(safe_name,"_Version2A.png")
    )
  )
}
#10張

sample_cols_v2A<-names(rgb_image_v2A)
stopifnot(length(sample_cols_v2A)==10)

#保存目前繪圖設定，避免後續圖形受到影響
old_par<-par(no.readonly=TRUE)

#將10位樣本排成2列×5欄，oma預留整張圖的總標題空間
par(mfrow=c(2,5),mar=c(0.5,0.5,2,0.5),oma=c(0,0,2,0))

for(x in sample_cols_v2A){
  #建立正方形繪圖區
  plot.new()
  plot.window(xlim=c(0,1),ylim=c(0,1),asp=1)
  #顯示目前樣本的Version 2A nucleotide影像
  rasterImage(
    as.raster(rgb_image_v2A[[x]]),
    0,0,1,1,
    interpolate=FALSE
  )
  #移除樣本名稱中Axiom_TPM之後的文字
  short_name<-sub("_\\(Axiom_TPM\\).*","",x)
  #加入樣本名稱
  title(main=short_name,cex.main=0.8,line=0.3)
}

#加入整張圖的總標題
mtext(
  "Version 2A: Nucleotide-based encoding",
  outer=TRUE,
  cex=1.2,
  font=2
)

#恢復原本繪圖設定
par(old_par)



#列出已輸出的10張Version 2A影像
list.files("Version2A_nucleotide_images",pattern="\\.png$")


#依Version 2的新image_order計算每個SNV的pixel列與欄
snp_v2_snv[,pixel_row_v2:=((image_order_v2-1L)%/%image_side_v2)+1L]
snp_v2_snv[,pixel_col_v2:=((image_order_v2-1L)%%image_side_v2)+1L]

#建立Version 2的pixel與SNV對照表
pixel_snp_map_v2<-snp_v2_snv[
  ,.(image_order_v2,
     original_image_order,
     pixel_row_v2,
     pixel_col_v2,
     probeset_id,
     dbSNP_RS_ID,
     Chr_id,
     Start,
     Allele_A,
     Allele_B)
]

#version2b----
#建立Version 2B顏色表；NoCall與Y_not_applicable使用不同code與顏色
color_encoding_v2B<-data.table(
  genotype_class=c("Padding","A/A","A/C","A/G","A/T","C/C","C/G","C/T","G/G","G/T","T/T","NoCall","Y_not_applicable"),
  code=c(0L,1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,12L),
  R=c(0,255,128,128,255,0,0,128,0,128,255,128,255),
  G=c(0,0,0,128,82,0,128,82,255,210,165,128,255),
  B=c(0,0,128,0,0,255,128,128,0,0,0,128,0)
)

#查看Version 2B顏色表
color_encoding_v2B

#取得10位樣本實際出現的nucleotide genotype類別
observed_nucleotide_values<-sort(
  unique(
    unlist(
      snp_v2_snv[,..nucleotide_cols],
      use.names=FALSE
    )
  )
)
#建立genotype文字對應Version 2B數值code的命名向量
code_lookup_v2B<-setNames(
  color_encoding_v2B$code,
  color_encoding_v2B$genotype_class
)

#查看lookup內容
code_lookup_v2B

#逐一將10位樣本的nucleotide genotype轉成Version 2B code
for(col in nucleotide_cols){
  #建立新欄位名稱
  new_col<-sub("_nucleotide$","_code_v2B",col)
  #依genotype文字查找對應code並存入新欄位
  snp_v2_snv[,(new_col):=as.integer(
    unname(code_lookup_v2B[get(col)])
  )]
}

#建立10個Version 2B code欄位名稱
code_cols_v2B<-sub(
  "_nucleotide$",
  "_code_v2B",
  nucleotide_cols
)

#確認欄位數量與是否全部成功建立
length(code_cols_v2B)
all(code_cols_v2B%in%names(snp_v2_snv))

#合併10位樣本的Version 2B code，用於QC
all_codes_v2B<-unlist(
  snp_v2_snv[,..code_cols_v2B],
  use.names=FALSE
)

#從code欄位名稱取得樣本名稱
sample_names_v2B<-sub(
  "_code_v2B$",
  "",
  code_cols_v2B
)

#每位樣本補329個Padding code 0，再轉成795×795矩陣
image_matrix_v2B<-setNames(
  lapply(code_cols_v2B,function(col){
    code_vector<-as.integer(snp_v2_snv[[col]])       #取出631,696個SNV code
    code_vector<-c(code_vector,rep(0L,padding_n_v2)) #序列最後加入329個Padding
    matrix(
      code_vector,
      nrow=image_side_v2,
      ncol=image_side_v2,
      byrow=TRUE
    )
  }),
  sample_names_v2B
)

#確認影像矩陣數量與尺寸
length(image_matrix_v2B)
dim(image_matrix_v2B[[1]])


#確認每位樣本的矩陣尺寸、SNV與Padding數量
matrix_qc_v2B<-rbindlist(
  lapply(seq_along(image_matrix_v2B),function(i){
    m<-image_matrix_v2B[[i]]
    data.table(
      sample=names(image_matrix_v2B)[i],
      matrix_row=nrow(m),
      matrix_col=ncol(m),
      SNV_pixel=sum(m!=0L),
      padding_pixel=sum(m==0L),
      y_not_applicable_pixel=sum(m==12L),
      check=nrow(m)==image_side_v2&
        ncol(m)==image_side_v2&
        sum(m!=0L)==nrow(snp_v2_snv)&
        sum(m==0L)==padding_n_v2
    )
  })
)

matrix_qc_v2B
all(matrix_qc_v2B$check)


#建立Version 2B唯一的code與RGB對照表
rgb_lookup_v2B<-color_encoding_v2B[
  ,.(code,R,G,B)
]

#依code排列
setorder(rgb_lookup_v2B,code)

#確認code完整包含0至12
rgb_lookup_v2B
identical(rgb_lookup_v2B$code,0:12)

#將795×795 code矩陣轉成795×795×3 RGB array
code_to_rgb_v2B<-function(code_matrix,lookup=rgb_lookup_v2B){
  #將每個pixel code配對到lookup中的列
  code_index<-match(code_matrix,lookup$code)
  #若有未定義code則停止
  if(anyNA(code_index)){
    stop("Code matrix contains undefined Version 2B color codes.")
  }
  #建立RGB三通道array
  rgb_array<-array(
    0,
    dim=c(nrow(code_matrix),ncol(code_matrix),3),
    dimnames=list(NULL,NULL,c("R","G","B"))
  )
  #將0至255的RGB轉成0至1
  rgb_array[,,1]<-matrix(
    lookup$R[code_index]/255,
    nrow=nrow(code_matrix),
    ncol=ncol(code_matrix)
  )
  rgb_array[,,2]<-matrix(
    lookup$G[code_index]/255,
    nrow=nrow(code_matrix),
    ncol=ncol(code_matrix)
  )
  rgb_array[,,3]<-matrix(
    lookup$B[code_index]/255,
    nrow=nrow(code_matrix),
    ncol=ncol(code_matrix)
  )
  return(rgb_array)
}
#將10位樣本的Version 2B code矩陣轉成RGB array
rgb_image_v2B<-lapply(
  image_matrix_v2B,
  code_to_rgb_v2B
)

#保留樣本名稱
names(rgb_image_v2B)<-names(image_matrix_v2B)

#確認數量與影像維度
length(rgb_image_v2B)
dim(rgb_image_v2B[[1]])

#取得Version 2B的10位樣本名稱
sample_cols_v2B<-names(rgb_image_v2B)

#保存目前繪圖參數
old_par<-par(no.readonly=TRUE)

#將10張圖排成2列×5欄
par(mfrow=c(2,5),mar=c(0.5,0.5,2,0.5),oma=c(0,0,2,0))

for(x in sample_cols_v2B){
  #建立正方形繪圖區
  plot.new()
  plot.window(xlim=c(0,1),ylim=c(0,1),asp=1)
  #顯示目前樣本的Version 2B影像
  rasterImage(
    as.raster(rgb_image_v2B[[x]]),
    0,0,1,1,
    interpolate=FALSE
  )
  #縮短樣本名稱
  short_name<-sub("_\\(Axiom_TPM\\).*","",x)
  #加入樣本名稱
  title(main=short_name,cex.main=0.8,line=0.3)
}

#加入總標題
mtext(
  "Version 2B: Nucleotide encoding with separated missing states",
  outer=TRUE,
  cex=1.1,
  font=2
)

#恢復原始繪圖參數
par(old_par)
#其中一位
#選擇第一位樣本進行視覺檢查
sample_view_v2B<-names(rgb_image_v2B)[1]

#利用rasterImage顯示，interpolate=FALSE避免像素被平滑
old_par<-par(mar=c(0,0,2,0))
plot.new()
plot.window(xlim=c(0,1),ylim=c(0,1),asp=1)
rasterImage(
  as.raster(rgb_image_v2B[[sample_view_v2B]]),
  0,0,1,1,
  interpolate=FALSE
)
title(main=paste0(sample_view_v2B,"Version 2B nucleotide encoding"))
par(old_par)


#


#version 3 加入cytoband (snp_v2_snv->snp_v3_cytoband)----
#先確genotype_data中的cytoband缺失值
#有787 snp(0.12%) 沒有
cytoband<-genotype_data[,.(probeset_id,Chr_id,Start,Cytoband=as.character(Cytoband))]
cytoband[,cytoband_clean:=trimws(Cytoband)] #移除Cytoband前後空白
cytoband[is.na(cytoband_clean)|cytoband_clean%in%c("",".","NA","N/A","Null","null"),cytoband_clean:=NA_character_]#空白 缺失轉的NA
cytoband_s<-cytoband[,.(total_SNP=.N,missing_cytoband=sum(is.na(cytoband_clean)),percentage=round(mean(is.na(cytoband_clean))*100,2))]

#取出所有缺少Cytoband的SNP
missing_cytoband_<-cytoband[is.na(cytoband_clean)]
#依染色體統計Cytoband缺失數量與比例
cytoband_snp<-cytoband[,.(total_SNP=.N,missing_cytoband=sum(is.na(cytoband_clean)),
                        percentage=round(mean(is.na(cytoband_clean))*100,2)
                        ),
                        by=Chr_id
                     ][order(Chr_id)]
#看每條染色體缺失情況
chr_levels<-c(as.character(1:22),"X","Y","MT")

#主要染色體依1–22、X、Y、MT排序；其他特殊Chr放在最後
cytoband[,chr_order:=match(Chr_id,chr_levels)]
cytoband[is.na(chr_order),chr_order:=length(chr_levels)+1L]

#依染色體統計Cytoband缺失數與比例
cytoband_snp<-cytoband[
  ,.(total_SNP=.N,
     missing_cytoband=sum(is.na(cytoband_clean)),
     available_cytoband=sum(!is.na(cytoband_clean)),
     percentage=round(mean(is.na(cytoband_clean))*100,2)),
  by=.(Chr_id,chr_order)
][order(chr_order,Chr_id)]

#移除只用於排序的chr_order欄位
cytoband_snp[,chr_order:=NULL]

#利用snp_v2_snv去看缺失
#把cytobamd加入snp_v2_snv
#確認同一個probeset_id是否對應超過1個非缺失Cytoband
cytoband_map_qc<-cytoband[,.(n_cytoband=uniqueN(na.omit(cytoband_clean))),by=probeset_id][n_cytoband>1]
nrow(cytoband_map_qc)

#每個probeset_id保留唯一的非缺失Cytoband；若完全沒有則保留NA
cytoband_map<-cytoband[
  ,.(Cytoband={
    x<-unique(na.omit(cytoband_clean))
    if(length(x)==0) NA_character_ else x[1]
  }),
  by=probeset_id
]
#確認probeset_id在對照表中沒有重複
anyDuplicated(cytoband_map$probeset_id)

#記錄加入前的SNV筆數
n_before_cytoband<-nrow(snp_v2_snv)

#依probeset_id將清理後的Cytoband加入snp_v2_snv
snp_v2_snv[cytoband_map,on="probeset_id",Cytoband:=i.Cytoband]

#確認加入後資料列數沒有改變
nrow(snp_v2_snv)
nrow(snp_v2_snv)==n_before_cytoband
#確認是否有SNV的probeset_id完全不存在於cytoband_map
sum(!snp_v2_snv$probeset_id%in%cytoband_map$probeset_id)

#統計631,696個SNV中Cytoband的完整情況
cytoband_v2_summary<-snp_v2_snv[,.(total_SNV=.N,
     missing_cytoband=sum(is.na(Cytoband)),
     available_cytoband=sum(!is.na(Cytoband)),
     percentage=round(mean(is.na(Cytoband))*100,2))
]

#看每條染色體缺失情況
head(snp_v2_snv)
chr_level<-c(as.character(1:22),"X","Y","MT")
snp_v2_snv[,chr_order_v2:=match(Chr_id,chr_level)]

#依染色體統計Version 2的Cytoband缺失數量與比例
cytoband_v2_by_chr<-snp_v2_snv[
  ,.(total_SNV=.N,
     missing_cytoband=sum(is.na(Cytoband)),
     available_cytoband=sum(!is.na(Cytoband)),
     percentage=round(mean(is.na(Cytoband))*100,4)),
  by=.(Chr_id,chr_order_v2)
][order(chr_order_v2)]

#移除只用於排序的欄位
cytoband_v2_by_chr[,chr_order_v2:=NULL]

#snp_v3_cytoband n=630941----
#cytoband共755個缺失其中 754都是粒線體 但歷險體本來就不適用這個變編號全數排除

#盤除第17條染色體中1個cytoband缺失值
removed_cytoband_snv<-snp_v2_snv[Chr_id!="MT"&is.na(Cytoband),.(probeset_id,dbSNP_RS_ID,Chr_id,Start,Allele_A,Allele_B)]
#排除的SNV數量為1
nrow(removed_cytoband_snv)

#建立Version 3資料，cytoband共755個缺失其中 754都是粒線體 但歷險體本來就不適用這個變編號 排除Chr17中缺少Cytoband的1個SNV
snp_v3_cytoband<-copy(snp_v2_snv[!is.na(Cytoband)])

#保留Version 2原本的影像順序，方便後續追蹤
snp_v3_cytoband[,original_image_order_v2:=image_order_v2]
sum(is.na(snp_v3_cytoband$Cytoband))

##查看前30種Cytoband格式
head(sort(unique(snp_v3_cytoband$Cytoband)),100)

#查看每條染色體的Cytoband種類數
cytoband_count_by_chr<-snp_v3_cytoband[,.(SNV_N=.N,cytoband_N=uniqueN(Cytoband)),by=Chr_id
][order(match(Chr_id,c(as.character(1:22),"X","Y")))]

snp_v3_cytoband[,Start:=as.numeric(Start)]#數值型
sum(is.na(snp_v3_cytoband$Start))

#染色體自然順序
#vhr:1–22、X、Y
chr_levels_v3<-c(as.character(1:22),"X","Y")

#Chr_id轉成整數=>Chr1=1、X=23、Y=24
snp_v3_cytoband[,chr_order_v3:=match(Chr_id,chr_levels_v3)]

#取Cytoband第一個字母->建立p短臂或q長臂欄位
snp_v3_cytoband[,cytoband_arm:=substr(Cytoband,1,1)]
table(snp_v3_cytoband$cytoband_arm,useNA="ifany")

#先合併完整的cytoband名字：染色體+Cytoband，例如Chr17+p13.1=17p13.1
snp_v3_cytoband[,cytoband_id:=paste0(Chr_id,Cytoband)]
head(sort(unique(snp_v3_cytoband$cytoband_id)),10)

#排序
#依染色體、實際基因組位置及probeset_id排序
#p端 → p arm → 著絲粒 → q arm → q端
setorder(snp_v3_cytoband,chr_order_v3,Start,probeset_id)

#取得排序後每條染色體實際出現的arm順序
arm_order_qc<-unique(snp_v3_cytoband[,.(Chr_id,chr_order_v3,cytoband_arm)])[order(chr_order_v3)]

#將每條染色體的arm順序合併成p→q文字
arm_order_qc<-arm_order_qc[,.(arm_sequence=paste(cytoband_arm,collapse="→")),by=.(Chr_id,chr_order_v3)][order(chr_order_v3)]

#確認順序是否為p→q；只有單一arm時也先視為可接受
arm_order_qc[,check:=arm_sequence%in%c("p→q","p","q")]

arm_order_qc
all(arm_order_qc$check)

#統計每個Cytoband中最小Start、最大Start及SNV數量
cytoband_order_table<-snp_v3_cytoband[,.(band_start=min(Start),band_end=max(Start),SNV_N=.N),by=.(Chr_id,chr_order_v3,cytoband_id,Cytoband,cytoband_arm)]

#依染色體和band_start排列Cytoband
setorder(cytoband_order_table,chr_order_v3,band_start,band_end)

#在每條染色體內建立Cytoband連續順序
cytoband_order_table[,cytoband_order:=seq_len(.N),by=Chr_id]


#查看第1號染色體所有Cytoband的實際排列
cytoband_order_table[Chr_id=="1",.(cytoband_order,cytoband_id,cytoband_arm,band_start,band_end,SNV_N)]


#依cytoband_id將Cytoband順序加入每一個SNV
snp_v3_cytoband[cytoband_order_table,on=.(Chr_id,cytoband_id),cytoband_order:=i.cytoband_order]

#確認所有SNV都成功取得Cytoband順序
sum(is.na(snp_v3_cytoband$cytoband_order))
#依染色體、Cytoband實際順序、Start及probeset_id排列
setorder(snp_v3_cytoband,chr_order_v3,cytoband_order,Start,probeset_id)

#保留Version 2原始順序，並建立Version 3新順序
snp_v3_cytoband[,image_order_v3:=.I]
head(snp_v3_cytoband)

#加入Cytoband_separator snp_v3_layout----
#找出10個樣本的基因型欄位
sample_cols_v3<-grep("_nucleotide$",names(snp_v3_layout),value=TRUE)
#在separator資料列的10個nucleotide欄位填入Cytoband_separator
snp_v3_layout[row_type=="separator",(sample_cols_v3):="Cytoband_separator"]
#查看Version 3核苷酸欄位中的所有可能值
sort(unique(unlist(snp_v3_layout[,..sample_cols_v3],use.names=FALSE)))
#確認每個樣本的Cytoband_separator數量都等於Cytoband區塊數
separator_qc_v3<-data.table(sample=sample_cols_v3,separator_N=sapply(sample_cols_v3,function(x) sum(snp_v3_layout[[x]]=="Cytoband_separator",na.rm=TRUE)),expected_N=n_cytoband_v3)

#判斷每個樣本是否都正確
separator_qc_v3[,check:=separator_N==expected_N]
separator_qc_v3

#複製Version 2A編碼表，避免修改原本物件
color_encoding_v3A<-copy(color_encoding_v2A)
#ytoband_separator，新增白色separator
if(!"Cytoband_separator"%in%color_encoding_v3A$genotype_class) color_encoding_v3A<-rbind(color_encoding_v3A,data.table(genotype_class="Cytoband_separator",code=max(color_encoding_v3A$code)+1L,R=255,G=255,B=255))
#完整編碼表
color_encoding_v3A

#建立核苷酸基因型對應數字code的查找表
code_lookup_v3A<-setNames(color_encoding_v3A$code,color_encoding_v3A$genotype_class)

#將核苷酸基因型轉成數字編碼矩陣
#將10個nucleotide欄位依code_lookup_v3A轉成數字code矩陣
image_matrix_v3A<-sapply(sample_cols_v3,function(x) unname(code_lookup_v3A[as.character(snp_v3_layout[[x]])]))

#保留原本的樣本欄位名稱並將矩陣轉成整數
colnames(image_matrix_v3A)<-sample_cols_v3
storage.mode(image_matrix_v3A)<-"integer"

#確認矩陣列數等於layout列數、欄數等於10個樣本
dim(image_matrix_v3A)
nrow(image_matrix_v3A)==nrow(snp_v3_layout)
ncol(image_matrix_v3A)==length(sample_cols_v3)

#確認所有基因型都成功轉成code，正常應為0
sum(is.na(image_matrix_v3A))

#取得Cytoband separator的數字code
separator_code_v3A<-unname(code_lookup_v3A["Cytoband_separator"])

#查看separator使用的code
separator_code_v3A
#統計每個樣本影像中的separator數量
separator_qc_v3A<-data.table(sample=sample_cols_v3,separator_N=colSums(image_matrix_v3A==separator_code_v3A))

#加入預期數量；每個Cytoband後面應有1個separator
separator_qc_v3A[,expected_N:=uniqueN(snp_v3_cytoband$cytoband_id)]

#確認實際與預期數量是否相同
separator_qc_v3A[,check:=separator_N==expected_N]
separator_qc_v3A

#查看第一位樣本各基因型code及separator的數量
data.table(code=image_matrix_v3A[,1])[,.N,by=code][order(code)]

#image
layout_n_v3A<-nrow(image_matrix_v3A)

#計算可容納所有SNP與separator的最小正方形邊長
image_side_v3A<-ceiling(sqrt(layout_n_v3A))

#計算填滿正方形影像所需的黑色padding數量
padding_v3A<-image_side_v3A^2-layout_n_v3A

#整理影像尺寸資訊
image_size_v3A<-data.table(layout_pixels=layout_n_v3A,image_side=image_side_v3A,total_pixels=image_side_v3A^2,padding_pixels=padding_v3A)

image_size_v3A
#依layout順序填入數字code，最後不足的部分補上code 0作為黑色padding
code_image_v3A<-lapply(seq_along(sample_cols_v3),function(i) matrix(c(image_matrix_v3A[,i],rep(0L,padding_v3A)),nrow=image_side_v3A,ncol=image_side_v3A,byrow=TRUE))

#將10張影像依樣本名稱命名
names(code_image_v3A)<-sample_cols_v3

#確認每張影像尺寸皆相同
sapply(code_image_v3A,dim)

#確認第一張影像中的黑色padding數量正確
sum(code_image_v3A[[1]]==0L)==padding_v3A

#確認第一張影像中的separator數量等於Cytoband區塊數量
sum(code_image_v3A[[1]]==separator_code_v3A)==uniqueN(snp_v3_cytoband$cytoband_id)

#確認同一個code沒有對應到不同RGB顏色，正常應全部為1
color_encoding_v3A[,uniqueN(paste(R,G,B,sep="_")),by=code]

#每個code只保留一組RGB，建立code對應顏色表
rgb_table_v3A<-unique(color_encoding_v3A[,.(code,R,G,B)],by="code")

#建立索引為code+1的RGB查找矩陣，因為R矩陣索引從1開始而Padding code為0
rgb_lookup_v3A<-matrix(0,nrow=max(rgb_table_v3A$code)+1L,ncol=3)
rgb_lookup_v3A[rgb_table_v3A$code+1L,]<-as.matrix(rgb_table_v3A[,.(R,G,B)])/255

#將單張code影像轉為高度×寬度×3的RGB陣列
code_to_rgb_v3A<-function(code_matrix){rgb_array<-array(0,dim=c(nrow(code_matrix),ncol(code_matrix),3));for(k in 1:3) rgb_array[,,k]<-matrix(rgb_lookup_v3A[as.integer(code_matrix)+1L,k],nrow=nrow(code_matrix),ncol=ncol(code_matrix));rgb_array}

#轉換10位樣本的Version 3A影像
rgb_image_v3A<-lapply(code_image_v3A,code_to_rgb_v3A)

#保留樣本名稱
names(rgb_image_v3A)<-names(code_image_v3A)

#找出第一個separator在layout中的位置
first_separator_v3A<-which(snp_v3_layout$row_type=="separator")[1]

#轉換成影像中的row與column
separator_row_v3A<-ceiling(first_separator_v3A/image_side_v3A)
separator_col_v3A<-(first_separator_v3A-1L)%%image_side_v3A+1L

#查看該pixel的RGB，預期為1、1、1
rgb_image_v3A[[1]][separator_row_v3A,separator_col_v3A,]
#查看影像最後一個pixel，若有padding則預期為0、0、0
rgb_image_v3A[[1]][image_side_v3A,image_side_v3A,]

#選擇要比較的樣本編號，1代表第一位樣本
sample_i<-1

#確認Version 2A與Version 3A選到的樣本名稱
names(rgb_image_v2A)[sample_i]
names(rgb_image_v3A)[sample_i]

#將Version 2A與Version 3A並排顯示
par(mfrow=c(1,2),mar=c(1,1,3,1))
plot.new()
plot.window(xlim=c(0,1),ylim=c(0,1),asp=1)
rasterImage(rgb_image_v2A[[sample_i]],0,0,1,1,interpolate=FALSE)
title(main="Version 2A\nNucleotide encoding")
plot.new()
plot.window(xlim=c(0,1),ylim=c(0,1),asp=1)
rasterImage(rgb_image_v3A[[sample_i]],0,0,1,1,interpolate=FALSE)
title(main="Version 3A\nNucleotide + Cytoband separator")
par(mfrow=c(1,1))

#取得第一個separator在一維layout中的位置
first_sep<-which(snp_v3_layout$row_type=="separator")[1]

#將一維位置換算為影像的row與column
first_sep_row<-ceiling(first_sep/image_side_v3A)
first_sep_col<-(first_sep-1L)%%image_side_v3A+1L

#設定separator周圍的放大範圍
row_range<-max(1,first_sep_row-10):min(image_side_v3A,first_sep_row+10)
col_range<-max(1,first_sep_col-20):min(image_side_v3A,first_sep_col+20)

#顯示separator附近區域
par(mar=c(1,1,3,1))
plot.new()
plot.window(xlim=c(0,1),ylim=c(0,1),asp=1)
rasterImage(rgb_image_v3A[[sample_i]][row_range,col_range,,drop=FALSE],0,0,1,1,interpolate=FALSE)
title(main=paste0("First Cytoband boundary: ",snp_v3_layout$cytoband_id[first_sep]))

#查看第一個separator前後的layout內容
snp_v3_layout[(first_sep-2):(first_sep+2),.(layout_order_v3,row_type,Chr_id,cytoband_id,Start)]

#將完整欄位名稱簡化成樣本名稱，例如DM-002與EOAD_P10
sample_titles_v3<-sub("_\\(Axiom_TPM\\).*","",names(rgb_image_v3A))

#確認共有10張影像及10個樣本名稱
length(rgb_image_v3A)
sample_titles_v3


#Mac開啟較大的繪圖視窗
quartz(width=15,height=7)

#將10張Version 3A影像畫在新視窗中
par(mfrow=c(2,5),mar=c(1,1,2.5,1),oma=c(0,0,2,0))
for(i in seq_along(rgb_image_v3A)){
  plot.new()
  plot.window(xlim=c(0,1),ylim=c(0,1),asp=1)
  rasterImage(rgb_image_v3A[[i]],0,0,1,1,interpolate=FALSE)
  title(main=sample_titles_v3[i],cex.main=0.9)
}
mtext("Version 3A: Nucleotide Encoding with Cytoband Separators",side=3,outer=TRUE,cex=1.2)
par(mfrow=c(1,1))

#輸出10位樣本的Version 3A影像為高解析度PNG
png("Version3A_10_samples.png",width=5000,height=2400,res=300)
par(mfrow=c(2,5),mar=c(1,1,2.5,1),oma=c(0,0,2,0))
for(i in seq_along(rgb_image_v3A)){
  plot.new()
  plot.window(xlim=c(0,1),ylim=c(0,1),asp=1)
  rasterImage(rgb_image_v3A[[i]],0,0,1,1,interpolate=FALSE)
  title(main=sample_titles_v3[i],cex.main=0.9)
}
mtext("Version 3A: Nucleotide Encoding with Cytoband Separators",side=3,outer=TRUE,cex=1.2)
dev.off()




#version4 利用LD (使用snp_v3_cytoband)----
#LD+排除Y

#排除Y n=617808----
#建立Version 4 SNP基礎資料，排除Y染色體但暫時保留X染色體
snp_v4_ld<-copy(snp_v3_cytoband[Chr_id!="Y"])

#保留Version 3原始順序，方便後續追蹤
snp_v4_ld[,original_image_order_v3:=image_order_v3]

#重新建立排除Y染色體後的連續SNP順序
snp_v4_ld[,snp_order_v4:=.I]
table(snp_v4_ld$Chr_id,useNA = "ifany")

#確認Y染色體筆數為0
sum(snp_v4_ld$Chr_id=="Y")

#確認新的SNP順序連續且沒有重複
range(snp_v4_ld$snp_order_v4)
anyDuplicated(snp_v4_ld$snp_order_v4)

#確認Version 4只包含1–22與X
setdiff(unique(snp_v4_ld$Chr_id),c(as.character(1:22),"X"))





