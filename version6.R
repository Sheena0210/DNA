#version6

#snp_v2_snv:631696
library(data.table)
snp_v2_snv<-readRDS("/Users/sheenazhengmei/Desktop/summerintern/snp_v2_snv.rds")

#snp_v6_del_y先排除chromosome y(13133):618563
snp_v6_del_y<-copy(snp_v2_snv[Chr_id!="Y"])


#Version 6A----
#step1:依 Start 建立 1 Mb windows----
#設定每個window大小為1,000,000 bp
window_size_bp_v6A<-1000000L

snp_v6_del_y[,Start_num:=as.numeric(Start)]
snp_v6_del_y[,chr_order_v6A:=match(Chr_id,c(as.character(1:22),"X","MT"))]

#計算每個snp所屬的window
snp_v6_del_y[,window_id_v6A:=
               as.integer(floor((Start_num-1)/window_size_bp_v6A)+1L)]

#起點終點
snp_v6_del_y[,window_bp_start_v6A:=(window_id_v6A-1L)*window_size_bp_v6A+1L]
snp_v6_del_y[,window_bp_end_v6A:=window_id_v6A*window_size_bp_v6A]

#window name
snp_v6_del_y[,window_key_v6A:=paste0("Chr",Chr_id,"_W",window_id_v6A)]
#排序
setorder(snp_v6_del_y,chr_order_v6A,Start_num,probeset_id)


#step2:統計每個window的snp----
window_snp_count_v6A<-snp_v6_del_y[,.(SNP_n=.N,observed_start_bp=min(Start_num),
                                      observed_end_bp=max(Start_num),
                                      window_bp_start_v6A=first(window_bp_start_v6A),
                                      window_bp_end_v6A=first(window_bp_end_v6A)),
                                   by=.(chr_order_v6A,Chr_id,window_id_v6A,window_key_v6A)]
#依染色體與window順序排列
setorder(window_snp_count_v6A,chr_order_v6A,window_id_v6A)

#查看前20個window
window_snp_count_v6A[1:20]

#全基因組中單一1 Mb window的最大SNP數:3551->hilbert curve:64*64
max_snp_per_window_v6A<-max(window_snp_count_v6A$SNP_n)
#列出所有並列最大SNP數的window
max_snp_window_v6A<-window_snp_count_v6A[SNP_n==max_snp_per_window_v6A]

max_snp_window_v6A

#每條染色體內SNP數最多的1 Mb window
max_window_each_chr_v6A<-window_snp_count_v6A[
  SNP_n==max(SNP_n),
  by=.(chr_order_v6A,Chr_id)
]

setorder(max_window_each_chr_v6A,chr_order_v6A)
max_window_each_chr_v6A


sum(window_snp_count_v6A$SNP_n)
nrow(snp_v6_del_y)

sum(window_snp_count_v6A$SNP_n)==nrow(snp_v6_del_y)

summary(window_snp_count_v6A$SNP_n)

max_snp_per_window_v6A


#step3 64 × 64 Hilbert template----

#設定Hilbert方格邊長與容量
hilbert_side_v6A<-64L
hilbert_capacity_v6A<-hilbert_side_v6A^2

#Version 6A Step 3：從頭建立Hilbert curve函數----------------------------

#旋轉或翻轉Hilbert子區塊
hilbert_rotate_v6A<-function(side,x,y,rx,ry){
  if(ry==0L){
    if(rx==1L){
      x<-side-1L-x
      y<-side-1L-y
    }
    temp<-x
    x<-y
    y<-temp
  }
  list(x=x,y=y)
}

#將一維Hilbert index轉成二維x、y座標
hilbert_d2xy_v6A<-function(side,index){
  
  #檢查side必須是2的次方，例如8、16、32、64
  if(side<=0L||bitwAnd(side,side-1L)!=0L){
    stop("side必須是2的次方")
  }
  
  #檢查index是否位於方格容量範圍內
  if(any(index<0L)||any(index>=side^2)){
    stop("Hilbert index超出方格容量")
  }
  
  index<-as.integer(index)
  x_result<-integer(length(index))
  y_result<-integer(length(index))
  
  for(i in seq_along(index)){
    x<-0L
    y<-0L
    temp_index<-index[i]
    current_side<-1L
    
    while(current_side<side){
      
      #決定目前Hilbert子區塊的方向
      rx<-bitwAnd(bitwShiftR(temp_index,1L),1L)
      ry<-bitwAnd(bitwXor(temp_index,rx),1L)
      
      #依方向旋轉或翻轉座標
      rotated<-hilbert_rotate_v6A(
        side=current_side,
        x=x,
        y=y,
        rx=rx,
        ry=ry
      )
      
      x<-rotated$x+current_side*rx
      y<-rotated$y+current_side*ry
      
      #每次處理Hilbert index的兩個bit
      temp_index<-bitwShiftR(temp_index,2L)
      current_side<-current_side*2L
    }
    
    x_result[i]<-x
    y_result[i]<-y
  }
  
  data.table(
    hilbert_index_v6A=index,
    x=x_result,
    y=y_result
  )
}

#設定所有1 Mb windows共用的Hilbert方格
hilbert_side_v6A<-64L
hilbert_capacity_v6A<-hilbert_side_v6A^2

#確認容量足以容納最大window的3551個SNP
hilbert_capacity_v6A
hilbert_capacity_v6A>=max_snp_per_window_v6A

#建立index 0–4095的共同Hilbert座標
hilbert_template_v6A<-hilbert_d2xy_v6A(
  side=hilbert_side_v6A,
  index=0:(hilbert_capacity_v6A-1L)
)

#將0-based座標轉為R matrix使用的1-based row與column
hilbert_template_v6A[,`:=`(
  patch_row_v6A=y+1L,
  patch_col_v6A=x+1L
)]

#只保留後續需要的欄位
hilbert_template_v6A<-hilbert_template_v6A[
  ,
  .(
    hilbert_index_v6A,
    patch_row_v6A,
    patch_col_v6A
  )
]

#step4每個window內排序並配置Hilbert座標----
#再次固定排序：染色體→ window→ Start→ probeset_id
setorder(snp_v6_del_y,chr_order_v6A,window_id_v6A,
  Start_num,probeset_id)
#每個window內，依Start由小到大建立SNP順序
snp_v6_del_y[,snp_order_within_window_v6A:=seq_len(.N),
             by=.(Chr_id,window_id_v6A)]

#Hilbert index從0開始
snp_v6_del_y[,hilbert_index_v6A:=
               snp_order_within_window_v6A-1L]
#依共同Hilbert template取得64×64 patch內的row與column
snp_v6_del_y[,patch_row_v6A:=
               hilbert_template_v6A$patch_row_v6A[
                 match(hilbert_index_v6A,
                       hilbert_template_v6A$hilbert_index_v6A)
               ]]

snp_v6_del_y[,patch_col_v6A:=
               hilbert_template_v6A$patch_col_v6A[
                 match(hilbert_index_v6A,
                       hilbert_template_v6A$hilbert_index_v6A)
               ]]



#start5建立每格snp的chromosome----

#資料中最大的window編號，決定整張atlas有幾個window欄位
max_window_id_v6A<-max(snp_v6_del_y$window_id_v6A,na.rm=TRUE)

#確認共有幾條染色體：1–22、X、MT，預期為24
n_chr_v6A<-uniqueN(snp_v6_del_y$chr_order_v6A)

#計算整張atlas的長與寬
atlas_height_v6A<-n_chr_v6A*hilbert_side_v6A
atlas_width_v6A<-max_window_id_v6A*hilbert_side_v6A

max_window_id_v6A
n_chr_v6A
atlas_height_v6A
atlas_width_v6A

#染色體決定垂直位置，patch_row決定染色體區塊內的位置
snp_v6_del_y[,atlas_row_v6A:=
               (chr_order_v6A-1L)*hilbert_side_v6A+patch_row_v6A]

#window決定水平位置，patch_col決定window區塊內的位置
snp_v6_del_y[,atlas_col_v6A:=
               (window_id_v6A-1L)*hilbert_side_v6A+patch_col_v6A]

#記錄每條染色體在整張atlas中的row範圍
chr_atlas_map_v6A<-unique(
  snp_v6_del_y[,.(chr_order_v6A,Chr_id)]
)

chr_atlas_map_v6A[,`:=`(
  atlas_row_start_v6A=(chr_order_v6A-1L)*hilbert_side_v6A+1L,
  atlas_row_end_v6A=chr_order_v6A*hilbert_side_v6A
)]

setorder(chr_atlas_map_v6A,chr_order_v6A)
chr_atlas_map_v6A

snp_v6_del_y[
  Chr_id=="1"&window_id_v6A<=2,
  .(
    probeset_id,
    Start_num,
    window_id_v6A,
    hilbert_index_v6A,
    patch_row_v6A,
    patch_col_v6A,
    atlas_row_v6A,
    atlas_col_v6A
  )
][1:30]


#step6 nucleotide genotype 合併成七類----
#0 = padding（下一步建立）
#1 = AA/TT
#2 = CC/GG
#3 = AG/CT
#4 = AC/GT
#5 = AT
#6 = CG
#7 = NoCall


#抓取10位受試者的nucleotide genotype欄位
sample_nucleotide_cols_v6A<-grep("_nucleotide$",names(snp_v6_del_y),value=TRUE)

#整理受試者名稱
sample_ids_v6A<-sub("_nucleotide$","",sample_nucleotide_cols_v6A)
sample_ids_v6A<-sub("_\\(Axiom_TPM\\).*","",sample_ids_v6A)

#確認欄位數與受試者名稱
length(sample_nucleotide_cols_v6A)
data.table(
  sample_id=sample_ids_v6A,
  nucleotide_col=sample_nucleotide_cols_v6A
)

#將diploid nucleotide genotype合併成七類
collapse_nucleotide_v6A<-function(x){
  x<-toupper(gsub("\\s+","",as.character(x)))
  x<-gsub("\\|","/",x)
  
  fcase(
    x%chin%c("A/A","AA","T/T","TT"),"AA/TT",
    x%chin%c("C/C","CC","G/G","GG"),"CC/GG",
    x%chin%c("A/G","G/A","AG","GA","C/T","T/C","CT","TC"),"AG/CT",
    x%chin%c("A/C","C/A","AC","CA","G/T","T/G","GT","TG"),"AC/GT",
    x%chin%c("A/T","T/A","AT","TA"),"AT",
    x%chin%c("C/G","G/C","CG","GC"),"CG",
    default="NoCall"
  )
}

#設定七類genotype的數值編碼；0保留給padding
genotype_code_map_v6A<-c(
  "AA/TT"=1L,
  "CC/GG"=2L,
  "AG/CT"=3L,
  "AC/GT"=4L,
  "AT"=5L,
  "CG"=6L,
  "NoCall"=7L
)

#逐位受試者新增七類文字欄位與數值欄位
for(i in seq_along(sample_nucleotide_cols_v6A)){
  genotype7_now<-collapse_nucleotide_v6A(
    snp_v6_del_y[[sample_nucleotide_cols_v6A[i]]]
  )
  
  text_col_now<-paste0(sample_ids_v6A[i],"_genotype7_v6A")
  code_col_now<-paste0(sample_ids_v6A[i],"_genotype_code_v6A")
  
  set(snp_v6_del_y,j=text_col_now,value=genotype7_now)
  set(snp_v6_del_y,j=code_col_now,
      value=unname(genotype_code_map_v6A[genotype7_now]))
}

#記錄後續會使用的欄位名稱
genotype7_cols_v6A<-paste0(sample_ids_v6A,"_genotype7_v6A")
genotype_code_cols_v6A<-paste0(sample_ids_v6A,"_genotype_code_v6A")


genotype7_qc_v6A<-rbindlist(lapply(seq_along(sample_ids_v6A),function(i){
  snp_v6_del_y[,.(SNP_n=.N),
               by=.(genotype7=get(genotype7_cols_v6A[i]))][
                 ,sample_id:=sample_ids_v6A[i]
               ][,.(sample_id,genotype7,SNP_n)]
}))

setorder(genotype7_qc_v6A,sample_id,genotype7)
genotype7_qc_v6A

#step7:建立matrix----


#選擇第一位受試者
sample_index_v6A<-1L
sample_id_now_v6A<-sample_ids_v6A[sample_index_v6A]
code_col_now_v6A<-genotype_code_cols_v6A[sample_index_v6A]

#建立全為0的atlas；0代表padding
atlas_matrix_test_v6A<-matrix(as.raw(0L),
  nrow=atlas_height_v6A,
  ncol=atlas_width_v6A
)

#將每個真實SNP的genotype code填入共同座標
atlas_matrix_test_v6A[
  cbind(
    snp_v6_del_y$atlas_row_v6A,
    snp_v6_del_y$atlas_col_v6A
  )
]<-as.raw(snp_v6_del_y[[code_col_now_v6A]])

dim(atlas_matrix_test_v6A)
atlas_height_v6A
atlas_width_v6A

#非padding pixel數應等於SNP總數
real_snp_pixel_n_v6A<-sum(
  atlas_matrix_test_v6A!=as.raw(0L)
)

real_snp_pixel_n_v6A
nrow(snp_v6_del_y)

real_snp_pixel_n_v6A==nrow(snp_v6_del_y)

#padding 數量
#整張atlas總pixel數
total_pixel_n_v6A<-atlas_height_v6A*atlas_width_v6A

#沒有真實SNP的位置皆為padding
padding_pixel_n_v6A<-total_pixel_n_v6A-nrow(snp_v6_del_y)

data.table(
  total_pixel_n=total_pixel_n_v6A,
  real_SNP_pixel_n=nrow(snp_v6_del_y),
  padding_pixel_n=padding_pixel_n_v6A,
  padding_rate=padding_pixel_n_v6A/total_pixel_n_v6A
)

#依SNP座標從atlas中取回genotype code
code_readback_v6A<-as.integer(
  atlas_matrix_test_v6A[
    cbind(
      snp_v6_del_y$atlas_row_v6A,
      snp_v6_del_y$atlas_col_v6A
    )
  ]
)

#和原始genotype code比較
all(
  code_readback_v6A==
    snp_v6_del_y[[code_col_now_v6A]]
)

#原始SNP資料中的七類數量
table(snp_v6_del_y[[code_col_now_v6A]])

#從atlas真實SNP座標讀回後的七類數量
table(code_readback_v6A)

#Chr6的atlas row範圍
chr6_row_start_v6A<-(6L-1L)*hilbert_side_v6A+1L

chr6_row_end_v6A<-6L*hilbert_side_v6A

#Window 32的atlas column範圍
window32_col_start_v6A<-(32L-1L)*hilbert_side_v6A+1L

window32_col_end_v6A<-32L*hilbert_side_v6A

#擷取Chr6_W32的64×64 Hilbert patch
chr6_w32_patch_test_v6A<-atlas_matrix_test_v6A[
  chr6_row_start_v6A:chr6_row_end_v6A,
  window32_col_start_v6A:window32_col_end_v6A
]

dim(chr6_w32_patch_test_v6A) #64*64=4096(Chr6_W32 snp:3551 padding:545)
table(as.integer(chr6_w32_patch_test_v6A))


#step8----
#Chr6_W32


#指定要查看的染色體與window
view_chr_v6A<-"6"
view_window_v6A<-32L

#找出該染色體在atlas中的順序
view_chr_order_v6A<-unique(snp_v6_del_y[Chr_id==view_chr_v6A,chr_order_v6A])

#計算該patch在完整atlas中的row與column範圍
view_row_start_v6A<-(view_chr_order_v6A-1L)*hilbert_side_v6A+1L
view_row_end_v6A<-view_chr_order_v6A*hilbert_side_v6A
view_col_start_v6A<-(view_window_v6A-1L)*hilbert_side_v6A+1L
view_col_end_v6A<-view_window_v6A*hilbert_side_v6A

#擷取64×64 patch
patch_test_v6A<-atlas_matrix_test_v6A[
  view_row_start_v6A:view_row_end_v6A,
  view_col_start_v6A:view_col_end_v6A
]

dim(patch_test_v6A)

#查看0–7各類pixel數量
patch_code_count_v6A<-table(factor(as.integer(patch_test_v6A),levels=0:7))

patch_code_count_v6A

#真實SNP與padding數量
sum(patch_test_v6A!=as.raw(0L))
sum(patch_test_v6A==as.raw(0L))

#0為padding；1–7為七類genotype
genotype_color_v6A<-c(
  "white",     #0 Padding
  "#E41A1C",   #1 AA/TT
  "#377EB8",   #2 CC/GG
  "#4DAF4A",   #3 AG/CT
  "#984EA3",   #4 AC/GT
  "#FF7F00",   #5 AT
  "#A65628",   #6 CG
  "grey40"     #7 NoCall
)

genotype_label_v6A<-c(
  "Padding",
  "AA/TT",
  "CC/GG",
  "AG/CT",
  "AC/GT",
  "AT",
  "CG",
  "NoCall"
)

#轉為一般integer matrix
patch_plot_v6A<-matrix(
  as.integer(patch_test_v6A),
  nrow=hilbert_side_v6A,
  ncol=hilbert_side_v6A
)

#轉置並上下翻轉，使matrix第1列顯示在圖片上方
patch_plot_display_v6A<-t(
  patch_plot_v6A[nrow(patch_plot_v6A):1,]
)

#畫出64×64 Hilbert genotype patch
par(mar=c(2,2,4,8),xpd=TRUE)
image(
  x=1:hilbert_side_v6A,
  y=1:hilbert_side_v6A,
  z=patch_plot_display_v6A,
  col=genotype_color_v6A,
  breaks=seq(-0.5,7.5,by=1),
  axes=FALSE,
  asp=1,
  useRaster=TRUE,
  main=paste0(sample_id_now_v6A,"－Chr6 Window 32")
)
box()
legend(
  "right",
  inset=c(-0.38,0),
  legend=genotype_label_v6A,
  fill=genotype_color_v6A,
  border=NA,
  bty="n",
  cex=0.85
)

snp_v6_del_y[
  Chr_id=="6"&window_id_v6A==32,
  .(
    probeset_id,
    Start_num,
    snp_order_within_window_v6A,
    hilbert_index_v6A,
    patch_row_v6A,
    patch_col_v6A,
    genotype=get(genotype7_cols_v6A[1]),
    genotype_code=get(genotype_code_cols_v6A[1])
  )
][1:20]

chr6_w32_order_check_v6A<-snp_v6_del_y[
  Chr_id=="6"&window_id_v6A==32
]

all(diff(chr6_w32_order_check_v6A$Start_num)>=0)
all(
  chr6_w32_order_check_v6A$hilbert_index_v6A==
    0:(nrow(chr6_w32_order_check_v6A)-1L)
)

#step9 10人的matrix----

#所有SNP共同使用的atlas座標
atlas_pixel_index_v6A<-cbind(
  snp_v6_del_y$atlas_row_v6A,
  snp_v6_del_y$atlas_col_v6A
)

#建立存放10張atlas的list
atlas_matrix_all_v6A<-setNames(vector("list",length(sample_ids_v6A)),sample_ids_v6A)


#逐位受試者建立raw matrix
for(i in seq_along(sample_ids_v6A)){
  atlas_now<-matrix(
    as.raw(0L),
    nrow=atlas_height_v6A,
    ncol=atlas_width_v6A
  )
  
  atlas_now[atlas_pixel_index_v6A]<-as.raw(snp_v6_del_y[[genotype_code_cols_v6A[i]]])
  
  atlas_matrix_all_v6A[[i]]<-atlas_now
  cat(i,"/",length(sample_ids_v6A),sample_ids_v6A[i],"完成\n")
}

#每張atlas 7種genotype數量
atlas_code_count_v6A<-rbindlist(lapply(seq_along(sample_ids_v6A),function(i){
  code_count<-table(
    factor(
      as.integer(atlas_matrix_all_v6A[[i]]),
      levels=0:7
    )
  )
  
  data.table(
    sample_id=sample_ids_v6A[i],
    genotype_code=0:7,
    pixel_n=as.integer(code_count)
  )
}))

atlas_code_count_v6A
#add genotype_class name
atlas_code_count_v6A[,genotype_class:=genotype_label_v6A[genotype_code+1L]]
setcolorder(atlas_code_count_v6A,c("sample_id","genotype_code","genotype_class","pixel_n"))

saveRDS(atlas_matrix_all_v6A,file="atlas_matrix_all_v6A.rds",compress="xz")

#step10----
#Version 6A Step 10：繪製單一受試者完整atlas（最終修正版）----------------

#選擇第一位受試者
view_sample_id_v6A<-sample_ids_v6A[1]

#取出該受試者的raw atlas matrix
atlas_view_raw_v6A<-atlas_matrix_all_v6A[[view_sample_id_v6A]]

#將raw matrix轉成0–7的integer matrix
atlas_view_plot_v6A<-matrix(
  as.integer(atlas_view_raw_v6A),
  nrow=atlas_height_v6A,
  ncol=atlas_width_v6A
)

#設定padding與七類genotype顏色
genotype_color_v6A<-c(
  "white",     #0 Padding
  "#E41A1C",   #1 AA/TT
  "#377EB8",   #2 CC/GG
  "#4DAF4A",   #3 AG/CT
  "#984EA3",   #4 AC/GT
  "#FF7F00",   #5 AT
  "#A65628",   #6 CG
  "grey40"     #7 NoCall
)

#設定圖例名稱
genotype_label_v6A<-c(
  "Padding",
  "AA/TT",
  "CC/GG",
  "AG/CT",
  "AC/GT",
  "AT",
  "CG",
  "NoCall"
)

#建立染色體順序與標籤位置
chr_label_v6A<-unique(
  snp_v6_del_y[,.(chr_order_v6A,Chr_id)]
)

setorder(chr_label_v6A,chr_order_v6A)

#因為Chr1放在圖片最上方，所以需要反向計算Y軸位置
chr_label_v6A[,display_y_v6A:=
                atlas_height_v6A-
                (chr_order_v6A-0.5)*hilbert_side_v6A
]

#建立X軸window刻度：第1個window後，每20 Mb標示一次
window_tick_id_v6A<-sort(unique(c(
  1L,
  seq(20L,max_window_id_v6A,by=20L)
)))

#若最大window不是20的倍數，將最後一個window也加入
window_tick_id_v6A<-sort(unique(c(
  window_tick_id_v6A,
  max_window_id_v6A
)))

#將window編號換算成atlas中的pixel位置
window_tick_x_v6A<-
  (window_tick_id_v6A-0.5)*hilbert_side_v6A

#輸出檔案位置
atlas_png_file_v6A<-file.path(
  getwd(),
  paste0(
    view_sample_id_v6A,
    "_Version6A_full_atlas_v3.png"
  )
)

#設定圖片解析度
atlas_png_res_v6A<-120

#設定四周留白，單位為inch
plot_margin_bottom_v6A<-0.55
plot_margin_left_v6A<-0.80
plot_margin_top_v6A<-1.10
plot_margin_right_v6A<-0.15

#依atlas實際pixel數計算輸出圖片大小
atlas_png_width_v6A<-as.integer(
  atlas_width_v6A+
    (plot_margin_left_v6A+plot_margin_right_v6A)*atlas_png_res_v6A
)

atlas_png_height_v6A<-as.integer(
  atlas_height_v6A+
    (plot_margin_bottom_v6A+plot_margin_top_v6A)*atlas_png_res_v6A
)

#開啟PNG繪圖裝置
png(
  filename=atlas_png_file_v6A,
  width=atlas_png_width_v6A,
  height=atlas_png_height_v6A,
  units="px",
  res=atlas_png_res_v6A
)

#固定四周留白，xpd=NA允許圖例畫在黑框上方
par(
  mai=c(
    plot_margin_bottom_v6A,
    plot_margin_left_v6A,
    plot_margin_top_v6A,
    plot_margin_right_v6A
  ),
  xpd=NA
)

#將matrix上下反轉並轉置
#使Chr1顯示於最上方，ChrMT顯示於最下方
atlas_display_v6A<-t(
  atlas_view_plot_v6A[
    atlas_height_v6A:1L,
    ,
    drop=FALSE
  ]
)

#繪製完整atlas
image(
  x=seq(0.5,atlas_width_v6A-0.5,by=1),
  y=seq(0.5,atlas_height_v6A-0.5,by=1),
  z=atlas_display_v6A,
  col=genotype_color_v6A,
  breaks=seq(-0.5,7.5,by=1),
  xlim=c(0,atlas_width_v6A),
  ylim=c(0,atlas_height_v6A),
  axes=FALSE,
  xlab="",
  ylab="",
  xaxs="i",
  yaxs="i",
  asp=1,
  useRaster=TRUE
)

#畫每條染色體之間的水平分隔線
abline(
  h=seq(
    hilbert_side_v6A,
    atlas_height_v6A-hilbert_side_v6A,
    by=hilbert_side_v6A
  ),
  col="grey65",
  lwd=0.4
)

#畫每個1 Mb window之間的垂直分隔線
abline(
  v=seq(
    hilbert_side_v6A,
    atlas_width_v6A-hilbert_side_v6A,
    by=hilbert_side_v6A
  ),
  col="grey88",
  lwd=0.2
)

#畫X軸，只顯示window編號
axis(
  side=1,
  at=window_tick_x_v6A,
  labels=window_tick_id_v6A,
  cex.axis=0.65,
  las=1,
  lwd=0.7,
  lwd.ticks=0.7
)

#畫Y軸，只顯示24條染色體名稱
axis(
  side=2,
  at=chr_label_v6A$display_y_v6A,
  labels=paste0("Chr",chr_label_v6A$Chr_id),
  las=1,
  cex.axis=0.65,
  lwd=0.7,
  lwd.ticks=0.7
)

#加入X軸與Y軸標題
mtext(
  "1-Mb genomic window",
  side=1,
  line=2.3,
  cex=0.8
)

mtext(
  "Chromosome",
  side=2,
  line=3.5,
  cex=0.8
)

#畫atlas黑色外框
rect(
  xleft=0,
  ybottom=0,
  xright=atlas_width_v6A,
  ytop=atlas_height_v6A,
  border="black",
  lwd=1
)

#橫向圖例，放在黑框上方
legend(
  x=atlas_width_v6A/2,
  y=atlas_height_v6A+hilbert_side_v6A*0.40,
  legend=genotype_label_v6A,
  fill=genotype_color_v6A,
  border=NA,
  bty="n",
  horiz=TRUE,
  xjust=0.5,
  yjust=0.5,
  cex=0.75,
  xpd=NA
)

#加入主標題
mtext(
  paste0(
    view_sample_id_v6A,
    " - Version 6A chromosome × window Hilbert atlas"
  ),
  side=3,
  line=4.2,
  cex=1,
  font=2
)

#完成並關閉圖片
dev.off()

#確認圖片成功輸出
atlas_png_file_v6A
file.exists(atlas_png_file_v6A)
file.info(atlas_png_file_v6A)$size/1024^2

#Mac直接開啟圖片
system(paste("open",shQuote(atlas_png_file_v6A)))

#刪除僅供繪圖使用的大型matrix
rm(
  atlas_view_plot_v6A,
  atlas_display_v6A
)

gc()

#step11. `10個人----


#建立輸出資料夾
atlas_output_dir_v6A<-"Version6A_full_atlas_v3"
dir.create(atlas_output_dir_v6A,showWarnings=FALSE,recursive=TRUE)

#建立染色體標籤位置
chr_label_v6A<-unique(
  snp_v6_del_y[,.(chr_order_v6A,Chr_id)]
)
setorder(chr_label_v6A,chr_order_v6A)

#Chr1顯示在上方，ChrMT顯示在下方
chr_label_v6A[,display_y_v6A:=
                atlas_height_v6A-
                (chr_order_v6A-0.5)*hilbert_side_v6A
]

#建立X軸刻度：Window 1、每20個window，以及最後一個window
window_tick_id_v6A<-sort(unique(c(
  1L,
  seq(20L,max_window_id_v6A,by=20L),
  max_window_id_v6A
)))

#將window編號轉成atlas pixel位置
window_tick_x_v6A<-
  (window_tick_id_v6A-0.5)*hilbert_side_v6A

#設定輸出解析度與四周留白
atlas_png_res_v6A<-120
plot_margin_bottom_v6A<-0.55
plot_margin_left_v6A<-0.80
plot_margin_top_v6A<-1.10
plot_margin_right_v6A<-0.15

#依atlas實際尺寸決定PNG大小
atlas_png_width_v6A<-as.integer(
  atlas_width_v6A+
    (plot_margin_left_v6A+plot_margin_right_v6A)*atlas_png_res_v6A
)

atlas_png_height_v6A<-as.integer(
  atlas_height_v6A+
    (plot_margin_bottom_v6A+plot_margin_top_v6A)*atlas_png_res_v6A
)
#建立單一受試者完整atlas繪圖函數
plot_full_atlas_v6A<-function(sample_id,atlas_raw){
  
  #將raw matrix轉成0–7 integer matrix
  atlas_plot<-matrix(
    as.integer(atlas_raw),
    nrow=atlas_height_v6A,
    ncol=atlas_width_v6A
  )
  
  #上下反轉並轉置，使Chr1在最上方、ChrMT在最下方
  atlas_display<-t(
    atlas_plot[
      atlas_height_v6A:1L,
      ,
      drop=FALSE
    ]
  )
  
  #設定輸出檔名
  output_file<-file.path(
    atlas_output_dir_v6A,
    paste0(sample_id,"_Version6A_full_atlas_v3.png")
  )
  
  #開啟PNG繪圖裝置
  png(
    filename=output_file,
    width=atlas_png_width_v6A,
    height=atlas_png_height_v6A,
    units="px",
    res=atlas_png_res_v6A
  )
  
  #設定四周留白，允許圖例畫在黑框外
  par(
    mai=c(
      plot_margin_bottom_v6A,
      plot_margin_left_v6A,
      plot_margin_top_v6A,
      plot_margin_right_v6A
    ),
    xpd=NA
  )
  
  #繪製完整atlas
  image(
    x=seq(0.5,atlas_width_v6A-0.5,by=1),
    y=seq(0.5,atlas_height_v6A-0.5,by=1),
    z=atlas_display,
    col=genotype_color_v6A,
    breaks=seq(-0.5,7.5,by=1),
    xlim=c(0,atlas_width_v6A),
    ylim=c(0,atlas_height_v6A),
    axes=FALSE,
    xlab="",
    ylab="",
    xaxs="i",
    yaxs="i",
    asp=1,
    useRaster=TRUE
  )
  
  #染色體之間的水平分隔線
  abline(
    h=seq(
      hilbert_side_v6A,
      atlas_height_v6A-hilbert_side_v6A,
      by=hilbert_side_v6A
    ),
    col="grey65",
    lwd=0.4
  )
  
  #1 Mb window之間的垂直分隔線
  abline(
    v=seq(
      hilbert_side_v6A,
      atlas_width_v6A-hilbert_side_v6A,
      by=hilbert_side_v6A
    ),
    col="grey88",
    lwd=0.2
  )
  
  #X軸：1 Mb window編號
  axis(
    side=1,
    at=window_tick_x_v6A,
    labels=window_tick_id_v6A,
    cex.axis=0.65,
    las=1,
    lwd=0.7,
    lwd.ticks=0.7
  )
  
  #Y軸：染色體名稱
  axis(
    side=2,
    at=chr_label_v6A$display_y_v6A,
    labels=paste0("Chr",chr_label_v6A$Chr_id),
    las=1,
    cex.axis=0.65,
    lwd=0.7,
    lwd.ticks=0.7
  )
  
  #座標軸標題
  mtext(
    "1-Mb genomic window",
    side=1,
    line=2.3,
    cex=0.8
  )
  
  mtext(
    "Chromosome",
    side=2,
    line=3.5,
    cex=0.8
  )
  
  #完整atlas黑色外框
  rect(
    xleft=0,
    ybottom=0,
    xright=atlas_width_v6A,
    ytop=atlas_height_v6A,
    border="black",
    lwd=1
  )
  
  #橫向圖例放在黑框上方
  legend(
    x=atlas_width_v6A/2,
    y=atlas_height_v6A+hilbert_side_v6A*0.40,
    legend=genotype_label_v6A,
    fill=genotype_color_v6A,
    border=NA,
    bty="n",
    horiz=TRUE,
    xjust=0.5,
    yjust=0.5,
    cex=0.75,
    xpd=NA
  )
  
  #主標題
  mtext(
    paste0(
      sample_id,
      " - Version 6A chromosome × window Hilbert atlas"
    ),
    side=3,
    line=4.2,
    cex=1,
    font=2
  )
  
  #關閉PNG裝置
  dev.off()
  
  #刪除本次繪圖暫存matrix
  rm(atlas_plot,atlas_display)
  invisible(gc())
  
  return(output_file)
}
#建立儲存輸出檔案路徑的向量
atlas_png_files_v6A<-setNames(
  character(length(sample_ids_v6A)),
  sample_ids_v6A
)

#逐位受試者輸出完整atlas
for(i in seq_along(sample_ids_v6A)){
  
  sample_id_now<-sample_ids_v6A[i]
  
  atlas_png_files_v6A[i]<-plot_full_atlas_v6A(
    sample_id=sample_id_now,
    atlas_raw=atlas_matrix_all_v6A[[sample_id_now]]
  )
  
  cat(
    i,"/",length(sample_ids_v6A),
    sample_id_now,
    "完成\n"
  )
}
#整理10張圖片的輸出資訊
atlas_png_qc_v6A<-data.table(
  sample_id=sample_ids_v6A,
  file_path=unname(atlas_png_files_v6A),
  file_exists=file.exists(atlas_png_files_v6A),
  file_size_MB=round(
    file.info(atlas_png_files_v6A)$size/1024^2,
    3
  )
)

atlas_png_qc_v6A

system(paste("open",shQuote(atlas_output_dir_v6A)))



#step12 SNP-to-pixel 對照表----

#建立真實SNP的pixel對照表；padding不屬於任何SNP，因此不會出現在這張表
pixel_snp_map_v6A<-snp_v6_del_y[,.(
  
  #原始SNP資訊
  probeset_id,
  dbSNP_RS_ID,
  Chr_id,
  Start=Start_num,
  Allele_A,
  Allele_B,
  
  #1 Mb window資訊
  chr_order_v6A,
  window_id_v6A,
  window_key_v6A,
  window_bp_start_v6A,
  window_bp_end_v6A,
  
  #window內排列資訊
  snp_order_within_window_v6A,
  hilbert_index_v6A,
  patch_row_v6A,
  patch_col_v6A,
  
  #完整atlas座標
  atlas_row_v6A,
  atlas_col_v6A
)]

#固定對照表順序
setorder(
  pixel_snp_map_v6A,
  chr_order_v6A,
  window_id_v6A,
  snp_order_within_window_v6A
)

#查看前20筆
pixel_snp_map_v6A[1:20]
saveRDS(
  pixel_snp_map_v6A,
  file="pixel_snp_map_v6A.rds",
  compress="xz"
)




#Version 6A：視覺化Hilbert curve路徑順序----------------


#複製Hilbert template；轉成繪圖座標
hilbert_path_plot_v6A<-copy(hilbert_template_v6A)
hilbert_path_plot_v6A[,plot_x:=patch_col_v6A]
hilbert_path_plot_v6A[,plot_y:=hilbert_side_v6A-patch_row_v6A+1L]

#前64步
n_show_v6A<-64L
first_path_v6A<-hilbert_path_plot_v6A[1:n_show_v6A]

#起點終點座標
start_x<-hilbert_path_plot_v6A$plot_x[1]
start_y<-hilbert_path_plot_v6A$plot_y[1]
end_x<-hilbert_path_plot_v6A$plot_x[nrow(hilbert_path_plot_v6A)]
end_y<-hilbert_path_plot_v6A$plot_y[nrow(hilbert_path_plot_v6A)]

png("Hilbert_path_v6A_fix.png",width=1800,height=900,res=150)
par(mfrow=c(1,2),mar=c(4,4,4,2))

#左圖：完整Hilbert path
plot(hilbert_path_plot_v6A$plot_x,hilbert_path_plot_v6A$plot_y,
     type="n",
     xlim=c(1,hilbert_side_v6A+4),
     ylim=c(1,hilbert_side_v6A+2),
     xlab="Column",
     ylab="Row",
     main=paste0("Version 6A Hilbert path (",hilbert_side_v6A,"×",hilbert_side_v6A,")"),
     asp=1,
     xaxt="n",
     yaxt="n")

abline(v=seq(0.5,hilbert_side_v6A+0.5,by=1),col="grey92",lwd=0.3)
abline(h=seq(0.5,hilbert_side_v6A+0.5,by=1),col="grey92",lwd=0.3)
axis(1,at=seq(1,hilbert_side_v6A,by=8))
axis(2,at=seq(1,hilbert_side_v6A,by=8),las=1)

lines(hilbert_path_plot_v6A$plot_x,hilbert_path_plot_v6A$plot_y,col="grey35",lwd=1)
points(hilbert_path_plot_v6A$plot_x,hilbert_path_plot_v6A$plot_y,pch=16,cex=0.12,col="grey35")

#起點終點
points(start_x,start_y,pch=16,cex=1.6,col="red")
points(end_x,end_y,pch=16,cex=1.6,col="blue")

text(start_x+1.5,start_y+0.5,labels="Start",col="red",cex=0.9,xpd=NA,adj=c(0,0.5))
text(end_x-1.5,end_y+0.5,labels="End",col="blue",cex=0.9,xpd=NA,adj=c(1,0.5))
box()

#右圖：前64步放大
plot(first_path_v6A$plot_x,first_path_v6A$plot_y,
     type="n",
     xlim=c(min(first_path_v6A$plot_x)-1,max(first_path_v6A$plot_x)+2),
     ylim=c(min(first_path_v6A$plot_y)-1,max(first_path_v6A$plot_y)+1.5),
     xlab="Column",
     ylab="Row",
     main=paste0("First ",n_show_v6A," Hilbert steps"),
     asp=1)

abline(v=seq(floor(min(first_path_v6A$plot_x))-1,ceiling(max(first_path_v6A$plot_x))+1,by=1),col="grey92",lwd=0.3)
abline(h=seq(floor(min(first_path_v6A$plot_y))-1,ceiling(max(first_path_v6A$plot_y))+1,by=1),col="grey92",lwd=0.3)

lines(first_path_v6A$plot_x,first_path_v6A$plot_y,col="grey35",lwd=1.2)
points(first_path_v6A$plot_x,first_path_v6A$plot_y,pch=16,cex=0.8,col="grey35")

#每一步的index
text(first_path_v6A$plot_x,first_path_v6A$plot_y,
     labels=first_path_v6A$hilbert_index_v6A,
     pos=3,
     cex=0.5)

#起點終點
points(first_path_v6A$plot_x[1],first_path_v6A$plot_y[1],pch=16,cex=1.6,col="red")
points(first_path_v6A$plot_x[nrow(first_path_v6A)],first_path_v6A$plot_y[nrow(first_path_v6A)],pch=16,cex=1.6,col="blue")

text(first_path_v6A$plot_x[1]-0.4,first_path_v6A$plot_y[1]-0.2,labels="Start",col="red",cex=0.9,xpd=NA,adj=c(0,1))
text(first_path_v6A$plot_x[nrow(first_path_v6A)]-0.4,first_path_v6A$plot_y[nrow(first_path_v6A)]-0.2,
     labels="End",col="blue",cex=0.9,xpd=NA,adj=c(0,1))
box()
dev.off()

system("open Hilbert_path_v6A_fix.png")

#Version 6A Step 13：視覺化實際window在Hilbert curve中的填入方式----------------

#選擇SNP最多的window：Chr6_W32
view_window_key_v6A<-"Chr6_W32"

#抓出該window的SNP資料
view_window_snp_v6A<-snp_v6_del_y[
  window_key_v6A==view_window_key_v6A
]

#確認SNP數
nrow(view_window_snp_v6A)

#複製完整64×64 Hilbert template
hilbert_window_view_v6A<-copy(hilbert_template_v6A)

#0-based Hilbert座標轉成畫圖座標
hilbert_window_view_v6A[,plot_x:=patch_col_v6A]
hilbert_window_view_v6A[,plot_y:=hilbert_side_v6A-patch_row_v6A+1L]

#標示每個Hilbert位置是真實SNP還是padding
hilbert_window_view_v6A[,pixel_type:=
                          ifelse(
                            hilbert_index_v6A<nrow(view_window_snp_v6A),
                            "SNP",
                            "Padding"
                          )
]

table(hilbert_window_view_v6A$pixel_type)


#輸出圖片
png(
  "Chr6_W32_Hilbert_fill_v6A.png",
  width=1000,
  height=1000,
  res=150
)

par(mar=c(4,4,4,2))

#建立空白64×64畫布
plot(
  hilbert_window_view_v6A$plot_x,
  hilbert_window_view_v6A$plot_y,
  type="n",
  xlim=c(0.5,64.5),
  ylim=c(0.5,64.5),
  xlab="Column",
  ylab="Row",
  main="Chr6_W32: SNP filling along the 64×64 Hilbert curve",
  asp=1,
  xaxt="n",
  yaxt="n"
)

#格線
abline(
  v=seq(0.5,64.5,by=1),
  h=seq(0.5,64.5,by=1),
  col="grey93",
  lwd=0.3
)

axis(1,at=seq(1,64,by=8))
axis(2,at=seq(1,64,by=8),las=1)

#完整Hilbert路徑
lines(
  hilbert_window_view_v6A$plot_x,
  hilbert_window_view_v6A$plot_y,
  col="grey75",
  lwd=0.8
)

#真實SNP位置
points(
  hilbert_window_view_v6A[pixel_type=="SNP",plot_x],
  hilbert_window_view_v6A[pixel_type=="SNP",plot_y],
  pch=16,
  cex=0.25,
  col="black"
)

#padding位置
points(
  hilbert_window_view_v6A[pixel_type=="Padding",plot_x],
  hilbert_window_view_v6A[pixel_type=="Padding",plot_y],
  pch=16,
  cex=0.25,
  col="grey80"
)

#起點
points(
  hilbert_window_view_v6A$plot_x[1],
  hilbert_window_view_v6A$plot_y[1],
  pch=16,
  cex=1.3,
  col="red"
)

text(
  hilbert_window_view_v6A$plot_x[1]+1,
  hilbert_window_view_v6A$plot_y[1],
  "Start",
  col="red",
  pos=4
)

#最後一個真實SNP
last_snp_index_v6A<-nrow(view_window_snp_v6A)

points(
  hilbert_window_view_v6A$plot_x[last_snp_index_v6A],
  hilbert_window_view_v6A$plot_y[last_snp_index_v6A],
  pch=16,
  cex=1.3,
  col="blue"
)

text(
  hilbert_window_view_v6A$plot_x[last_snp_index_v6A]-1,
  hilbert_window_view_v6A$plot_y[last_snp_index_v6A],
  "Last SNP",
  col="blue",
  pos=2
)

legend(
  "bottomright",
  legend=c("SNP","Padding","Start","Last SNP"),
  pch=16,
  col=c("black","grey80","red","blue"),
  bty="n"
)

box()
dev.off()

#Mac直接開啟
system("open Chr6_W32_Hilbert_fill_v6A.png")



#Version 6B--------------
#Version 6B Step 1：建立 chromosome sector 基礎表,確認每條染色體在圓形 atlas 中應該占多少個 1 Mb 單位
##整理每條染色體的SNP數、最大window與實際涵蓋位置
chr_sector_v6b<-snp_v6_del_y[,.(SNP_n=.N,
                                max_window_v6b=max(window_id_v6A),
                                first_snp_bp=min(Start_num),
                                last_snp_bp=max(Start_num)),#染色體越長、涵蓋的 1 Mb windows 越多，在圓周上占的 sector 就越長
                                by=.(chr_order_v6A,Chr_id)]

#order by chr1-22 x mt
setorder(chr_sector_v6b,chr_order_v6A)
#rename
setnames(chr_sector_v6b,"chr_order_v6A","chr_order_v6b")

#把整個 genome 分配到 360°,所以計算all chr 共有多少mb->3043
#不能代表實際有snp的window
total_window_v6b<-sum(chr_sector_v6b$max_window_v6b)

#check 實際有snp的window
observe_window_v6b<-unique(snp_v6_del_y[,.(Chr_id,window_id_v6A)
                                        ]
                          )[,.(observe_window_n=.N),
                            by=Chr_id]
#合併到chr sector
chr_sector_v6b[observe_window_v6b,observe_window_n:=i.observe_window_n,on="Chr_id"]



#Step 2：加入 chromosome gap，計算每條染色體在圓上的起始與終止角度
#Chr1 從 12 點鐘方向開始 每兩條染色體間留 1° gap
#所以一個 1 Mb window 在所有染色體上都具有相同角度

#每條染色體之間保留1度gap
chr_gap_deg_v6b<-1

#共有24條染色體，因此整圈共有24個gap
n_chr_v6b<-nrow(chr_sector_v6b)
total_gap_deg_v6b<-n_chr_v6b*chr_gap_deg_v6b

#剩餘角度用來放真正的genomic windows
genome_available_deg_v6b<-360-total_gap_deg_v6b


#計算每1 Mb window固定占多少角度->0.1104174
angle_per_window_v6b<-genome_available_deg_v6b/total_window_v6b

#每條染色體需要的角度 = window數 × 每個window固定角度
chr_sector_v6b[,sector_span_deg_v6b:=max_window_v6b*angle_per_window_v6b]

#依序累積前面染色體sector與gap，得到每條染色體起始角度
chr_sector_v6b[,sector_start_clock_deg_v6b:=cumsum(shift(sector_span_deg_v6b+chr_gap_deg_v6b,fill=0))]
#染色體sector結束角度
chr_sector_v6b[,sector_end_clock_deg_v6b:=sector_start_clock_deg_v6b+sector_span_deg_v6b]
#該染色體後gap的範圍
chr_sector_v6b[,gap_start_clock_deg_v6b:=sector_end_clock_deg_v6b]
chr_sector_v6b[,gap_end_clock_deg_v6b:=sector_end_clock_deg_v6b+chr_gap_deg_v6b]

chr_sector_v6b[,.(chr_order_v6b,Chr_id,max_window_v6b,sector_span_deg_v6b,sector_start_clock_deg_v6b,sector_end_clock_deg_v6b,gap_start_clock_deg_v6b,gap_end_clock_deg_v6b)]


#建立弧度
#x = r × cos(theta)
#y = r × sin(theta)

#將12點鐘起始、順時針角度轉成R標準radian
chr_sector_v6b[,sector_start_theta_v6b:=
                 pi/2-sector_start_clock_deg_v6b*pi/180]

chr_sector_v6b[,sector_end_theta_v6b:=
                 pi/2-sector_end_clock_deg_v6b*pi/180]

#Step 3：把每條染色體的 1 Mb windows 映射到固定角度
#每個 1 Mb window 在整個圓上都占完全相同的角度

#建立每條染色體完整的physical window
window_map_v6b<-rbindlist(
  lapply(seq_len(nrow(chr_sector_v6b)),function(i){
    
    data.table(
      chr_order_v6b=chr_sector_v6b$chr_order_v6b[i],
      Chr_id=chr_sector_v6b$Chr_id[i],
      window_id_v6A=seq_len(chr_sector_v6b$max_window_v6b[i]),
      sector_start_clock_deg_v6b=chr_sector_v6b$sector_start_clock_deg_v6b[i]
    )
    
  })
)

#建立window名稱
window_map_v6b[,window_key_v6b:=
                 paste0("Chr",Chr_id,"_W",window_id_v6A)]

setorder(window_map_v6b,chr_order_v6b,window_id_v6A)

#step3-2計算每個 window 的起點與終點角度
#每 1 Mb 占幾度:angle_per_window_v6b
#每個1 Mb window的順時針起始角度
window_map_v6b[,window_start_clock_deg_v6b:=sector_start_clock_deg_v6b+(window_id_v6A-1L)*angle_per_window_v6b]

#每個window的終點角度
window_map_v6b[,window_end_clock_deg_v6b:= window_start_clock_deg_v6b+angle_per_window_v6b]

#window中央角度，之後可用來標示位置
window_map_v6b[,window_mid_clock_deg_v6b:=(window_start_clock_deg_v6b+window_end_clock_deg_v6b)/2]

#step3-3
#轉成R使用的radian
window_map_v6b[,window_start_theta_v6b:=
                 pi/2-window_start_clock_deg_v6b*pi/180]

window_map_v6b[,window_end_theta_v6b:=
                 pi/2-window_end_clock_deg_v6b*pi/180]

window_map_v6b[,window_mid_theta_v6b:=
                 pi/2-window_mid_clock_deg_v6b*pi/180]

#Step 3-4：加入每個 window 實際 SNP 數
#SNP_n > 0 → 真實有SNP的window
#SNP_n = 0 → physical window存在，但array沒有SNP

#把每個window實際SNP數加入circular window map
window_map_v6b[
  window_snp_count_v6A,
  SNP_n:=i.SNP_n,
  on=.(Chr_id,window_id_v6A)
]

#沒有SNP的physical window設為0
window_map_v6b[is.na(SNP_n),SNP_n:=0L]


window_map_v6b[
  1:20,
  .(
    Chr_id,
    window_id_v6A,
    window_key_v6b,
    SNP_n,
    window_start_clock_deg_v6b,
    window_end_clock_deg_v6b,
    window_mid_clock_deg_v6b
  )
]




#step4 如何排列
#角度方向固定 64 格；不足的 SNP 往下一個 radial row 繼續排
#一蘭最多3552個snp 排序下一個換行 最多需要888行

#step4每個1 Mb window沿角度方向切成4格
angular_slots_v6b<-4L
#根據最大window SNP數決定需要多少radial rows->888
radial_slots_v6b<-ceiling(max_snp_per_window_v6A/angular_slots_v6b)

#每個window總容量->3552
window_capacity_v6b<-angular_slots_v6b*radial_slots_v6b

#Step 4-1：每個 window 內重新依 Start 排序
#固定排序：染色體→window→Start→probeset_id
setorder(
  snp_v6_del_y,
  chr_order_v6A,
  window_id_v6A,
  Start_num,
  probeset_id
)
#Version6B window內SNP順序
snp_v6_del_y[,snp_order_within_window_v6b:=
               seq_len(.N),
             by=.(Chr_id,window_id_v6A)]

#Step 4-2：把 SNP 放入 64 × 56 的 sector grid
#蛇行排列 ：一行滿了換下一行
#radial row 1：1 → 2 → 3 → ... → 64
#radial row 2：128 ← 127 ← ... ← 65

#先建立0-based sequence index
snp_v6_del_y[,sector_index_v6b:=snp_order_within_window_v6b-1L]

#計算位於第幾個radial row，0-based
snp_v6_del_y[,radial_index_v6b:=sector_index_v6b%/%angular_slots_v6b]

#計算原始angular位置，0–63
snp_v6_del_y[,angular_index_raw_v6b:=sector_index_v6b%%angular_slots_v6b]

#snake排列：
#偶數radial row由左到右
#奇數radial row由右到左
snp_v6_del_y[,angular_index_v6b:=
               fifelse(
                 radial_index_v6b%%2L==0L,
                 angular_index_raw_v6b,
                 angular_slots_v6b-1L-angular_index_raw_v6b
               )]
#Step 4-3：計算每個 SNP 在自己 window 中的角度
#把Version6B window角度資訊加入每個SNP
snp_v6_del_y[
  window_map_v6b,
  `:=`(
    window_start_clock_deg_v6b=i.window_start_clock_deg_v6b,
    window_end_clock_deg_v6b=i.window_end_clock_deg_v6b
  ),
  on=.(Chr_id,window_id_v6A)
]

#每個window內，每個angular slot占多少度
angle_per_snp_slot_v6b<-angle_per_window_v6b/angular_slots_v6b

#每個SNP所在angular cell的中央角度
snp_v6_del_y[,snp_clock_deg_v6b:=window_start_clock_deg_v6b+(angular_index_v6b+0.5)*angle_per_snp_slot_v6b]

#轉成R使用的radian
snp_v6_del_y[,snp_theta_v6b:=
               pi/2-snp_clock_deg_v6b*pi/180]

#Step 4-4：目前 radial 座標先保留為 index
#轉為1-based，方便之後解釋
snp_v6_del_y[,radial_slot_v6b:=
               radial_index_v6b+1L]

#angular位置也轉成1-based
snp_v6_del_y[,angular_slot_v6b:=
               angular_index_v6b+1L]
#check
snp_v6_del_y[
  window_key_v6A=="Chr6_W32",
  .(
    probeset_id,
    Start_num,
    snp_order_within_window_v6b,
    radial_slot_v6b,
    angular_slot_v6b,
    snp_clock_deg_v6b
  )
][1:70]

#Step 5：建立polar coordinate----

#每個radial slot暫定1 pixel厚
radial_step_px_v6b<-1

#整個genotype track需要的徑向厚度
genotype_track_width_v6b<-radial_slots_v6b*radial_step_px_v6b

#每個SNP angular slot所占的弧度
angle_per_snp_slot_rad_v6b<-angle_per_snp_slot_v6b*pi/180

#讓最內圈每個angular slot至少約1 pixel寬->2076
target_arc_px_v6b<-1
inner_radius_v6b<-ceiling(target_arc_px_v6b/angle_per_snp_slot_rad_v6b)

#Step 5-3：決定 outer radius外半徑
#外半徑 = 內半徑 + 888個radial slots
#中心先保留空白 未來可放其他東西
outer_radius_v6b<-inner_radius_v6b+genotype_track_width_v6b

#Step 5-4：為每個 SNP 建立實際 radius
#radial slot 1   → 最靠內
#radial slot 888 → 最靠外
#每個 SNP 放在該 radial cell 的中央
#每個SNP的實際radius
snp_v6_del_y[,snp_radius_v6b:=inner_radius_v6b+(radial_index_v6b+0.5)*radial_step_px_v6b]

#Step 5-5：polar coordinate → Cartesian coordinate
#polar coordinate轉成以(0,0)為圓心的Cartesian coordinate
snp_v6_del_y[,x_cont_v6b:=
               snp_radius_v6b*cos(snp_theta_v6b)]

snp_v6_del_y[,y_cont_v6b:=
               snp_radius_v6b*sin(snp_theta_v6b)]

#影像外圍留白
circle_margin_px_v6b<-20L

#從中心到圖片邊界需要的距離
circle_canvas_radius_v6b<-
  ceiling(outer_radius_v6b)+circle_margin_px_v6b

#完整正方形影像邊長
circle_image_side_v6b<-2L*circle_canvas_radius_v6b+1L

#Step 5-7：建立以 image 左上角為基準的「暫定 pixel coordinate」

#圓心在正方形影像中央
circle_center_v6b<-circle_canvas_radius_v6b+1L

#將連續x、y座標轉成暫定pixel座標
snp_v6_del_y[,pixel_col_temp_v6b:=round(circle_center_v6b+x_cont_v6b)]

#matrix row方向向下增加，因此y需要反轉
snp_v6_del_y[,pixel_row_temp_v6b:=round(circle_center_v6b-y_cont_v6b)]

#檢查有多少組SNP經過round後落到相同pixel
pixel_collision_v6b<-snp_v6_del_y[
  ,
  .N,
  by=.(pixel_row_temp_v6b,pixel_col_temp_v6b)
][N>1]

nrow(pixel_collision_v6b)

pixel_collision_v6b[
  ,
  .(
    duplicated_pixel_n=.N,
    affected_SNP_n=sum(N)
  )
]


inner_radius_v6b
outer_radius_v6b
circle_image_side_v6b
nrow(pixel_collision_v6b)


#Version 6B Step 6：尋找可避免pixel collision的最小放大倍率----------------

#建立一個函數：給定放大倍率後，重新計算pixel座標並統計collision數量
check_collision_v6b<-function(scale_factor){
  
  #放大後每個radial slot的厚度
  radial_step_test<-scale_factor
  
  #放大後希望最內圈每個angular slot的弧長
  target_arc_test<-scale_factor
  
  #重新計算inner radius
  inner_radius_test<-ceiling(
    target_arc_test/angle_per_snp_slot_rad_v6b
  )
  
  #重新計算outer radius
  outer_radius_test<-inner_radius_test+
    radial_slots_v6b*radial_step_test
  
  #重新計算canvas半徑與影像邊長
  circle_canvas_radius_test<-ceiling(outer_radius_test)+circle_margin_px_v6b
  circle_image_side_test<-2L*circle_canvas_radius_test+1L
  circle_center_test<-circle_canvas_radius_test+1L
  
  #重新計算每個SNP的radius
  snp_radius_test<-inner_radius_test+
    (snp_v6_del_y$radial_index_v6b+0.5)*radial_step_test
  
  #polar轉Cartesian
  x_cont_test<-snp_radius_test*cos(snp_v6_del_y$snp_theta_v6b)
  y_cont_test<-snp_radius_test*sin(snp_v6_del_y$snp_theta_v6b)
  
  #轉成pixel座標
  pixel_col_test<-round(circle_center_test+x_cont_test)
  pixel_row_test<-round(circle_center_test-y_cont_test)
  
  #統計collision
  collision_n<-data.table(
    pixel_row_test=pixel_row_test,
    pixel_col_test=pixel_col_test
  )[
    ,
    .N,
    by=.(pixel_row_test,pixel_col_test)
  ][
    N>1,
    .N
  ]
  
  #若沒有collision，回傳0
  if(length(collision_n)==0) collision_n<-0L
  
  #整理結果
  data.table(
    scale_factor=scale_factor,
    inner_radius_v6b=inner_radius_test,
    outer_radius_v6b=outer_radius_test,
    circle_image_side_v6b=circle_image_side_test,
    collision_n=collision_n
  )
}

#先測試一組倍率範圍
scale_test_result_v6b<-rbindlist(
  lapply(c(1,2,3,4,5,6,8,10),check_collision_v6b)
)

scale_test_result_v6b


#Version 6B Step 6：固定最終解析度與pixel座標----------------

#根據collision測試結果，選用最小且collision=0的倍率
scale_factor_v6b<-2

#每個radial slot的厚度
radial_step_px_v6b<-scale_factor_v6b

#最內圈希望每個angular slot至少有2 pixel弧長
target_arc_px_v6b<-scale_factor_v6b

#重新計算最終inner radius
inner_radius_v6b<-ceiling(target_arc_px_v6b/angle_per_snp_slot_rad_v6b)

#重新計算最終outer radius
outer_radius_v6b<-inner_radius_v6b+radial_slots_v6b*radial_step_px_v6b

#外圍留白維持20 pixel
circle_margin_px_v6b<-20L

#整體canvas半徑與正方形影像邊長
circle_canvas_radius_v6b<-ceiling(outer_radius_v6b)+circle_margin_px_v6b
circle_image_side_v6b<-2L*circle_canvas_radius_v6b+1L

#圓心位置
circle_center_v6b<-circle_canvas_radius_v6b+1L

#每個SNP的最終radius
snp_v6_del_y[,snp_radius_v6b:=inner_radius_v6b+(radial_index_v6b+0.5)*radial_step_px_v6b]

#polar轉Cartesian連續座標
snp_v6_del_y[,x_cont_v6b:=snp_radius_v6b*cos(snp_theta_v6b)]
snp_v6_del_y[,y_cont_v6b:=snp_radius_v6b*sin(snp_theta_v6b)]

#轉成最終pixel座標
snp_v6_del_y[,pixel_col_v6b:=round(circle_center_v6b+x_cont_v6b)]
snp_v6_del_y[,pixel_row_v6b:=round(circle_center_v6b-y_cont_v6b)]

#查看最終設定
scale_factor_v6b
inner_radius_v6b
outer_radius_v6b
circle_image_side_v6b

#檢查最終pixel座標是否仍有重疊
pixel_collision_final_v6b<-snp_v6_del_y[
  ,
  .N,
  by=.(pixel_row_v6b,pixel_col_v6b)
][N>1]

nrow(pixel_collision_final_v6b)

#確認row與col都落在影像範圍內
range(snp_v6_del_y$pixel_row_v6b)
range(snp_v6_del_y$pixel_col_v6b)

all(snp_v6_del_y$pixel_row_v6b>=1L & snp_v6_del_y$pixel_row_v6b<=circle_image_side_v6b)
all(snp_v6_del_y$pixel_col_v6b>=1L & snp_v6_del_y$pixel_col_v6b<=circle_image_side_v6b)


#Version 6B Step 7：先畫第一位受試者circular atlas----------------

#選第一位受試者
view_sample_id_v6b<-sample_ids_v6A[1]
view_code_col_v6b<-genotype_code_cols_v6A[1]
#Step 7-2：確認 genotype code
#確認仍然只有1–7，0保留給padding
sort(unique(snp_v6_del_y[[view_code_col_v6b]]))

#確認沒有NA
sum(is.na(snp_v6_del_y[[view_code_col_v6b]]))

#Step 7-3：先直接用 SNP 座標畫圓形圖
#建立每個SNP的顏色
snp_color_v6b<-genotype_color_v6A[snp_v6_del_y[[view_code_col_v6b]]+1L ]

#輸出圖片
circular_png_v6b<-paste0(view_sample_id_v6b,"Version6B_circular_atlas.png")

png(
  circular_png_v6b,
  width=5000,
  height=5000,
  res=300
)

par(
  mar=c(1,1,5,1),
  xpd=NA
)

#建立圓形畫布
plot(
  NA,
  xlim=c(-outer_radius_v6b-100,
         outer_radius_v6b+100),
  ylim=c(-outer_radius_v6b-100,
         outer_radius_v6b+100),
  asp=1,
  axes=FALSE,
  xlab="",
  ylab="",
  main=""
)

#畫所有SNP
points(
  snp_v6_del_y$x_cont_v6b,
  snp_v6_del_y$y_cont_v6b,
  pch=15,
  cex=0.12,
  col=snp_color_v6b
)

#Step 7-4：畫 chromosome gaps
#取得每條染色體sector結束角度
chr_boundary_clock_v6b<-
  chr_sector_v6b$sector_end_clock_deg_v6b

chr_boundary_theta_v6b<-
  pi/2-chr_boundary_clock_v6b*pi/180

#畫染色體分隔線
for(theta_now in chr_boundary_theta_v6b){
  
  segments(
    x0=inner_radius_v6b*cos(theta_now),
    y0=inner_radius_v6b*sin(theta_now),
    x1=outer_radius_v6b*cos(theta_now),
    y1=outer_radius_v6b*sin(theta_now),
    col="white",
    lwd=1
  )
}

#Step 7-5：加入 chromosome labels
#計算每條染色體sector中間角度
chr_sector_v6b[,sector_mid_clock_deg_v6b:=
                 (sector_start_clock_deg_v6b+
                    sector_end_clock_deg_v6b)/2]

chr_sector_v6b[,sector_mid_theta_v6b:=
                 pi/2-sector_mid_clock_deg_v6b*pi/180]

#label放在外圈再往外一點
chr_label_radius_v6b<-outer_radius_v6b+120

chr_sector_v6b[,`:=`(
  label_x_v6b=
    chr_label_radius_v6b*cos(sector_mid_theta_v6b),
  label_y_v6b=
    chr_label_radius_v6b*sin(sector_mid_theta_v6b)
)]

#標示Chr名稱
text(
  chr_sector_v6b$label_x_v6b,
  chr_sector_v6b$label_y_v6b,
  labels=paste0("Chr",chr_sector_v6b$Chr_id),
  cex=0.55
)

#Step 7-6：加入橫向 genotype legend
legend(
  x=0,
  y=outer_radius_v6b+430,
  legend=genotype_label_v6A,
  fill=genotype_color_v6A,
  border=NA,
  horiz=TRUE,
  bty="n",
  xjust=0.5,
  cex=0.6
)

title(
  main=paste0(
    view_sample_id_v6b,
    " - Version 6B Circular Genotype Atlas"
  ),
  line=2.5,
  cex.main=1
)

dev.off()

system(paste("open",shQuote(circular_png_v6b)))


#正式版 Version 6B Step 1：建立所有 1 Mb Hilbert block 的全域排列順序，並在染色體之間插入 separator block----
#1個 1 Mb window= 1個 64 × 64 Hilbert patch

#Step 1-1：取得每條染色體需要多少個 physical windows
#確認每條染色體從W1保留到最後一個physical window
chr_window_span_v6b<-snp_v6_del_y[
  ,.(max_window_v6b=max(window_id_v6A)),
  by=.(chr_order_v6A,Chr_id)
]

setorder(chr_window_span_v6b,chr_order_v6A)

chr_window_span_v6b



#Step 1-2：建立完整的 1 Mb window block
#每條染色體建立W1到最後一個window
window_blocks_v6b<-rbindlist(
  lapply(seq_len(nrow(chr_window_span_v6b)),function(i){
    data.table(
      chr_order_v6b=chr_window_span_v6b$chr_order_v6A[i],
      Chr_id=chr_window_span_v6b$Chr_id[i],
      window_id_v6b=seq_len(chr_window_span_v6b$max_window_v6b[i])
    )
  })
)

#建立window名稱
window_blocks_v6b[,window_key_v6b:=
                    paste0("Chr",Chr_id,"_W",window_id_v6b)]

#標記這是真實physical window block
window_blocks_v6b[,block_type_v6b:="window"]

#Step 1-3：加入每個 window 實際 SNP 數
#整理原始資料中每個window真正有多少SNP
window_snp_n_v6b<-snp_v6_del_y[
  ,.(SNP_n=.N),
  by=.(Chr_id,window_id_v6A)
]

#合併進block表
window_blocks_v6b[
  window_snp_n_v6b,
  SNP_n:=i.SNP_n,
  on=.(Chr_id,window_id_v6b=window_id_v6A)
]

#沒有SNP的physical window設為0
window_blocks_v6b[is.na(SNP_n),SNP_n:=0L]

#Step 1-4：在 chromosome 之間插入 separator block
#建立最終block順序：每條染色體後加入1個separator，最後ChrMT除外
block_layout_v6b<-rbindlist(
  lapply(seq_len(nrow(chr_window_span_v6b)),function(i){
    
    chr_now<-chr_window_span_v6b$Chr_id[i]
    chr_order_now<-chr_window_span_v6b$chr_order_v6A[i]
    
    #抓出這條染色體所有window blocks
    chr_blocks<-copy(
      window_blocks_v6b[Chr_id==chr_now]
    )
    
    #不是最後一條染色體就加separator
    if(i<nrow(chr_window_span_v6b)){
      separator_block<-data.table(
        chr_order_v6b=chr_order_now,
        Chr_id=NA_character_,
        window_id_v6b=NA_integer_,
        window_key_v6b=paste0("SEP_after_Chr",chr_now),
        block_type_v6b="separator",
        SNP_n=0L
      )
      
      chr_blocks<-rbind(
        chr_blocks,
        separator_block,
        fill=TRUE
      )
    }
    
    chr_blocks
  })
) 

#Step 1-5：建立 Version 6B 全域 block order
#依實際排列順序建立global block index
block_layout_v6b[,global_block_order_v6b:=seq_len(.N)]

#64×64 Hilbert block尺寸
block_side_v6b<-hilbert_side_v6A

block_layout_v6b[1:30]

#Version 6b Step 2：建立同心圓block layout----------------

#每個Hilbert block為64×64
block_side_v6b<-hilbert_side_v6A

#block之間額外保留8 pixel空隙
block_gap_px_v6b<-8L

#為避免固定方向的正方形在斜角位置互相重疊
#中心最小距離使用block對角線+gap
block_pitch_v6b<-ceiling(
  sqrt(2)*block_side_v6b
)+block_gap_px_v6b

block_side_v6b
block_pitch_v6b

#不同ring之間也使用相同安全距離
ring_spacing_v6b<-block_pitch_v6b

#最內圈半徑先設定為一個block pitch
inner_ring_radius_v6b<-block_pitch_v6b


#越外面的 ring 半徑越大，所以能放的 block 也越多。
#給定ring radius，計算該ring最多可安全放多少個block
calc_ring_capacity_v6b<-function(radius){
  
  ratio<-block_pitch_v6b/(2*radius)
  
  if(ratio>=1){
    return(1L)
  }
  
  as.integer(
    floor(pi/asin(ratio))
  )
}


#包含所有1 Mb windows與23個chromosome separators
total_block_n_v6b<-nrow(block_layout_v6b)

total_block_n_v6b

#從內圈開始建立ring，直到總容量足以容納所有block
ring_list_v6b<-list()

current_radius_v6b<-inner_ring_radius_v6b
current_capacity_v6b<-0L
ring_id_v6b<-1L

while(current_capacity_v6b<total_block_n_v6b){
  
  capacity_now<-calc_ring_capacity_v6b(
    current_radius_v6b
  )
  
  ring_list_v6b[[ring_id_v6b]]<-data.table(
    ring_id_inner_to_outer_v6b=ring_id_v6b,
    ring_radius_v6b=current_radius_v6b,
    ring_capacity_v6b=capacity_now
  )
  
  current_capacity_v6b<-
    current_capacity_v6b+capacity_now
  
  current_radius_v6b<-
    current_radius_v6b+ring_spacing_v6b
  
  ring_id_v6b<-ring_id_v6b+1L
}

ring_layout_v6b<-rbindlist(ring_list_v6b)

#正式排列由外圈往內圈
setorder(ring_layout_v6b,-ring_radius_v6b)

#外圈=ring 1
ring_layout_v6b[,ring_order_v6b:=seq_len(.N)]


ring_layout_v6b

#累積每個ring容量
ring_layout_v6b[,cumulative_capacity_v6b:=
                  cumsum(ring_capacity_v6b)]

#每個ring負責的第一個global block
ring_layout_v6b[,block_start_order_v6b:=
                  shift(cumulative_capacity_v6b,fill=0L)+1L]

#每個ring負責的最後一個global block
ring_layout_v6b[,block_end_order_v6b:=
                  pmin(cumulative_capacity_v6b,total_block_n_v6b)]

#實際有使用的block數量
ring_layout_v6b[,used_block_n_v6b:=
                  pmax(
                    0L,
                    block_end_order_v6b-block_start_order_v6b+1L
                  )]



#2-7 建立每個 block 的圓形中心座標
#將global block分配到各個ring
block_position_map_v6b<-rbindlist(
  lapply(seq_len(nrow(ring_layout_v6b)),function(i){
    
    n_use<-ring_layout_v6b$used_block_n_v6b[i]
    
    if(n_use<=0){
      return(NULL)
    }
    
    ring_capacity_now<-
      ring_layout_v6b$ring_capacity_v6b[i]
    
    position_now<-seq_len(n_use)
    
    #12點鐘=0度，順時針增加
    clock_deg_now<-
      (position_now-1L)*
      360/ring_capacity_now
    
    #轉成R標準theta
    theta_now<-
      pi/2-clock_deg_now*pi/180
    
    radius_now<-
      ring_layout_v6b$ring_radius_v6b[i]
    
    data.table(
      global_block_order_v6b=
        ring_layout_v6b$block_start_order_v6b[i]:
        ring_layout_v6b$block_end_order_v6b[i],
      
      ring_order_v6b=
        ring_layout_v6b$ring_order_v6b[i],
      
      ring_radius_v6b=radius_now,
      
      ring_capacity_v6b=ring_capacity_now,
      
      position_in_ring_v6b=position_now,
      
      block_clock_deg_v6b=clock_deg_now,
      
      block_theta_v6b=theta_now,
      
      block_center_x_v6b=
        radius_now*cos(theta_now),
      
      block_center_y_v6b=
        radius_now*sin(theta_now)
    )
  })
)


#將圓形位置加入每個1 Mb window / separator block
block_layout_v6b[
  block_position_map_v6b,
  `:=`(
    ring_order_v6b=i.ring_order_v6b,
    ring_radius_v6b=i.ring_radius_v6b,
    ring_capacity_v6b=i.ring_capacity_v6b,
    position_in_ring_v6b=i.position_in_ring_v6b,
    block_clock_deg_v6b=i.block_clock_deg_v6b,
    block_theta_v6b=i.block_theta_v6b,
    block_center_x_v6b=i.block_center_x_v6b,
    block_center_y_v6b=i.block_center_y_v6b
  ),
  on="global_block_order_v6b"
]


#Step 3-1：建立每個 block 的四個邊界

#Version 6b Step 3：建立空的64×64 Hilbert block layout----------------

#每個block的一半邊長
block_half_v6b<-block_side_v6b/2

#計算每個固定方向正方形的邊界
block_layout_v6b[,`:=`(
  block_xmin_v6b=block_center_x_v6b-block_half_v6b,
  block_xmax_v6b=block_center_x_v6b+block_half_v6b,
  block_ymin_v6b=block_center_y_v6b-block_half_v6b,
  block_ymax_v6b=block_center_y_v6b+block_half_v6b
)]



#建立繪圖類別
block_layout_v6b[,plot_type_v6b:=fcase(
  block_type_v6b=="separator","Chromosome separator",
  block_type_v6b=="window"&SNP_n==0L,"Empty window",
  block_type_v6b=="window"&SNP_n>0L,"Window with SNP"
)]

table(block_layout_v6b$plot_type_v6b)


#QC圖使用的顏色
block_fill_color_v6b<-c(
  "Window with SNP"="grey85",
  "Empty window"="white",
  "Chromosome separator"="grey25"
)


#找每條染色體第一個window的位置
chr_label_map_v6b<-block_layout_v6b[
  block_type_v6b=="window",
  .SD[1],
  by=chr_order_v6b
]

#標籤往該block所在半徑方向再向外移
chr_label_radius_v6b<-
  max(block_layout_v6b$ring_radius_v6b)+
  block_side_v6b*1.8

chr_label_map_v6b[,label_theta_v6b:=
                    atan2(block_center_y_v6b,block_center_x_v6b)]

chr_label_map_v6b[,`:=`(
  label_x_v6b=
    chr_label_radius_v6b*cos(label_theta_v6b),
  label_y_v6b=
    chr_label_radius_v6b*sin(label_theta_v6b)
)]



#所有block最外側位置
plot_limit_v6b<-max(
  abs(c(
    block_layout_v6b$block_xmin_v6b,
    block_layout_v6b$block_xmax_v6b,
    block_layout_v6b$block_ymin_v6b,
    block_layout_v6b$block_ymax_v6b
  ))
)

#再多留一些空間給chromosome label
plot_limit_v6b<-plot_limit_v6b+
  block_side_v6b*3

#輸出空block layout QC圖
block_layout_png_v6b<-
  "Version6b_empty_Hilbert_block_layout.png"

png(
  filename=block_layout_png_v6b,
  width=3000,
  height=3000,
  res=200
)

par(
  mar=c(1,1,4,1),
  xpd=NA
)

#建立正方形畫布
plot(
  NA,
  xlim=c(-plot_limit_v6b,plot_limit_v6b),
  ylim=c(-plot_limit_v6b,plot_limit_v6b),
  asp=1,
  axes=FALSE,
  xlab="",
  ylab="",
  main="Version 6b - Concentric 1-Mb Hilbert Block Layout"
)

#逐個畫64×64固定方向正方形
for(i in seq_len(nrow(block_layout_v6b))){
  
  rect(
    xleft=block_layout_v6b$block_xmin_v6b[i],
    ybottom=block_layout_v6b$block_ymin_v6b[i],
    xright=block_layout_v6b$block_xmax_v6b[i],
    ytop=block_layout_v6b$block_ymax_v6b[i],
    
    col=block_fill_color_v6b[
      block_layout_v6b$plot_type_v6b[i]
    ],
    
    border="grey45",
    lwd=0.35
  )
}

#加染色體名稱
text(
  x=chr_label_map_v6b$label_x_v6b,
  y=chr_label_map_v6b$label_y_v6b,
  labels=paste0("Chr",chr_label_map_v6b$Chr_id),
  cex=0.55
)

#加入圖例
legend(
  "topright",
  legend=c(
    "Window with SNP",
    "Empty window",
    "Chromosome separator"
  ),
  fill=block_fill_color_v6b[
    c(
      "Window with SNP",
      "Empty window",
      "Chromosome separator"
    )
  ],
  border="grey45",
  bty="n",
  cex=0.7
)

dev.off()

#Mac直接開啟
system(paste("open",shQuote(block_layout_png_v6b)))

#Version 6b Step 4：建立正式circular Hilbert genotype atlas----------------

#每個64×64 block外圍再保留一些canvas空間
canvas_margin_v6b<-30L

#根據最外圈block決定整張正方形canvas大小
max_block_radius_v6b<-max(
  sqrt(
    block_layout_v6b$block_center_x_v6b^2+
      block_layout_v6b$block_center_y_v6b^2
  )
)

canvas_radius_v6b<-ceiling(
  max_block_radius_v6b+
    block_side_v6b/2+
    canvas_margin_v6b
)

canvas_side_v6b<-2L*canvas_radius_v6b+1L
canvas_center_v6b<-canvas_radius_v6b+1L

canvas_side_v6b

#每個block左上角pixel位置
block_layout_v6b[,block_col_start_v6b:=
                   round(
                     canvas_center_v6b+
                       block_center_x_v6b-
                       (block_side_v6b-1)/2
                   )
]

block_layout_v6b[,block_row_start_v6b:=
                   round(
                     canvas_center_v6b-
                       block_center_y_v6b-
                       (block_side_v6b-1)/2
                   )
]

#每個block右下角pixel位置
block_layout_v6b[,block_col_end_v6b:=
                   block_col_start_v6b+block_side_v6b-1L]

block_layout_v6b[,block_row_end_v6b:=
                   block_row_start_v6b+block_side_v6b-1L]


#只取真正的1 Mb window blocks，不包含separator
window_position_map_v6b<-block_layout_v6b[
  block_type_v6b=="window",
  .(
    Chr_id,
    window_id_v6b,
    global_block_order_v6b,
    ring_order_v6b,
    position_in_ring_v6b,
    block_row_start_v6b,
    block_col_start_v6b
  )
]

#把block位置合併回每個SNP
snp_v6_del_y[
  window_position_map_v6b,
  `:=`(
    global_block_order_v6b=i.global_block_order_v6b,
    ring_order_v6b=i.ring_order_v6b,
    position_in_ring_v6b=i.position_in_ring_v6b,
    block_row_start_v6b=i.block_row_start_v6b,
    block_col_start_v6b=i.block_col_start_v6b
  ),
  on=.(Chr_id,window_id_v6A=window_id_v6b)
]

#Version6b最終SNP pixel座標
snp_v6_del_y[,atlas_row_v6b:=
               block_row_start_v6b+patch_row_v6A-1L]

snp_v6_del_y[,atlas_col_v6b:=
               block_col_start_v6b+patch_col_v6A-1L]

#第一位受試者
view_sample_id_v6b<-sample_ids_v6A[1]
view_code_col_v6b<-genotype_code_cols_v6A[1]

#建立全padding canvas
atlas_matrix_test_v6b<-matrix(
  as.raw(0L),
  nrow=canvas_side_v6b,
  ncol=canvas_side_v6b
)

#將618,563個真實SNP填到Version6b新座標
atlas_matrix_test_v6b[
  cbind(
    snp_v6_del_y$atlas_row_v6b,
    snp_v6_del_y$atlas_col_v6b
  )
]<-as.raw(
  snp_v6_del_y[[view_code_col_v6b]]
)
real_snp_pixel_n_v6b<-sum(
  atlas_matrix_test_v6b!=as.raw(0L)
)

real_snp_pixel_n_v6b
nrow(snp_v6_del_y)

real_snp_pixel_n_v6b==nrow(snp_v6_del_y)


code_readback_v6b<-as.integer(
  atlas_matrix_test_v6b[
    cbind(
      snp_v6_del_y$atlas_row_v6b,
      snp_v6_del_y$atlas_col_v6b
    )
  ]
)

original_code_v6b<-
  snp_v6_del_y[[view_code_col_v6b]]

all(code_readback_v6b==original_code_v6b)


#Chr6在Version6a atlas中的row
chr6_row_start_v6A<-(6L-1L)*hilbert_side_v6A+1L
chr6_row_end_v6A<-6L*hilbert_side_v6A

#Window32在Version6a atlas中的column
w32_col_start_v6A<-(32L-1L)*hilbert_side_v6A+1L
w32_col_end_v6A<-32L*hilbert_side_v6A

#Version6a的Chr6_W32
patch_chr6_w32_v6A<-
  atlas_matrix_all_v6A[[view_sample_id_v6b]][
    chr6_row_start_v6A:chr6_row_end_v6A,
    w32_col_start_v6A:w32_col_end_v6A
  ]


chr6_w32_block_v6b<-block_layout_v6b[
  Chr_id=="6"&window_id_v6b==32
]

patch_chr6_w32_v6b<-atlas_matrix_test_v6b[
  chr6_w32_block_v6b$block_row_start_v6b:
    chr6_w32_block_v6b$block_row_end_v6b,
  
  chr6_w32_block_v6b$block_col_start_v6b:
    chr6_w32_block_v6b$block_col_end_v6b
]

identical(
  patch_chr6_w32_v6A,
  patch_chr6_w32_v6b
)

 
#Version 6b Step 5：畫第一位受試者正式circular Hilbert atlas----------------

#確認目前查看的sample
view_sample_id_v6b<-sample_ids_v6A[1]
view_sample_id_v6b

#設定輸出檔名
atlas_png_file_v6b<-paste0(
  view_sample_id_v6b,
  "_Version6b_concentric_Hilbert_atlas.png"
)

#圖片大小
png(
  filename=atlas_png_file_v6b,
  width=3500,
  height=3500,
  res=250
)

#上方留空間放橫向legend與title
par(
  mar=c(1,1,5,1),
  xpd=NA
)

#將raw matrix轉為integer
atlas_plot_v6b<-matrix(
  as.integer(atlas_matrix_test_v6b),
  nrow=canvas_side_v6b,
  ncol=canvas_side_v6b
)

#畫完整Version6b atlas
image(
  x=seq_len(canvas_side_v6b),
  y=seq_len(canvas_side_v6b),
  z=t(atlas_plot_v6b[nrow(atlas_plot_v6b):1L,,drop=FALSE]),
  col=genotype_color_v6A,
  breaks=seq(-0.5,7.5,by=1),
  axes=FALSE,
  xlab="",
  ylab="",
  xaxs="i",
  yaxs="i",
  asp=1,
  useRaster=TRUE
)

#抓出23個chromosome separator blocks
separator_blocks_v6b<-block_layout_v6b[
  block_type_v6b=="separator"
]

#在圖上標示separator block
for(i in seq_len(nrow(separator_blocks_v6b))){
  
  #image的y方向與matrix row方向相反，所以重新轉換
  xleft<-separator_blocks_v6b$block_col_start_v6b[i]-0.5
  xright<-separator_blocks_v6b$block_col_end_v6b[i]+0.5
  
  ybottom<-canvas_side_v6b-
    separator_blocks_v6b$block_row_end_v6b[i]+0.5
  
  ytop<-canvas_side_v6b-
    separator_blocks_v6b$block_row_start_v6b[i]+1.5
  
  rect(
    xleft=xleft,
    ybottom=ybottom,
    xright=xright,
    ytop=ytop,
    border="black",
    lwd=1.2
  )
}

#抓每條染色體第一個window
chr_first_block_v6b<-block_layout_v6b[
  block_type_v6b=="window",
  .SD[1],
  by=chr_order_v6b
]

#計算每個label相對於整張圖中心的方向
chr_first_block_v6b[,dx_v6b:=
                      block_center_x_v6b]

chr_first_block_v6b[,dy_v6b:=
                      block_center_y_v6b]

chr_first_block_v6b[,label_distance_v6b:=
                      sqrt(dx_v6b^2+dy_v6b^2)+90]

chr_first_block_v6b[,label_theta_v6b:=
                      atan2(dy_v6b,dx_v6b)]

#換成image座標
chr_first_block_v6b[,label_x_v6b:=
                      canvas_center_v6b+
                      label_distance_v6b*cos(label_theta_v6b)]

chr_first_block_v6b[,label_y_v6b:=
                      canvas_center_v6b+
                      label_distance_v6b*sin(label_theta_v6b)]

#畫Chr1–22、X、MT
text(
  x=chr_first_block_v6b$label_x_v6b,
  y=chr_first_block_v6b$label_y_v6b,
  labels=paste0("Chr",chr_first_block_v6b$Chr_id),
  cex=0.45,
  font=2
)

#圖例放在整個圓形atlas上方
legend(
  x=canvas_side_v6b/2,
  y=canvas_side_v6b+canvas_side_v6b*0.035,
  legend=genotype_label_v6A,
  fill=genotype_color_v6A,
  border=NA,
  horiz=TRUE,
  bty="n",
  xjust=0.5,
  yjust=0.5,
  cex=0.65,
  xpd=NA
)

#主標題
mtext(
  paste0(
    view_sample_id_v6b,
    " - Version 6b Concentric Hilbert Genotype Atlas"
  ),
  side=3,
  line=2.8,
  cex=1,
  font=2
)

dev.off()

#開啟圖片
system(paste("open",shQuote(atlas_png_file_v6b)))

















