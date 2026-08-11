#Step 1：建立每個 SNP 所屬的 1 Mb window 把線性結構整理好
# Step 1：依照 Start 把每個 SNP 分配到 1 Mb window ----

#設定每個 physical window 的大小為 1,000,000 bp
window_size_bp_v6 <- 1000000L

#複製一份資料，避免直接改到原始資料
snp_v6_chr <- copy(snp_v6_del_y)

#把 Start 轉成數值，後續才能做 window 切割
snp_v6_chr[, Start_num := as.numeric(Start)]

#建立染色體順序，後續畫圖時會依照 1~22、X、MT 排列
snp_v6_chr[, chr_order_v6 := match(Chr_id, c(as.character(1:22), "X", "MT"))]

#依照 Start 所在位置計算第幾個 1 Mb window
#例如 1~1,000,000 為第1個window；1,000,001~2,000,000 為第2個window
snp_v6_chr[, window_id_v6 := floor((Start_num - 1) / window_size_bp_v6) + 1L]

#建立每個 window 的實際 bp 起點
#例如第1個window起點=1；第2個window起點=1,000,001
snp_v6_chr[, window_bp_start_v6 := (window_id_v6 - 1L) * window_size_bp_v6 + 1L]

#建立每個 window 的實際 bp 終點
#例如第1個window終點=1,000,000；第2個window終點=2,000,000
snp_v6_chr[, window_bp_end_v6 := window_id_v6 * window_size_bp_v6]

#建立 window 名稱，之後做統計與轉圖會比較方便
#例如 Chr1_W1、Chr1_W2、ChrX_W3、ChrMT_W1
snp_v6_chr[, window_key_v6 := paste0("Chr", Chr_id, "_W", window_id_v6)]

#依照染色體順序、Start位置、probeset_id 排序
#確保後續的線性順序符合生物位置
setorder(snp_v6_chr, chr_order_v6, Start_num, probeset_id)

#先檢查前幾筆資料
head(snp_v6_chr[, .(
  Chr_id, Start, Start_num, chr_order_v6,
  window_id_v6, window_bp_start_v6, window_bp_end_v6, window_key_v6
)], 10)

#查看每條染色體前幾個 SNP 是否有被正確分到 window
snp_v6_chr[, .(
  first_start = min(Start_num, na.rm = TRUE),
  last_start = max(Start_num, na.rm = TRUE),
  n_window = uniqueN(window_id_v6)
), by = .(chr_order_v6, Chr_id)][order(chr_order_v6)]


#Step 2：統計每條染色體每個 window 有多少 SNP，並找出該染色體的最大 window 高度
#width  = 該染色體的1 Mb window數
#height = 該染色體中「SNP最多的window」的SNP數

#沿用Step 1資料，但後續統一使用小寫v6b
snp_v6b_chr<-copy(snp_v6_chr)

setnames(
  snp_v6b_chr,
  old=c(
    "chr_order_v6",
    "window_id_v6",
    "window_bp_start_v6",
    "window_bp_end_v6",
    "window_key_v6"
  ),
  new=c(
    "chr_order_v6b",
    "window_id_v6b",
    "window_bp_start_v6b",
    "window_bp_end_v6b",
    "window_key_v6b"
  )
)

#統計每個實際有SNP的1 Mb window
window_snp_count_v6b<-snp_v6b_chr[
  ,
  .(
    SNP_n_v6b=.N,
    observed_start_min_v6b=min(Start_num,na.rm=TRUE),
    observed_start_max_v6b=max(Start_num,na.rm=TRUE)
  ),
  by=.(
    chr_order_v6b,
    Chr_id,
    window_id_v6b,
    window_bp_start_v6b,
    window_bp_end_v6b,
    window_key_v6b
  )
]

setorder(
  window_snp_count_v6b,
  chr_order_v6b,
  window_id_v6b
)

head(window_snp_count_v6b)


#每條染色體最後一個SNP的位置與最後一個1 Mb window
chr_window_span_v6b<-snp_v6b_chr[
  ,
  .(
    max_Start_v6b=max(Start_num,na.rm=TRUE),
    max_window_id_v6b=max(window_id_v6b),
    total_SNP_v6b=.N
  ),
  by=.(chr_order_v6b,Chr_id)
]

setorder(chr_window_span_v6b,chr_order_v6b)

chr_window_span_v6b

#建立每條染色體從W1到最後一個window的完整window map
complete_window_map_v6b<-rbindlist(
  lapply(seq_len(nrow(chr_window_span_v6b)),function(i){
    
    chr_now<-chr_window_span_v6b$Chr_id[i]
    chr_order_now<-chr_window_span_v6b$chr_order_v6b[i]
    max_window_now<-chr_window_span_v6b$max_window_id_v6b[i]
    
    data.table(
      chr_order_v6b=chr_order_now,
      Chr_id=chr_now,
      window_id_v6b=seq_len(max_window_now)
    )
  })
)

#建立每個window的physical range
complete_window_map_v6b[,window_bp_start_v6b:=
                          (window_id_v6b-1L)*1000000L+1L]

complete_window_map_v6b[,window_bp_end_v6b:=
                          window_id_v6b*1000000L]

complete_window_map_v6b[,window_key_v6b:=
                          paste0("Chr",Chr_id,"_W",window_id_v6b)]


#把observed SNP count放回完整physical window map
complete_window_map_v6b[
  window_snp_count_v6b,
  `:=`(
    SNP_n_v6b=i.SNP_n_v6b,
    observed_start_min_v6b=i.observed_start_min_v6b,
    observed_start_max_v6b=i.observed_start_max_v6b
  ),
  on=.(Chr_id,window_id_v6b)
]

#完全沒有SNP的window設成0
complete_window_map_v6b[
  is.na(SNP_n_v6b),
  SNP_n_v6b:=0L
]

setorder(
  complete_window_map_v6b,
  chr_order_v6b,
  window_id_v6b
)


#每條染色體有多少window沒有任何SNP
empty_window_summary_v6b<-complete_window_map_v6b[
  ,
  .(
    total_window_n_v6b=.N,
    window_with_SNP_n_v6b=sum(SNP_n_v6b>0),
    empty_window_n_v6b=sum(SNP_n_v6b==0)
  ),
  by=.(chr_order_v6b,Chr_id)
]

empty_window_summary_v6b


#每條染色體的window SNP密度摘要
chr_matrix_size_v6b<-complete_window_map_v6b[
  ,
  .(
    window_n_v6b=.N,
    max_SNP_per_window_v6b=max(SNP_n_v6b),
    mean_SNP_per_window_v6b=as.numeric(mean(SNP_n_v6b)),
    median_SNP_per_window_v6b=as.numeric(median(SNP_n_v6b)),
    empty_window_n_v6b=sum(SNP_n_v6b==0L),
    total_SNP_v6b=sum(SNP_n_v6b)
  ),
  by=.(chr_order_v6b,Chr_id)
]

setorder(
  chr_matrix_size_v6b,
  chr_order_v6b
)

chr_matrix_size_v6b

#找每條染色體SNP最密集的1 Mb window
max_window_each_chr_v6b<-complete_window_map_v6b[
  ,
  .SD[which.max(SNP_n_v6b)],
  by=.(chr_order_v6b,Chr_id)
][
  ,
  .(
    chr_order_v6b,
    Chr_id,
    window_id_v6b,
    window_key_v6b,
    SNP_n_v6b,
    window_bp_start_v6b,
    window_bp_end_v6b
  )
]

max_window_each_chr_v6b

#所有window SNP加總
sum(complete_window_map_v6b$SNP_n_v6b)

#原始SNP總數
nrow(snp_v6b_chr)

#兩者必須完全相同
sum(complete_window_map_v6b$SNP_n_v6b)==
  nrow(snp_v6b_chr)

chr_matrix_size_v6b[
  Chr_id=="6"
]

complete_window_map_v6b[
  Chr_id=="6"&
    window_id_v6b>=28&
    window_id_v6b<=36,
  .(
    window_key_v6b,
    SNP_n_v6b,
    observed_start_min_v6b,
    observed_start_max_v6b
  )
]
chr_matrix_size_v6b[
  Chr_id=="MT"
]

complete_window_map_v6b[
  Chr_id=="MT"
]

#Step 3：正式建立單一染色體的 rectangular genotype matrix
#windows中用snake排
#Version 6b Step 3：建立單條染色體的64×56 window strip----------------

#每個1 Mb window固定64 columns × 56 rows
window_width_v6b<-64L
window_height_v6b<-56L
window_capacity_v6b<-window_width_v6b*window_height_v6b

window_capacity_v6b

#抓Chr6
chr_snp_v6b<-copy(
  snp_v6b_chr[Chr_id==test_chr_v6b]
)

#依1 Mb window、Start排序
#Start相同時用probeset_id固定順序
setorder(
  chr_snp_v6b,
  window_id_v6b,
  Start_num,
  probeset_id
)

#每個window內重新從1開始排序
chr_snp_v6b[,snp_order_within_window_v6b:=
              seq_len(.N),
            by=window_id_v6b
]
#每64個SNP換到下一row
chr_snp_v6b[,local_row_v6b:=
              ((snp_order_within_window_v6b-1L)%/%window_width_v6b)+1L
]
#每個row內的位置1–64
chr_snp_v6b[,position_in_row_v6b:=
              ((snp_order_within_window_v6b-1L)%%window_width_v6b)+1L
]

#奇數row左→右；偶數row右→左
chr_snp_v6b[,local_col_v6b:=
              fifelse(
                local_row_v6b%%2L==1L,
                position_in_row_v6b,
                window_width_v6b-position_in_row_v6b+1L
              )
]

#把local column轉成整條染色體的global column
chr_snp_v6b[,global_col_v6b:=
              (window_id_v6b-1L)*window_width_v6b+
              local_col_v6b
]

#row不用變
chr_snp_v6b[,global_row_v6b:=local_row_v6b]

#Chr6 Start-based window數
chr_window_n_v6b<-chr_window_span_v6b[
  Chr_id==test_chr_v6b,
  max_window_id_v6b
]

#整條Chr6 image的寬度
chr_strip_width_v6b<-
  chr_window_n_v6b*window_width_v6b

#高度固定56
chr_strip_height_v6b<-
  window_height_v6b

chr_window_n_v6b
chr_strip_width_v6b
chr_strip_height_v6b

#0代表padding
chr_strip_matrix_v6b<-matrix(
  0L,
  nrow=chr_strip_height_v6b,
  ncol=chr_strip_width_v6b
)

dim(chr_strip_matrix_v6b)

#抓第一位sample genotype code
chr_snp_v6b[,genotype_code_v6b:=
              as.integer(get(test_code_col_v6b))
]

#填進長條矩陣
chr_strip_matrix_v6b[
  cbind(
    chr_snp_v6b$global_row_v6b,
    chr_snp_v6b$global_col_v6b
  )
]<-chr_snp_v6b$genotype_code_v6b


sum(chr_strip_matrix_v6b!=0L)
nrow(chr_snp_v6b)

sum(chr_strip_matrix_v6b!=0L)==
  nrow(chr_snp_v6b)

#Chr6展開版QC圖片
chr_strip_file_v6b<-paste0(
  test_sample_id_v6b,
  "_Chr",
  test_chr_v6b,
  "_64x56_window_strip_v6b.png"
)

png(
  filename=chr_strip_file_v6b,
  width=12000,
  height=1200,
  res=300
)

par(
  mar=c(4,4,3,1)
)

image(
  x=seq_len(ncol(chr_strip_matrix_v6b)),
  y=seq_len(nrow(chr_strip_matrix_v6b)),
  z=t(
    chr_strip_matrix_v6b[
      nrow(chr_strip_matrix_v6b):1,
      ,
      drop=FALSE
    ]
  ),
  col=genotype_color_v6A,
  breaks=seq(-0.5,7.5,by=1),
  axes=FALSE,
  xlab="Genomic position",
  ylab="Snake row",
  main=paste0(
    test_sample_id_v6b,
    " - Chr",
    test_chr_v6b,
    " Genotype Strip"
  ),
  useRaster=TRUE,
  xaxs="i",
  yaxs="i"
)
#每64 pixels就是1 Mb
window_boundary_x_v6b<-
  seq(
    window_width_v6b+0.5,
    chr_strip_width_v6b,
    by=window_width_v6b
  )

abline(
  v=window_boundary_x_v6b,
  col="grey70",
  lwd=0.25
)

major_window_v6b<-seq(
  10,
  chr_window_n_v6b,
  by=10
)

major_x_v6b<-
  major_window_v6b*window_width_v6b+0.5

abline(
  v=major_x_v6b,
  col="black",
  lwd=0.8
)

axis(
  side=1,
  at=major_x_v6b,
  labels=paste0(
    major_window_v6b,
    " Mb"
  ),
  cex.axis=0.6
)

axis(
  side=2,
  at=c(1,14,28,42,56),
  labels=c(1,14,28,42,56),
  las=1,
  cex.axis=0.6
)

box()

dev.off()

system(
  paste(
    "open",
    shQuote(chr_strip_file_v6b)
  )
)
#Version 6b Step 4：將Chr6 rectangular strip捲成spiral----------------

#spiral相鄰兩圈之間額外保留12個cell寬度的空隙
spiral_gap_v6b<-12

#一整圈半徑增加量
#56 = strip本身厚度；+12 = 相鄰兩圈gap
spiral_pitch_v6b<-
  window_height_v6b+
  spiral_gap_v6b

#spiral中心線最開始的半徑
#確保最內側不會碰到圓心
spiral_start_radius_v6b<-
  window_height_v6b/2+25

#Archimedean spiral參數
#每轉完整一圈2pi，半徑增加spiral_pitch_v6b
spiral_b_v6b<-
  spiral_pitch_v6b/(2*pi)

spiral_pitch_v6b
spiral_start_radius_v6b
spiral_b_v6b

#整條Chr6總共有多少horizontal columns
spiral_width_v6b<-ncol(chr_strip_matrix_v6b)

#建立每個column edge的theta
#W個columns需要W+1個edge
spiral_theta_edge_v6b<-numeric(
  spiral_width_v6b+1L
)

spiral_theta_edge_v6b[1]<-0

#用近似arc-length=1的方式逐格往前
for(j in seq_len(spiral_width_v6b)){
  
  theta_now<-
    spiral_theta_edge_v6b[j]
  
  r_now<-
    spiral_start_radius_v6b+
    spiral_b_v6b*theta_now
  
  #第一次估計
  dtheta_now<-
    1/sqrt(r_now^2+spiral_b_v6b^2)
  
  #用midpoint再修正一次
  theta_mid<-
    theta_now+dtheta_now/2
  
  r_mid<-
    spiral_start_radius_v6b+
    spiral_b_v6b*theta_mid
  
  dtheta_now<-
    1/sqrt(r_mid^2+spiral_b_v6b^2)
  
  spiral_theta_edge_v6b[j+1L]<-
    theta_now+dtheta_now
}

#每個column edge的中心線半徑
spiral_center_radius_edge_v6b<-
  spiral_start_radius_v6b+
  spiral_b_v6b*spiral_theta_edge_v6b



#strip最外側edge
spiral_outer_radius_edge_v6b<-
  spiral_center_radius_edge_v6b+
  window_height_v6b/2

#strip最內側edge
spiral_inner_radius_edge_v6b<-
  spiral_center_radius_edge_v6b-
  window_height_v6b/2


#0 rad = 12點鐘，角度增加為順時針
spiral_outer_x_v6b<-
  spiral_outer_radius_edge_v6b*
  sin(spiral_theta_edge_v6b)

spiral_outer_y_v6b<-
  spiral_outer_radius_edge_v6b*
  cos(spiral_theta_edge_v6b)

spiral_inner_x_v6b<-
  spiral_inner_radius_edge_v6b*
  sin(spiral_theta_edge_v6b)

spiral_inner_y_v6b<-
  spiral_inner_radius_edge_v6b*
  cos(spiral_theta_edge_v6b)

#每個SNP所在column的左右theta edge
chr_snp_v6b[,spiral_theta_left_v6b:=
              spiral_theta_edge_v6b[global_col_v6b]]

chr_snp_v6b[,spiral_theta_right_v6b:=
              spiral_theta_edge_v6b[global_col_v6b+1L]]

#左右兩側spiral中心線半徑
chr_snp_v6b[,spiral_center_r_left_v6b:=
              spiral_start_radius_v6b+
              spiral_b_v6b*spiral_theta_left_v6b]

chr_snp_v6b[,spiral_center_r_right_v6b:=
              spiral_start_radius_v6b+
              spiral_b_v6b*spiral_theta_right_v6b]

#每一row的上、下邊界相對centerline的距離
chr_snp_v6b[,spiral_row_upper_offset_v6b:=
              window_height_v6b/2-
              (global_row_v6b-1L)]

chr_snp_v6b[,spiral_row_lower_offset_v6b:=
              window_height_v6b/2-
              global_row_v6b]

chr_snp_v6b[,spiral_x1_v6b:=
              (spiral_center_r_left_v6b+
                 spiral_row_upper_offset_v6b)*
              sin(spiral_theta_left_v6b)]

chr_snp_v6b[,spiral_y1_v6b:=
              (spiral_center_r_left_v6b+
                 spiral_row_upper_offset_v6b)*
              cos(spiral_theta_left_v6b)]
chr_snp_v6b[,spiral_x2_v6b:=
              (spiral_center_r_right_v6b+
                 spiral_row_upper_offset_v6b)*
              sin(spiral_theta_right_v6b)]

chr_snp_v6b[,spiral_y2_v6b:=
              (spiral_center_r_right_v6b+
                 spiral_row_upper_offset_v6b)*
              cos(spiral_theta_right_v6b)]
chr_snp_v6b[,spiral_x3_v6b:=
              (spiral_center_r_right_v6b+
                 spiral_row_lower_offset_v6b)*
              sin(spiral_theta_right_v6b)]

chr_snp_v6b[,spiral_y3_v6b:=
              (spiral_center_r_right_v6b+
                 spiral_row_lower_offset_v6b)*
              cos(spiral_theta_right_v6b)]
chr_snp_v6b[,spiral_x4_v6b:=
              (spiral_center_r_left_v6b+
                 spiral_row_lower_offset_v6b)*
              sin(spiral_theta_left_v6b)]

chr_snp_v6b[,spiral_y4_v6b:=
              (spiral_center_r_left_v6b+
                 spiral_row_lower_offset_v6b)*
              cos(spiral_theta_left_v6b)]


draw_chr_spiral_genotype_v6b<-function(
    dt,
    genotype_col,
    color_vec
){
  
  code_vec<-dt[[genotype_col]]
  
  for(code_now in 1:7){
    
    idx<-which(code_vec==code_now)
    
    if(length(idx)==0L){
      next
    }
    
    x_mat<-cbind(
      dt$spiral_x1_v6b[idx],
      dt$spiral_x2_v6b[idx],
      dt$spiral_x3_v6b[idx],
      dt$spiral_x4_v6b[idx],
      NA_real_
    )
    
    y_mat<-cbind(
      dt$spiral_y1_v6b[idx],
      dt$spiral_y2_v6b[idx],
      dt$spiral_y3_v6b[idx],
      dt$spiral_y4_v6b[idx],
      NA_real_
    )
    
    polygon(
      x=as.vector(t(x_mat)),
      y=as.vector(t(y_mat)),
      col=color_vec[code_now+1L],
      border=NA
    )
  }
}

#最外側最大半徑
spiral_max_radius_v6b<-
  max(spiral_outer_radius_edge_v6b)

#多留一些空間放Mb label
spiral_plot_limit_v6b<-
  spiral_max_radius_v6b+30

spiral_pdf_v6b<-paste0(
  test_sample_id_v6b,
  "_Chr",
  test_chr_v6b,
  "_spiral_genotype_v6b.pdf"
)

pdf(
  file=spiral_pdf_v6b,
  width=12,
  height=12,
  useDingbats=FALSE
)

par(
  mar=c(1,1,4,1),
  xpd=NA
)

plot(
  NA,
  xlim=c(-spiral_plot_limit_v6b,
         spiral_plot_limit_v6b),
  ylim=c(-spiral_plot_limit_v6b,
         spiral_plot_limit_v6b),
  asp=1,
  axes=FALSE,
  xlab="",
  ylab=""
)

polygon(
  x=c(
    spiral_outer_x_v6b,
    rev(spiral_inner_x_v6b)
  ),
  y=c(
    spiral_outer_y_v6b,
    rev(spiral_inner_y_v6b)
  ),
  col="white",
  border="grey75",
  lwd=0.5
)
draw_chr_spiral_genotype_v6b(
  dt=chr_snp_v6b,
  genotype_col="genotype_code_v6b",
  color_vec=genotype_color_v6A
)

#每一個1 Mb boundary的global column edge
window_edge_col_v6b<-seq(
  0L,
  chr_window_n_v6b*window_width_v6b,
  by=window_width_v6b
)

for(col_now in window_edge_col_v6b){
  
  edge_index<-col_now+1L
  
  theta_now<-
    spiral_theta_edge_v6b[edge_index]
  
  center_r_now<-
    spiral_start_radius_v6b+
    spiral_b_v6b*theta_now
  
  r_outer_now<-
    center_r_now+window_height_v6b/2
  
  r_inner_now<-
    center_r_now-window_height_v6b/2
  
  #目前是第幾Mb
  mb_now<-col_now/window_width_v6b
  
  #每10 Mb粗線，其餘1 Mb細線
  if(mb_now%%10==0){
    
    line_width_now<-1.2
    line_color_now<-"black"
    
  }else{
    
    line_width_now<-0.25
    line_color_now<-"grey65"
  }
  
  segments(
    x0=r_inner_now*sin(theta_now),
    y0=r_inner_now*cos(theta_now),
    x1=r_outer_now*sin(theta_now),
    y1=r_outer_now*cos(theta_now),
    col=line_color_now,
    lwd=line_width_now
  )
}

major_mb_v6b<-seq(
  10,
  chr_window_n_v6b,
  by=10
)

for(mb_now in major_mb_v6b){
  
  col_now<-
    mb_now*window_width_v6b
  
  theta_now<-
    spiral_theta_edge_v6b[col_now+1L]
  
  center_r_now<-
    spiral_start_radius_v6b+
    spiral_b_v6b*theta_now
  
  label_r_now<-
    center_r_now+
    window_height_v6b/2+
    5
  
  text(
    x=label_r_now*sin(theta_now),
    y=label_r_now*cos(theta_now),
    labels=paste0(mb_now," Mb"),
    cex=0.45
  )
}

text(
  x=0,
  y=0,
  labels=paste0("Chr",test_chr_v6b),
  cex=1.5,
  font=2
)

legend(
  "bottom",
  legend=genotype_label_v6A[-1],
  fill=genotype_color_v6A[-1],
  border=NA,
  horiz=TRUE,
  bty="n",
  cex=0.55,
  inset=-0.01
)

mtext(
  paste0(
    test_sample_id_v6b,
    " - Chr",
    test_chr_v6b,
    " Spiral Genotype Map"
  ),
  side=3,
  line=1.5,
  font=2
)

dev.off()

system(
  paste(
    "open",
    shQuote(spiral_pdf_v6b)
  )
)