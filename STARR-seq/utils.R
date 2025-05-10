# Utility files for motif.R script, created by Yuanzhong Pan (michapan93@gmail.com)
#
#
library(data.table)

# remove repeated peaks because of genome alt redundancy
peak_dedup <- function(peak_list){
  l <- length(peak_list)
  is_dup <- vector(length = l)
  for (i in 2:l){
    #print(paste("i",i,sep = " "))
    seq1 <- peak_list[i]$seq
    for (j in 2:i-1){
      #print(j)
      seq2 <- peak_list[i-j]$seq
      if (seq1 == seq2){
        is_dup[i] <- TRUE
        break
      }
    }
    if (i%%100 == 0){print(paste("Processed", i, sep = " "))}
  }
  return(peak_list[!is_dup])
}

# Count PWMs in a sequence
tf_mat <- function(peak_list, tflist, match.score = "90%"){
  rscale <- c(length(peak_list), length(tflist))
  cnt_mat <- matrix(nrow = rscale[1], ncol = rscale[2])
  one_count <- 0
  for (i in 1:rscale[1]){
    seq <- peak_list[i][[1]]
    for (j in 1:rscale[2]){
      one_count <- 0
      tf <- tflist[[j]]['motif']
      one_count <- countPWM(tf, seq, min.score = match.score)
      one_count <- one_count + countPWM(tf,reverseComplement(seq), min.score = match.score)
#      if (tflist[[j]]['extrainfo'] == "0"){
#        one_count <- one_count + countPWM(tf,complement(seq), min.score = match.score)
#        one_count <- one_count + countPWM(tf,reverseComplement(seq), min.score = match.score)
#      }
      cnt_mat[i,j] <-one_count      
    }
  }
  #for (tfname in 1:rscale[2]){colnames(cnt_mat)[tfname] <- tflist[[tfname]]['altname']}
  return(cnt_mat)
}


# Find empty cluster
count.empty.cluster <- function(count_mat){
  sum_vec <- rowSums(count_mat)
  num_empty <- sum(sum_vec == 0)
  return(num_empty)
}


# Calculate two-sided Fischer's exact test p value for co-occurrence 
co.occur <- function(count_mat, tflist){
  l <- nrow(count_mat)
  d <- length(tflist)
  name_vac <- list(length = d)
  occurance <- matrix(as.numeric(count_mat > 0), ncol = d)
  cooccur.p <- matrix(nrow = d, ncol = d)
  for (tf1 in 1:d){
    name_vac[tf1] <- tflist[[tf1]]['altname']
    for (tf2 in 1:d){
      if (sum(occurance[,tf1]) == 0 | sum(occurance[,tf2]) == 0 ){
        cooccur.p[tf1,tf2] <- 1
      }
      else {
        cooccur.p[tf1,tf2] <- fisher.test(occurance[,tf1],occurance[,tf2],alternative = "two.sided")$p
      }
    }
  }
  print("loops done")
  rownames(cooccur.p) <- name_vac
  colnames(cooccur.p) <- name_vac
  return(cooccur.p)
}

# For each CRE type, calculate a score that summarizes the pileup of peaks that belong to that type 
score_type <- function(tf_mat, type_mat, peak_mat){
  n_type <- dim(type_mat)[1]
  score_mat <- matrix(nrow = n_type, ncol = 1)
  n_tf <- dim(tf_mat)[1]
  for (i in 1:n_type){
    one_type <- 0
    for (j in 1:n_tf){
      tf_bin <- as.numeric(tf_mat[j,] > 0)
      if (identical(as.numeric(type_mat[i,]),tf_bin)){
        one_type <- one_type + peak_mat$pileup[j]
      }
    }
    score_mat[i] <- one_type
  }
  return(score_mat)
}

match_type <- function(tf_mat, type_mat){
  n_type <- dim(type_mat)[1]
  n_tf <- dim(tf_mat)[1]
  is_type <- vector(length = n_tf)
  for (i in 1:n_tf){
    tf_bin <- as.numeric(tf_mat[i,] > 0)
    for (j in 1:n_type){
      if (identical(as.numeric(type_mat[j,]),tf_bin)){
        is_type[i] <- TRUE
        break
      }
    }
  }
  return(is_type)
}

cooccur_vec_from_binary <- function(type_mat){
  l <- dim(type_mat)[1]
  w <- dim(type_mat)[2]
  n_vec <- (w**2 - w)/2
  cooccur_mat <- matrix(nrow = l, ncol = n_vec)
  for (i in 1:l){
    occur_vec <- type_mat[i,]
    cooccur_vec <- vector(length = n_vec)
    vec_pos <- 1
    for (a in 1:(w-1)){
      for (b in (a+1):w){
        if (occur_vec[a]==0){cooccur_vec[vec_pos] <- 0}
        else if (occur_vec[a] == occur_vec[b]){cooccur_vec[vec_pos] <- 1}
        vec_pos <- vec_pos +1
      }
    }
    cooccur_mat[i,] <- cooccur_vec
  }
  # Name the columns and rows of the new matrix
  name_vec <- vector(length = n_vec)
  vec_pos <- 1
  for (a in 1:(w-1)){
    for (b in (a+1):w){
      name_vec[vec_pos] <- paste(colnames(type_mat)[a],colnames(type_mat)[b],sep = "_")
      vec_pos <- vec_pos +1
    }
  }
  colnames(cooccur_mat)<-name_vec
  rownames(cooccur_mat)<-rownames(type_mat)
  #
  return(cooccur_mat)
}

cooccur_vec_from_mat <- function(co_mat){
  n <- dim(co_mat)[1]
  n_vec <- (n**2 - n)/2
  cooccur_vec <- vector(length = n_vec)
  name_vec <- vector(length = n_vec)
  vec_pos <- 1
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      cooccur_vec[vec_pos] <- co_mat[i,j]
      name_vec[vec_pos] <- paste(rownames(co_mat)[i],colnames(co_mat)[j],sep = "_")
      vec_pos <- vec_pos +1
    }
  }
  names(cooccur_vec) <- name_vec
  return(cooccur_vec)
}

annotate_type <- function(count_mat, type_mat){
  n_peaks <- dim(count_mat)[1]
  type_vec <- matrix(length = n_peaks)
  for (i in 1:n_peaks){
    annotated_types <- apply(df_peaks, 1, function(count_mat) {
    match_index <- which(apply(df_types, 1, function(type_row) identical(peak_row, type_row)))
    
    if (length(match_index) == 1) {
      return(rownames(df_types)[match_index])
    } else {
      return(NA)  # or "unmatched", depending on your use case
    }
  })
      if (identical(as.numeric(type_mat[i,]),tf_bin)){
        one_type <- one_type + peak_mat$pileup[j]
      }
    }
    score_mat[i] <- one_type
  }
  return(score_mat)
}


annotate_type <- function(count_mat, type_mat){
  # Convert to data.table
  dt_peaks <- as.data.table(count_mat)
  dt_peaks[, names(dt_peaks) := lapply(.SD, function(x) as.integer(x > 0))]
  dt_types <- as.data.table(type_mat)

  # Add a row identifier to types (row names as a column)
  dt_types[, type_id := rownames(type_mat)]
  
  # Create a unique string key for each row to speed up matching
  dt_peaks[, key := do.call(paste, c(.SD, sep = "_"))]
  dt_peaks[, original_order := .I]
  dt_types[, key := do.call(paste, c(.SD, sep = "_")), .SDcols = names(type_mat)]

  # Perform a fast join to annotate types
  dt_annotated <- merge(
    dt_peaks,
    dt_types[, .(key, type_id)],
    by = "key",
    all.x = TRUE,
    sort = FALSE
  )
  setorder(dt_annotated, original_order)
  return(dt_annotated$type_id)
}


welch_p <- function(m1, v1, n1, m2, v2, n2) {
  se <- sqrt(v1/n1 + v2/n2)
  t_stat <- (m1 - m2) / se
  df_num <- (v1/n1 + v2/n2)^2
  df_denom <- ((v1/n1)^2)/(n1 - 1) + ((v2/n2)^2)/(n2 - 1)
  df <- df_num / df_denom
  pval <- 2 * pt(-abs(t_stat), df)
  return(pval)
}
