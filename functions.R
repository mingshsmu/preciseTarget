find_crRNA <- function(input_seq,mutation_site,mutate_to=c("A","T/U","G","C"),
                       label="",Cas_protein=c("suCas12a2")){
  library(dplyr)
  library(stringr)
  # 1. check input ---------------------------------------------------------------
  mutation_site <- as.numeric(mutation_site)
  input_seq <- toupper(input_seq)
  mutate_to <- toupper(mutate_to)
  if(mutate_to == "T/U"){mutate_to <- "T"}
  input_seq_vector <- unlist(strsplit(input_seq,""))
  if(!mutate_to %in% c("A","T","U","G","C")){
    stop("Error! The mutate_to should be a character in 'A', 'T', 'U', 'G', 'C'. ")
  }
  Cas_proteins <- c("suCas12a2")
  if(!Cas_protein %in% Cas_proteins){
    stop("Error! Please input the right name of Cas_protein.")
  }
  # 1.1 check mutation_site
  input_length <- length(input_seq_vector)
  if(mutation_site < 25 | (input_length - mutation_site) < 25){
    stop("Input error! Please enter the correct sequence, and ensure that its two ends are 25 nucleotides away from the mutation site.")
  }
  # 1.2 check sequence
  ref <- setNames(c("A","T","U","G","C"),
                  c("A","T","U","G","C"))
  check_sequence <- ref[input_seq_vector]
  if(length(check_sequence) != length(input_seq_vector)){
    stop("Input error! The input sequence should only contains A/T/G/C (for DNA) or A/U/G/C (for RNA).")
  }
  # 2. detect is DNA or RNA ------------------------------------------------------
  T_contain_flag <- length(dplyr::contains(c("T"),vars = check_sequence))>0
  U_contain_flag <- length(dplyr::contains(c("U"),vars = check_sequence))>0
  if(all(T_contain_flag,U_contain_flag)){
    print("Input error!. The input sequence should only contains A/T/G/C (for DNA) or A/U/G/C (for RNA).")
  }
  
  input_type <- "DNA"
  if(U_contain_flag){
    input_type <- "RNA"
  }
  
  # 3. display the mutation site
  left_seq_index <- seq(mutation_site - 24,mutation_site -1,1)
  left_seq <- paste0(check_sequence[left_seq_index],collapse = "")
  right_seq_index <- seq(mutation_site + 1,mutation_site +24, 1)
  right_seq <- paste0(check_sequence[right_seq_index],collapse = "")
  mutate_base_from <- check_sequence[mutation_site]
  mutate_base_to <- mutate_to
  
  if(input_type=="RNA" & mutate_base_to == "T"){
    mutate_base_to <- "U"
  }
  
  input_sequence <- paste0(left_seq,mutate_base_from,right_seq)
  mutant_sequence <- paste0(left_seq,mutate_base_to,right_seq)
  seqs_wt_mut <- c(input_sequence,mutant_sequence)
  
  if(F){
    # for display
    left_seq_display <-  paste0("---",left_seq)
    right_seq_display <-  paste0(right_seq,"---")
    # 两端都带上'---'
    # 突变位点用蓝色（WT）和红色（Mut）展示，其余碱基用黑色展示
  }
  
  # 4. if is RNA, convert to DNA for calculation
  RNA2DNA_fun <- function(seq1){
    ref2RNA <- setNames(c("A","U","G","C","U"),
                        c("A","T","G","C","U"))
    seq1 <- ref2RNA[unlist(strsplit(seq1,""))]
    seq1 <- paste0(seq1,collapse = "")
    return(seq1)
  }
  
  # 5. searching candidates
  # functions
  find_crRNA_Cas12a2 <- function(){
    search_candidates_fun <- function(seq1){
      if(stringr::str_length(seq1)!=49){
        stop("Search candidates failed. The length of seq != 49.")
      }
      seq1 <- RNA2DNA_fun(seq1)
      candidates <- sapply(1:25,function(i){
        seq1 <- unlist(strsplit(seq1,""))
        candidate <- seq1[seq(i,(i+24),1)]
        candidate <- paste0(candidate,collapse = "")
        return(candidate)
      })
      return(candidates)
    }
    candidates_input <- search_candidates_fun(seq1 = input_sequence)
    candidates_mutant <- search_candidates_fun(seq1 = mutant_sequence)
    
    dat <- data.frame(target_seq = rep(c(mutant_sequence,input_sequence),
                                       c(length(candidates_mutant),length(candidates_input))),
                      type_seq = rep(c("Mutant sequence","Input sequence"),
                                     c(length(candidates_mutant),length(candidates_input))),
                      candidates = c(candidates_mutant,candidates_input),
                      candidate_index = rep(1:25,2),
                      mutation_location = c(rep(c("PFS","Protospacer"),c(5,20))),
                      mutation_type = c(rep("M0",25),
                                        rep("M0",5),paste0("M",1:20)),
                      mutation_site = c(rep(mutate_base_to,25),
                                        rep(mutate_base_from,25)))
    complement_RNA_fun <- function(seq1){
      RNA_complement_ref <- setNames(object = c("A","U","G","C"),nm = c("U","A","C","G"))
      sapply(1:length(seq1), function(i){
        seq_c <- unlist(strsplit(seq1[i],""))
        seq_c <- RNA_complement_ref[seq_c] %>% paste0(.,collapse = "")
        return(seq_c)
      })
    }
    dat$mutation_site_com <- complement_RNA_fun(dat$mutation_site)
    dat$mutation_site_crRNA <- paste0(dat$mutation_type,dat$mutation_site_com)
    dat[grep("^M0",dat$mutation_site_crRNA),"mutation_site_crRNA"] <- "M0"
    dat$label <- label
    dat$PFS <- substr(dat$candidates,21,25)
    dat$target <- substr(dat$candidates,1,20)
    seq_reverse_complement_fun <- function(seq1){
      # sequence reverse and complement
      seq1 <- sapply(seq1,function(seq_i){
        seq_i <- rev(unlist(strsplit(seq_i,"")))
        ref_complement <- setNames(c("U","A","C","G"),
                                   c("A","U","G","C"))
        seq_i <- paste0(ref_complement[seq_i],collapse = "")
        return(seq_i)
      })
      return(seq1)
    }
    dat$guide <- seq_reverse_complement_fun(dat$target)
    # seq_T2U_fun <- function(seq1){
    #   seq1 <- sapply(seq1,function(seq_i){
    #     seq_i <- unlist(strsplit(seq_i,""))
    #     ref_T2U <- setNames(c("A","U","G","C"),
    #                                c("A","T","G","C"))
    #     seq_i <- paste0(ref_T2U[seq_i],collapse = "")
    #     return(seq_i)
    #   })
    #   return(seq1)
    # }
    dat$handle <- "AAUUUCUACUAUUGUAGAU"
    # dat$guide_T2U <- seq_T2U_fun(dat$guide)
    dat$crRNA <- paste0(dat$handle,dat$guide)
    
    # 6. activity score
    dat$crRNA_index <- paste0("crRNA",dat$candidate_index)
    dat$crRNA_index <- factor(dat$crRNA_index,levels = unique(dat$crRNA_index))
    ref_PFS_score <- readRDS("./metadata/03_ref_PFS_score.rds")
    # dat$PFS_T2U <- seq_T2U_fun(dat$PFS)
    index_PFS <- match(dat$PFS,ref_PFS_score$PFS_site)
    dat$PFS_activity <- ref_PFS_score[index_PFS,"handle_activity",drop=T]
    
    ref_guide_score <- readRDS("./metadata/03_ref_guide_score.rds")
    index_mismatch <- match(dat$mutation_type,ref_guide_score$sample)
    dat$guide_activity <- ref_guide_score[index_mismatch,"target_activity",drop=T]
    
    # 2024-05-04 add
    ref_MM12314_score <- readRDS("./metadata/07_Cas12a2_crRNA_MM12314_activity_result.rds")
    index_mismatch_mm <- match(dat$mutation_site_crRNA,ref_MM12314_score$sample)
    dat$MM_activity <- ref_MM12314_score[index_mismatch_mm,"target_activity",drop=T]
    dat$guide_activity <- ifelse(is.na(dat$MM_activity),dat$guide_activity,dat$MM_activity)
    
    dat$crRNA_activity <- dat$PFS_activity * dat$guide_activity
    
    crRNA_activity_diff_score <- function(dat){
      res <- data.frame(label = dat$label[1],
                        mutant_sequence = dat[1:25,"target_seq"],
                        input_sequence = dat[26:50,"target_seq"],
                        crRNA = dat[1:25,"crRNA"],
                        crRNA_index =  dat[1:25,"crRNA_index"],
                        activity_input_seq = dat[26:50,"crRNA_activity"],
                        activity_mutant_seq = dat[1:25,"crRNA_activity"]) %>% 
        mutate(activity_diff = activity_mutant_seq - activity_input_seq) %>% 
        mutate(relative_diff = activity_diff / activity_input_seq * activity_mutant_seq)
      sigmod_fun <- function(values){
        sapply(values,function(x){
          return((1/(1+exp(-x))-0.5)*2)
        })
      }
      res$sigmod_diff <- sigmod_fun(res$relative_diff)
      res$diff_score <- res$sigmod_diff * res$activity_mutant_seq^(1/2)
      res$crRNA_index <- factor(res$crRNA_index,levels = unique(res$crRNA_index))
      res$category <- ifelse(res$diff_score > 0.5,"Excellent",
                             ifelse(res$diff_score > 0.15, "Moderate","Imprecise"))
      res[is.na(res$category),"category"] <- "Imprecise"
      res$category <- factor(res$category,levels = c("Excellent","Moderate","Imprecise"))
      return(res)
    }
    
    
    select_col <- c("label","crRNA_index","crRNA","type_seq","mutation_location","PFS",
                    "guide","mutation_type","mutation_site_crRNA","PFS_activity","guide_activity","crRNA_activity")
    res <- dat[,select_col]
    res[26:50,"crRNA"] <- res[1:25,"crRNA"]
    # write.csv(res,file = "./metadata/03_target_crRNA_searching&evaluation.csv")
    
    # diff_score -------------------------------------------------------------------
    dat_diff <- crRNA_activity_diff_score(dat = dat)
    return(list(result=res, dat_diff=dat_diff, seqs_wt_mut=seqs_wt_mut))
  }
  
  if(Cas_protein == "suCas12a2"){
    result_list <- find_crRNA_Cas12a2()
  }
  return(result_list)
}

find_crRNA_multi <- function(parameters_df){
  if(!is.data.frame(parameters_df)){
    stop("Error! The parameters_df should be a data.frame containing columns including input_seq, mutation_site, mutate_to and label.")
  }
  input_seq_m <- parameters_df$input_seq
  mutation_site_m <- parameters_df$mutation_site
  mutate_to_m <- parameters_df$mutate_to
  label_m <- parameters_df$label
  mission_length <- nrow(parameters_df)
  
  # progress bar
  pb <- txtProgressBar(min = 0,max = mission_length,style = 3)
  res <- data.frame()
  dat_diff <- data.frame()
  for (i in 1:mission_length) {
    crRNA_list <- find_crRNA(input_seq = input_seq_m[i],mutation_site = mutation_site_m[i],
                             mutate_to = mutate_to_m[i],label = label_m[i])
    res_tmp <- crRNA_list$result
    dat_diff_tmp <- crRNA_list$dat_diff
    res <- rbind(res,res_tmp)
    dat_diff <- rbind(dat_diff,dat_diff_tmp)
    setTxtProgressBar(pb,i)
  }
  
  res$label <- factor(res$label,levels = label_m)
  dat_diff$label <- factor(dat_diff$label,levels = label_m)
  close(pb)
  return(list(result=res,dat_diff=dat_diff))
}

plot_seq <- function(seq_wt_mut){
  csl <- make_col_scheme(chars = c("A","T", "C", "G"),
                         # groups = c("gr1","gr1", "gr2","gr2"),
                         cols = c("black","black","black","black"))
  p1 <- ggplot()+
    annotate("rect",xmin = 24+0.5,xmax = 25+0.5,ymin=0,ymax=2,fill ="blue",alpha=0.7)+
    geom_logo(seq_wt_mut[1],col_scheme=csl)+
    theme_bw()+
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          panel.grid = element_blank())
  
  p2 <- ggplot()+
    annotate("rect",xmin = 24+0.5,xmax = 25+0.5,ymin=0,ymax=2,fill ="red",alpha=0.7)+
    geom_logo(seq_wt_mut[2],col_scheme=csl)+
    theme_bw()+
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          panel.grid = element_blank())
  return(p1+p2+plot_layout(ncol = 1))
}

plot_crRNA_activity <- function(res,palette=c("#00468b","#ed0000"),palette_rev="No",alpha=1,label=NULL){
  library(ggplot2)
  library(patchwork)
  if(alpha<0){alpha <-  0}
  if(alpha>1){alpha <- 1}
  if(!is.data.frame(res)){
    stop("Error! The res parameter should be a `data.frame`, please use `find_crRNA$result` instead.")
  }
  if(length(palette)<2){
    stop("Error! The palette should contain at least 2 colors.")
  }
  if(is.null(label)){
    label <- unique(res$label)
  }
  # label <- label[1]
  res <- res[res$label %in% label,]
  # plot
  p1 <- ggplot(res,aes(x=type_seq,y=crRNA_activity,fill=type_seq))+
    geom_bar(stat = "identity",alpha=alpha)+
    facet_grid(label~crRNA_index)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          strip.background = element_rect(fill = "#eaeae0"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.background = element_rect(colour = "black"),
          legend.title = element_blank(),
          legend.position = "top",
          plot.title = element_text(hjust = 0.5,size = 15))+
    labs(x=NULL,fill=NULL,y="Estimated crRNA activity")+
    scale_y_continuous(limits = c(0,1),expand = c(0,0))
    
  
  p1 <- switch(palette_rev,
               "Yes" = p1 + scale_fill_manual(values = rev(palette)),
               "No" = p1 + scale_fill_manual(values = palette))
  
  return(p1)
}

plot_crRNA_diff <- function(dat_diff,palette=c("#dc0000","#4dbbd5","#00a087"),
                            palette_rev="No",
                            alpha=1,label=NULL){
  library(ggplot2)
  if(alpha<0){alpha <-  0}
  if(alpha>1){alpha <- 1}
  if(!is.data.frame(dat_diff)){
    stop("Error! The res parameter should be a `data.frame`, please use `find_crRNA$result` instead.")
  }
  if(length(palette)<2){
    stop("Error! The palette should contain at least 2 colors.")
  }
  if(is.null(label)){
    label <- unique(dat_diff$label)
  }
  # dat_diff <- dat_diff[dat_diff$label %in% label,]
  # plot
  p1 <- ggplot(dat_diff,aes(x=crRNA_index,y=diff_score,fill=category))+
    geom_bar(stat = "identity")+
    geom_hline(yintercept = 0.15,lty=2)+
    geom_text(aes(label=sprintf("%.2f", diff_score),y=diff_score+0.015))+
    theme_classic()+
    theme(legend.position = "right",
          legend.background = element_rect(color = "black"),
          axis.text.x = element_text(angle = 45,hjust = 1),
          plot.title = element_text(hjust = 0.5,size = 15))+
    scale_y_continuous(expand = c(0,0),limits = c(0,1.05))+
    labs(x=NULL,y="Ability score to distinguish mutant from input sequence",
         fill="Performance")
  
  p1 <- switch(palette_rev,
               "Yes" = p1 + scale_fill_manual(values = rev(palette)),
               "No" = p1 + scale_fill_manual(values = palette))
  return(p1)
}

plot_crRNA_diff_multi <- function(dat_diff,palette=c("#dc0000","#4dbbd5","#00a087"),
                                  palette_rev="No",
                                  alpha=1,label=NULL){
  library(ggplot2)
  if(alpha<0){alpha <-  0}
  if(alpha>1){alpha <- 1}
  if(!is.data.frame(dat_diff)){
    stop("Error! The res parameter should be a `data.frame`, please use `find_crRNA$result` instead.")
  }
  if(length(palette)<2){
    stop("Error! The palette should contain at least 2 colors.")
  }
  if(is.null(label)){
    label <- unique(dat_diff$label)
  }
  dat_diff <- dat_diff[dat_diff$label %in% label,]
  dat_diff$crRNA_index <- factor(dat_diff$crRNA_index,levels = rev(unique(dat_diff$crRNA_index)))
  # plot
  p1 <- ggplot(dat_diff,aes(x=crRNA_index,y=diff_score,fill=category))+
    geom_bar(stat = "identity")+
    geom_hline(yintercept = 0.15,lty=2)+
    geom_text(aes(label=sprintf("%.2f", diff_score),y=diff_score+0.015))+
    theme_bw()+
    facet_wrap(label~.,scales = "free_y")+
    theme(legend.position = "right",
          legend.background = element_rect(color = "black"),
          # axis.text.x = element_text(angle = 45,hjust = 1),
          plot.title = element_text(hjust = 0.5,size = 15),
          strip.background = element_rect(fill = "#eaeae0"))+
    scale_y_continuous(expand = c(0,0),limits = c(0,1.05))+
    labs(x=NULL,y="Ability to distinguish mutant sequence from input sequence",
         fill="Performance")+
    coord_flip()
  p1 <- switch(palette_rev,
               "Yes" = p1 + scale_fill_manual(values = rev(palette)),
               "No" = p1 + scale_fill_manual(values = palette))
  return(p1)
}


crRNA_primer_design <- function(crRNAs,trans_method=c("T7 RNA polymerase"),
                                Cas_protein=c("suCas12a2")){
  library(dplyr)
  # U2T
  seq_U2T_fun <- function(seq1){
    seq1 <- sapply(seq1,function(seq_i){
      seq_i <- unlist(strsplit(seq_i,""))
      ref_T2U <- setNames(c("A","T","G","C"),
                          c("A","U","G","C"))
      seq_i <- paste0(ref_T2U[seq_i],collapse = "")
      return(seq_i)
    })
    return(seq1)
  }
  
  crRNA_U2T <- seq_U2T_fun(crRNAs)
  dat <- NULL
  if(Cas_protein == "suCas12a2"){
    if(trans_method == "T7 RNA polymerase"){
      dat <- data.frame(index = 1:length(crRNAs),
                        crRNA = crRNAs,
                        crRNA_U2T = crRNA_U2T)
      handle_U2T <- unlist(strsplit(dat$crRNA_U2T[1],""))[1:19] %>% paste0(.,collapse = "")
      T7_F <- "TAATACGACTCACTATAGGG"
      primer_F <- paste0(T7_F,handle_U2T,collapse = "")
      dat$primer_Fu <- primer_F
      
      primer_R_design_fun <- function(crRNA_U2T){
        primer_R <- sapply(crRNA_U2T,function(seq_i){
          seq_i <- unlist(strsplit(seq_i,""))
          ref_complement <- setNames(c("T","A","C","G"),
                                     c("A","T","G","C"))
          seq_i <- paste0(c(ref_complement[rev(seq_i)],"C"),collapse = "")
          return(seq_i)
        })
        return(primer_R)
      }
      dat$primer_R <- primer_R_design_fun(crRNA_U2T = dat$crRNA_U2T)
      rownames(dat) <- 1:nrow(dat)
      dat <- dat[,c("crRNA","primer_Fu","primer_R")]
    }
  }
  return(dat)
}