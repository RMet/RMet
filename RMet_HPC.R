# Before running this code you should input the following data
#.................................................................
#A) Input the directory in which current code is copied to
current_dir <-                # e.g. /home/RMet_users/Dr.____ 
#B) Input number of metabolites & samples:
Number_of_metabolites <-      #Input just a number: e.g. 24
Number_of_samples     <-      #Input just a number: e.g. 8
#C) Input which constain to impose:
# warning: you should type the answer with capital letters  
Non_negativity <- TRUE        #Input TRUE or FALSE: e.g. TRUE 
Unimodality    <- TRUE        #Input TRUE or FALSE: e.g. FALSE 
Closure        <- TRUE        #Input TRUE or FALSE: e.g. FALSE 
#.................................................................
# Now you can sumbit the code to your HPC!
#_________________________________________________________________
#_________________________________________________________________
#_________________________________________________________________
#_________________________________________________________________
#_________________________________________________________________
#_________________________________________________________________
#_________________________________________________________________
#_________________________________________________________________
#_________________________________________________________________
#_________________________________________________________________
#_________________________________________________________________
#_________________________________________________________________
#_________________________________________________________________
#_________________________________________________________________
#_________________________________________________________________
#_________________________________________________________________
#_________________________________________________________________
#_________________________________________________________________
#_________________________________________________________________
#_________________________________________________________________
#_________________________________________________________________
#_________________________________________________________________
#_________________________________________________________________
#_________________________________________________________________
#_________________________________________________________________
#_________________________________________________________________
#_________________________________________________________________
#_________________________________________________________________
#_________________________________________________________________
#_________________________________________________________________
#_________________________________________________________________
#_________________________________________________________________
#_________________________________________________________________
#_________________________________________________________________
#_________________________________________________________________


n_comp <- Number_of_metabolites
mydata <-  load(file = paste(current_dir,"/HPC_intial.Rdata" , sep = "") ) 
saved_data <- mydata
#  n_comp <- as.integer(readline("input number of components: "))
#  last_comp <- as.numeric(readline(("input number of components for last analysi: ")))
last_comp <<- 1000

most_dis <- function(mydata,ref_data){
  row_numb <<- nnrow <- as.integer(nrow(mydata)) + 1 -1 
  det_list <- c()
  for (i in c(1:row_numb)){
    
    s_data <- data.matrix(matrix(unlist(mydata[i,]),nrow = 1, byrow=T))
    r_data <- rBind(ref_data,s_data)  
    d_data <- r_data %*% t(r_data)
    d_det <- det(d_data)
    det_list[length(det_list)+1 ] <- d_det
    
  }
  
  sort_det_list <- sort(det_list,decreasing = TRUE, index.return = TRUE)
  #sort_det_list_val <- sort(det_list)
  sort_index <- sort_det_list$ix
  max_index <- which.max(det_list) + 1 - 1
  
  #max_spec <- data.matrix(matrix(unlist(mydata[max_index,]),nrow = 1, byrow=T))
  return(list(first = max_index, second = sort_index))
  
}
mid_ref <- data.matrix(matrix(unlist(colMeans(mydata)),nrow = 1, byrow=T))
#main_func <- most_dis(mydata,mid_ref)
#max_index <- main_func$first
#sort_index0 <- main_func$second
max_index_list = c()
time_list = c()
#sort2 <- length(sort_index0) + 1 - 1
#sort1 <- as.integer(round(sort2 * 0.9)) +1 - 1
#if (sort1 > 1000){
#  mydata <- mydata[-sort_index0[sort1:sort2],]
#opa2
row_numb2 <- nnrow <- as.integer(nrow(mydata)) + 1 -1 
del_p <- (last_comp/row_numb2)**(1/n_comp)

for (i in c(1:n_comp)){
  #cat("\014")
  t1 <- proc.time()
  if (i == 1){
    ref_data <- mid_ref
    main_func <- most_dis(mydata, ref_data)
    max_index <- main_func$first
    max_index_list[length(max_index_list)+1 ] <- max_index
    ref_data <- data.matrix(matrix(unlist(mydata[max_index,]),nrow = 1, byrow=T))
    sort_index <- main_func$second
    sort2 <- length(sort_index) + 1 - 1
    sort1 <- as.integer(round(sort2 * del_p)) +1 - 1
    
    if (sort1 > 1000){
      mydata <- mydata[-sort_index[sort1 : (sort2-i)],]
    }
  }
  else {
    main_func <- most_dis(mydata, ref_data)
    max_index <- main_func$first
    max_index_list[length(max_index_list)+1 ] <- max_index
    ref_data <- rBind(ref_data, data.matrix(matrix(unlist(mydata[max_index,]),nrow = 1, byrow=T)))
    sort_index <- main_func$second
    sort2 <- length(sort_index) + 1 - 1
    sort1 <- as.integer(round(sort2 * del_p)) + 1 -1
    
    if (sort1 > 1000 ){
      mydata <- mydata[-sort_index[sort1:(sort2-i)],]
    }
  }
  
  svalue(win5_lab14) <- paste("Generating initial estimiation: ",round(100 * i / n_comp) , "%" ," Completed." )
  #    svalue(win5_lab9)  <<- "Generating initial estimiation:"
  #    svalue(win5_lab10)  <<- paste(round(100 * i / n_comp) , "%")
  #    svalue(win5_lab11)  <<- "Completed"
  #    font(win5_lab9) <<- list( weight = "bold" , color = "blue",size = 14)
  #    font(win5_lab10) <<- list(   color = "blue",size = 22)
  #    font(win5_lab11) <<- list( weight = "bold" , color = "blue",size = 14)
  
}
svalue(win5_lab14) <- paste("Initiating the MCR..." )

max_index_data <- data.matrix(max_index_list)

s_opa = c()

row_numb3 <-as.integer(nrow(saved_data)) + 1 -1
sum_list <<- c()
ref_data <<- ref_data


main_index <<- c()
row_numb4 <- as.integer(nrow(ref_data)) +1 -1

sum_row <- rowSums(saved_data)
ref_row <- rowSums(ref_data)


for ( i in ref_row){
  
  #  print(dim(ref_data))
  #    check_sum <- round(sum(ref_data[i,]),2)
  
  this_check <- which(sum_row == i)
  #  print(this_check)
  
  main_index[length(main_index) + 1] <- this_check
}



main_index_data <- data.matrix(main_index)


for (j in main_index_data){
  s_opa <- rBind(s_opa,saved_data[j,])
}

s_opa <<- data.matrix(s_opa)
#assign(paste("sopt",svalue(win5_aug_combo) , sep = "_"), s_opa, envir = .GlobalEnv)
assign(svalue(win5_con_name2), s_opa, envir = .GlobalEnv)

#  is_save <- readline("want to save opa result?(y/n) ")

#  is_t_plot <- readline("want process time plot?(y/n) ")
#  if (is_t_plot == "y"){
#    plot(time_list)
#  }


#  raw_data <- mydata 
s_estimate <- s_opa 

if ((nrow(saved_data) %% first_col_hash[[svalue(win5_aug_combo)]]) != 0){
  cur_x <<- first_col_hash[[svalue(win5_aug_combo)]]
  cur_nrow <<- nrow(saved_data)
  excess_x <<- cur_nrow - floor(cur_nrow / cur_x) * cur_x
  
  if (excess_x %% 2 == 0){
    haf <- excess_x / 2
    del_x <- c(1 : haf)
    append( del_x , c(cur_nrow - haf + 1 :cur_nrow))
    saved_data <<- saved_data[ - del_x , ]
    
  }
  else{
    haf <- floor(excess_x / 2)
    del_x <- c(1 : haf + 1)
    append( del_x , c(cur_nrow - haf + 1 :cur_nrow))
    saved_data <<- saved_data[ - del_x , ]
  }
  
  
}



raw_c_ini <<- t(ginv(t(s_estimate)) %*% t(saved_data))
n_mat <<- as.integer(seg_mat_hash[[svalue(win5_aug_combo)]])  
#as.integer(readline("input Number of matrices: ")) + 1 -1

data_row <- as.integer(nrow(raw_c_ini)) +1 - 1
data_col <- as.integer(ncol(raw_c_ini)) +1 - 1
s_row <<- as.integer(data_row / n_mat) +1 -1
ini_esti_c_list <<- list()
for (i in c(1:n_mat)){
  
  ini_numb <- (i-1)*s_row  + 1
  end_numb <- i * s_row
  ex_c_mydata <- raw_c_ini[ini_numb:end_numb,]
  
  ini_esti_c_list <<- append(ini_esti_c_list,list(ex_c_mydata))
  
}


mydata_list = list()

for (i in c(1:n_mat)){
  
  ini_numb <- (i-1)*s_row  + 1
  end_numb <- i * s_row
  ex_mydata <- saved_data[ini_numb:end_numb,]
  mydata_list <- append(mydata_list,list(ex_mydata))
  
  
}


#svalue(win5_lab9) <<- " "
#svalue(win5_lab10) <<- "  "
#svalue(win5_lab11) <<- "  "
#  t11 <- proc.time()
test0 <- saeed_als(CList=ini_esti_c_list,S=matrix(1,nrow= dim(saved_data)[2],ncol= svalue(win5_comp_numb )),
                   PsiList=mydata_list,normS = 2,optS1st=TRUE,thresh = as.numeric(svalue(win5_conv))  , nonnegC = TRUE,nonnegS = TRUE,maxiter = svalue(win5_svd_spin) , uniC = is_uni , uniS = is_uni )
#  t22 <- proc.time()
#  print(t22 - t11)

first_r_mcr_s <<- test0$S

svalue(win5_lab14) <<- "Generating profiles...."
num_seg <<- svalue(win5_samp_numb)

c_real <- t(ginv(first_r_mcr_s) %*% t(saved_data))
assign(svalue(win5_con_name1) , c_real , envir = .GlobalEnv)

c_real_row <<- nrow(c_real)
sample_row <- c_real_row / num_seg





for ( k in c(1:num_seg)){
  
  uaug <- c_real[(c(1 + (k-1)*sample_row) : (k*sample_row) ) ,]
  c_nrow <- as.numeric(nrow(uaug)) 
  c_ncol <- as.numeric(ncol(uaug)) 
  
  
  nnmat <- n_mat / num_seg
  m_nrow <- c_nrow / nnmat
  #cut_edge <- c_nrow / num_seg
  
  #c_ncol
  u1 = c()
  u2 = c()
  a1 <<- c()
  a2 <<- c()
  
  
  for( j in c(1: c_ncol)){
    
    m <- matrix(uaug[,j] , nrow = m_nrow , ncol = nnmat)
    m_svd <- svd(m , nu = m_nrow , nv = nnmat)
    u_m <- m_svd$u
    u_v <- m_svd$v
    u_s <- m_svd$d
    u1 <<- cbind(u1 , u_m[,1])
    u2 <<- cbind(u2 , as.vector(t(u_v[1, ])))
    a1 <<- cbind(a1 ,as.vector( t(colMeans(t(m)))) )
    a2 <<- cbind(a2 ,as.vector( t(colMeans(m))) )
    #  u2 <- t(u2)
    
  }
  
  
  first_col_con <<- paste("First_col","S",k , svalue(win5_aug_combo) , sep = "_" )
  second_col_con <<- paste("Sec_col","S",k ,svalue(win5_aug_combo), sep = "_" )
  
  sec_col_mat_list[length(sec_col_mat_list) + 1] <<- second_col_con
  fir_col_mat_list[length(fir_col_mat_list) + 1] <<- first_col_con
  
  assign(first_col_con, a1, envir = .GlobalEnv)
  assign(second_col_con, a2, envir = .GlobalEnv)  
  
}


svalue(win5_lab14) <<- "Analysis is complete."


for (t in fir_col_mat_list){
  s_tab[length(s_tab) + 1] <<- t
}

for (t in sec_col_mat_list){
  s_tab[length(s_tab) + 1] <<- t
}

s_tab[length(s_tab) + 1] <<- svalue(win5_con_name1)
s_tab[length(s_tab) + 1] <<- svalue(win5_con_name2)

s_tab <<- rev(s_tab)

})


addHandlerChanged(win5_next_but, function(h,...){
  # dev.off()
  visible(window5) <<- FALSE
  visible(window6) <<- TRUE
  plot(0,type = "n",xlab = "Retention time", ylab = "Intensity", main = "Chromatogram")
  
  prof_combo[] <<- s_tab
  win6_comp_combo[] <<- c(1: svalue(win5_comp_numb))
  
  
})
