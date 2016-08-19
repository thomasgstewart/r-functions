impute_boot <- function(
  data,
  g
){
# THE IMPUTE THEN BOOT FUNCTION
}

parameter_plot <- function(model, ...){
  UseMethod("parameter_plot", model)
}

parameter_plot.lm <- function(model, intercept = FALSE, ...){
  data <- data.frame(
    parameter = names(coef(model)),
    point = coef(model),
    B     = confint(model),
    stringsAsFactors = FALSE
  )
  names(data)[3:4] <- c("lb","ub") 
  
  if(!intercept) data <- data[-1, , drop = FALSE]
    
  plot_point_ci(
    point = 2,
    lb = 3,
    ub = 4,
    data = data,
    labels = 1,
    vertical_lines = 0,
    ...
  )
  invisible(data)
}

plot_point_ci <- function(
  point,
  lb,
  ub,
  data,
  labels,
  vertical_lines = NULL,
  col = NULL,
  ...
){
# BIG PICTURE: GENERATE A PLOT OF POINT ESTIMATES AND CONFIDENCE INTERVALS
  
  data <- data[nrow(data):1, , drop = FALSE]
  
  xlim <- range(data[,c(lb,ub)])
  ylim <- c(1,nrow(data))
  plot.new()
  plot.window(xlim = xlim, ylim = ylim + c(-.5, .5 + 2*strheight("X")))
  if(!is.null(vertical_lines)){
    abline(v = vertical_lines, col = "gray80")
  }
  if(is.null(col)){
    cols <- rep("black", nrow(data))
  }else{
    cols <- data[,col]
  }
  
  for(i in 1:nrow(data)){
    lines(xlim, rep(i,2), lty = 2, col = "gray80")
    lines(data[i,c(lb,ub)],rep(i,2), col = cols[i], ...)
    points(data[i,point], i, col = cols[i], ...)
  }

  text(
    rep(xlim[1],nrow(data)), 
    strheight("X") + 1:nrow(data), 
    data[,labels], 
    pos = 4
  )
  axis(1)
  box()  
}


capwords <- function(
  s, 
  strict = FALSE
){
# FROM R HELP
  cap <- function(s){
    paste(
      toupper(substring(s, 1, 1)),
      {s <- substring(s, 2); if(strict) tolower(s) else s},
      sep = "", 
      collapse = " " 
    )
  } 
  sapply(
    strsplit(s, split = " "), 
    cap, 
    USE.NAMES = !is.null(names(s))
  )
}

labelify <- function(x){
  x <- gsub("_", " ", x) # REPLACE _ WITH " "
  x <- capwords(x)       # CAPITALIZE WORDS
  return(x)
}


nona <- function(x) return(x[!is.na(x)])

boot_impute <- function(
  data, 
  boots = 100, 
  imputes = 1,
  strata = NULL,
  seed = 123,
  g,
  g_packages,
  ...
){
# OUTPUT = a list of lists (dimensions = boots + 1 X imputes)
  #
  # Each row corresponds to a bootstrapped sample, each
  # column corresponds to an imputation sample.  The 
  # boots + 1 row is the outcome when the input data set
  # is not bootstrapped sampled, i.e., it is the original
  # dataset.
  #
  # The contents of each cell is the output 
  # of the function g
  #
  # Bootstrapping is done in parallel (requires foreach and
  # doParallel). Function will use all cores except one for
  # calculations.
  #
  # Imputation is done via mice package
  #
# INPUTS
  # data    = data frame to be bootstrapped than imputed
  # boots   = integer, number of bootstrap samples
  # imputes = integer, number of imputed datasets
  # strata  = 3 options:
  #           NA - do not bootstrap samples (original dataset 
  #                is sent to each "row")
  #           NULL - bootstrap samples
  #           character string of variable name - when
  #             bootstrapping, sample independently 
  #             from within each strata level of variable.
  # seed    = random number generator seed
  # g       = function which calculates statistics 
  #           to be bootstrapped and imputed.  Required 
  #           inputs are: 
  #           x = a data.frame
  #           observed_call = TRUE/FALSE.  
  # g_packages = packages needed by R to perform the
  #              calculations in g
  # ... = arguments to g
  

  
  require(mice)
  require(foreach)
  require(doParallel)
  set.seed(seed)
  N <- nrow(data)
  
  if(boots > 0){
    cores_2_use <- detectCores() - 1
    cl <- makeCluster(cores_2_use)
    clusterSetRNGStream(cl, seed)
    registerDoParallel(cl)
    out <- foreach(
      no = 1:boots, 
      .packages = g_packages
    )%dopar%{
      if(is.null(strata)){
        idx <- sample.int(N, replace = TRUE)
      }else if(is.na(strata)){
        idx <- 1:nrow(data)
      }else{
        strata_levels <- unique(strata)
        id <- 1:nrow(data)
        idx <- NULL
        for(s in strata_levels){
          N_s <- sum(strata == s)
          idx <- c(idx, id[strata == s][sample.int(N_s, replace = TRUE)])
        }
      }
      DF <- data[idx,]
      
      # ONLY DO MI IF MISSING DATA IS AVAILABLE
      if(sum(is.na(DF))){
        tempData <- mice(
          data      = DF, 
          m         = imputes, 
          maxit     = 5,
          meth      = 'pmm',
          #seed      = seed,
          printFlag = FALSE
        )
        long <- complete(tempData, 'long')
        outline <- by(long, long[,'.imp'], g, observed_call = FALSE, ...)
        outline
      }else{
        outline <- list(g(x = DF, observed_call = FALSE, ...))
        outline
      }
    }
    stopCluster(cl)
  }else{
    out <- list()
  }
  
  # Calcuation with observed data
  # ONLY DO MI IF MISSING DATA IS AVAILABLE
  if(sum(is.na(data))){
    tempData <- mice(
      data      = data, 
      m         = imputes, 
      maxit     = 5,
      meth      = 'pmm',
      seed      = seed,
      printFlag = FALSE)
    long <- complete(tempData, 'long')
    outline <- by(long, long[,'.imp'], g, observed_call = TRUE, ...)
  }else{
    outline <- list(g(x = data, observed_call = FALSE, ...))
  }
  out[[length(out)+1]] <- outline
  return(out)  
}




missingness_info <- function(
  data, 
  upper_limit = 5, 
  max_vars    = 5
){
# OUTPUT: a list containing 2 data frames.  
#
#   The first data frame gives the frequency
#   distribution on the number of missing variables.
# 
#   The second data frame lists the number of missing
#   observations for each variable.
#
# INPUTS:
#   data        = data.frame for which a missing data 
#                 report will be generated
#   upper_limit = number.  The right tail of the frequency 
#                 distribution reported in the first 
#                 data frame is truncated to 
#                 [upper_limit, Inf) and is displayed as
#                 "upper_limit +".
#   max_vars    = number.  Limits the list of variables
#                 reported in the second data frame to 
#                 the first max_vars most frequently 
#                 missing variables.  If there are multiple
#                 variables with the same number of missing
#                 values, all such the variables will 
#                 be reported. (This means more than 
#                 max_vars variables can appear in the output)
  
  missing_count <- apply(is.na(data), 1, sum)
  missing_count_cat <- cut(
    missing_count, 
    breaks = c(-1:upper_limit, Inf),
    labels = c(0:upper_limit, upper_limit %|% "+" )
  )
  out <- xtabs(~missing_count_cat)
  
  count <- as.numeric(out)
  pct <- count/sum(count)*100
  cumpct <- cumsum(pct)
  out1 <- data.frame(
    `Number of Missing Variables<br>(within a participant)` = names(out),
    `Frequency`         = as.numeric(out),
    `Percent`           = sprintf("%5.2f", pct),
    `Cumulative Percent`= sprintf("%5.2f", cumpct),
    stringsAsFactors    = FALSE,
    check.names         = FALSE
  )
  is100 <- cumpct > 100 - 1e-15
  keep_row <- !is100 | !duplicated(is100)
  out1 <- out1[keep_row,,drop = FALSE]
  
  var_missing_count <- apply(is.na(data),2, sum)
  max_var_count <- sort(
    var_missing_count, 
    decreasing = TRUE)
  cutoff <- max_var_count[max_vars]
  max_var_count <- max_var_count[max_var_count >= cutoff]
  max_var_pct <- 100*as.numeric(max_var_count) / nrow(data)
  out2 <- data.frame(
    Variable = names(max_var_count),
    `Obs Available` = nrow(data) - as.numeric(max_var_count),
    `Percent Available` = sprintf("%5.2f",100 - max_var_pct),
    `Obs Missing` = as.numeric(max_var_count),
    `Percent Missing`   = sprintf("%5.2f", max_var_pct),
    stringsAsFactors    = FALSE,
    check.names         = FALSE
  )
  return(list(out1,out2))
}

get_new_data <- function(lm1){
    MF <- model.frame(lm1)
    ff <- formula(lm1)

    MF[[1]] <- NULL
    for(i in seq_along(MF)){
      if("rms" %in% class(MF[[i]])){
        knots <- attr(MF[[i]],"parms")
        new_ff <- as.formula(".~ + rcs(x = " %|% colnames(MF[[i]])[1] %|% ", parms = c(" %|% paste(knots, collapse = ",") %|% ")) -rcs(" %|% colnames(MF[[i]])[1] %|% ", " %|% length(knots) %|% ") + .")
        ff <- update.formula(ff,new_ff)
        names(MF)[i] <- colnames(MF[[i]])[1]
        MF[[i]] <- MF[[i]][,1] 
      }
    }
    out <- list(formula = ff, DF = MF)
    return(out)
  }

formatp <- function(
  x, 
  digits, 
  sig = 0.05, 
  sig_marker = c("","")
){
# CUSTOM p-value FORMAT FUNCTION
#
# OUTPUT: a character string.  
#
#   The user specifies the number of digits to report.
#   The input p-value is formated so that values
#   which round to 0 are reported as "< 0.001" or 
#   "< 0.0001" depending on the number of digits.
#
#   The user can also specify significance markers
#   to place before or after the formatted p-value
#   if the p-value is below the significance
#   threshold.  For example, if sig = 0.05 and p = 0.01,
#   the output could be "<b>0.01</b>" to bold the p-value
#   in HTML output.
#
# INPUTS:
#   x          = number to be formatted
#   digits     = number of digits to report.
#   sig        = number.  The significance threshold.
#   sig_marker = character vector of length 2.  First element
#                is the place before the formatted value.
#                The second element is placed after the 
#                formatted value.  Example: c("<b>","</b>")
    
  p <- if(round(x,digits) < 1/(10^digits)){
    "< " %|% (1/(10^digits))
  }else{
    sprintf("%4." %|% digits %|% "f", x)
  }
  
  if(x < sig) p <- sig_marker[1] %|% p %|% sig_marker[2]
  
  return(p)
}



explicit_rcs <- function(
  data, 
  formula
){
#
# BIG PICTURE: CHANGE 
#   y ~ rcs(x1, 3) + rcs(x2, 4) TO
#   y ~ rcs(x1, parms = c(1, 2, 3)) + rcs(x2, parms = c(10, 15, 16, 27))
#
# FUNCTION WILL MORPH A FORMULA WHICH INCLUDES RESTRICTED 
# CUBIC SPLINES PARAMETERIZED WITH THE NUMBER OF KNOTS
# TO A FORMULA IN WHICH THE KNOT PLACEMENT IS EXPLICIT
#
# WHY:
#   1. EXPLICIT KNOTS ALLOW THE FORMULA TO BE USED WITH NON-HMISC 
#      FUNCTIONS, LIKE predict.lm IN THE stats PACKAGE
#   2. ENSURE THAT KNOTS ARE STABLE
#
# OUTPUT: a formula in the form described above.
#
# INPUTS:
#  formula - a formula with rcs items to make explicit
#  data    - a data.frame from which the knots will be calculated

  ff <- gsub(" ", "", deparse(formula[[3]], width.cutoff = 500L))
  while(length(grep("^.*rcs\\(([[:alnum:]_\\.]+),[0-9]+\\).*$",ff))){
    variable <- sub("^.*rcs\\(([[:alnum:]_\\.]+),[0-9]+\\).*$", "\\1", ff)
    nknots <- sub("^.*rcs\\([[:alnum:]_\\.]+,([0-9]+)\\).*$", "\\1", ff)
    rcs_x <- Hmisc:::rcspline.eval(x  = data[,variable],
                                   nk = as.numeric(nknots))
    knot_locations <- attr(rcs_x,"knots")
    ff <- sub("rcs\\(" %|% variable %|% "," %|% nknots %|% "\\)", 
              "rcs\\(" %|% variable %|% ", parms = c(" %|% paste(knot_locations, collapse = ",") %|% ")\\)",
              ff)
  }
  formula <- as.formula(deparse(formula[[2]]) %|% "~" %|% ff)
  return(formula)
}

npct <- function(x) sprintf("%i (%2.0f%%)", sum(x), mean(x)*100)

character_to_numeric <- function(data, keep_character = NULL){
# 
# BIG PICTURE: CONVERT COLUMNS OF A data.frame THAT CAN BE NUMERIC
#              TO NUMERIC
#
# OUTPUT: A data.frame WITH COLUMNS CONVERTED TO NUMERIC IF POSSIBLE
#
# INPUTS: 
#   data           - data.frame
#   keep_character - vector of positions or variable names which
#                    should not be converted, even if potentially 
#                    numeric.  For example, study_id.
  
  # IDENTIFY COLUMNS TO CHECK
  if(is.character(keep_character)) keep_character <- which(names(data) %in% keep_character)
  columns_to_check <- setdiff(seq_along(data), keep_character)
  
  for(i in columns_to_check){
    
    # ONLY CONSIDER CHARACTER COLUMNS
    if(!is.character(data[,i])) next

    # MARK "NA" as NA
    na_idx <- data[ , i] %in% "NA"
    data[na_idx, i] <- NA_character_
    
    # CHECK THAT ALL ENTRIES ARE NUMERIC, i.e, numbers.numbers or numbers
    numeric_idx <- grep("^[ ]*[0-9]*[\\.]{0,1}[0-9]*[ ]*$", data[, i])
    
    # IF NUMERIC, CONVERT TO NUMERIC
    if(length(numeric_idx)==nrow(data)){
      tryCatch(
        expr    = {data[ , i] <- as.numeric(data[ , i])}, 
        warning = function(w) cat("problem column: ", i, "\n")
      )  
    } 
  }
  return(data)
}


lineplot_ci <- function(
  x, 
  y, 
  y_lb, 
  y_ub, 
  data     = NULL, 
  add      = FALSE,
  line_col = "darkblue",
  ci_col   = "gray80",
  type     = "l",  
  xlim     = NULL, 
  ylim     = NULL,
  log      = "", 
  xlab     = NULL, 
  ylab     = NULL,
  ann      = par("ann"), 
  axes     = TRUE, 
  frame.plot  = axes,
  panel.first = grid(), 
  panel.last  = NULL, 
  asp      = NA, 
  ...
){
#
# BIG PICTURE: GENERATE A LINE PLOT WITH CONFIDENCE BAND
#
# OUTPUT: See big picture
#  
# INPUTS:
#
#  x        - character string or number of x-axis variable in data
#  y        - character string or number of y-axis variable in data
#  y_lb     - "         "      "  "      "  y lower bound   "  "
#  x_lb     - "         "      "  "      "  y upper bound   "  "
#  data     - data.frame with x, y, y_lb, x_lb
#  
#  -OR-
#  
#  x        - numeric vector, x-axis values
#  y        - "       "     , y-axis values
#  y_lb     - "       "     , y lower bound
#  x_lb     - "       "     , y upper bound
#  data     - NULL
#
# -OTHER ARGUMENTS-
#  
#  add      - logical. TRUE - add to previous plot, FALSE - generate new plot
#  line_col - color of mean line.
#  ci_col   - color of confidence band
#  type     - type of mean line
#  see list of arguments
  
  if(!is.null(data)){
    x <- data[, x]
    y <- data[, y]
    y_lb <- data[, y_lb]
    y_ub <- data[, y_ub]
  }
  #browser()
  if(is.null(ylim)) ylim <- range(y_lb,y_ub)
  if(is.null(xlim)) xlim <- range(x)
  if(is.null(xlab)) xlab <- ""
  if(is.null(ylab)) ylab <- ""
  if(!add){
    plot.new()
    plot.window(ylim = ylim, xlim = xlim, log = log, asp = asp, ...)
  }
  if(!is.null(panel.first)) panel.first
  polygon(c(x, rev(x)),c(y_lb, rev(y_ub)), col = ci_col, border = ci_col)
  lines(x, y, col=line_col, lwd=3, type = type)
  if(axes){
    axis(1,...)
    axis(2,...)
    title(xlab = xlab, ylab = ylab)
  }
  if(frame.plot) box()
  if(!is.null(panel.last)) panel.last
}

`%|%` <- function(a,b) paste0(a,b)

findvarname <- function(var, data = NULL){
  if(is.null(data) | !is.data.frame(data)){
    message("Please specify a dataset")
  }else{
    return(sort(names(data)[grep(var,names(data),ignore.case=TRUE)]))
  }
}

# ALIAS FOR findvarname
fvn <- function(var, data = NULL) findvarname(var, data)


namify <- function(x){
  # CHANGE VARIABLE NAMES OR CHARACTER STRINGS TO BE
  # STYLE GUIDE APPROPRIATE VARIABLE NAMES
  # SEE namify.default FOR LIST OF CHANGES
  UseMethod("namify",x)
}

namify.default <- function(txt){
  txt <- gsub("\\.", "_", txt)          # Change . to _
  txt <- gsub(" +$|^ +|\\(|\\)","",txt) # Remove leading/trailing spaces, (, )
  txt <- gsub("%","pct",txt)            # Change % to pct
  txt <- gsub(" ","_",txt)              # Change interior space to underscore
  txt <- tolower(txt)                   # to lower case
  txt <- gsub("[^_[:alnum:]]","",txt)   # Remove punctuation
  return(txt)
}

namify.data.frame <- function(x){
  names(x) <- namify(names(x))
  return(x)
}


summaryM_to_df <- function(tbl2,html_space=TRUE,...){
o <- capture.output(print(tbl2,...))
o <- o[-grep("+---",o)]
o <- o[grep("^\\|",o)]
o <- gsub("^\\||\\|$","",o)
if(html_space) o <- gsub(" ","&nbsp;",o)
nc <- length(strsplit(o[1],"\\|")[[1]])
out <- as.data.frame(array(NA_character_,dim=c(length(o),nc)),stringsAsFactors=FALSE)

for(i in 1:length(o)) out[i,] <- strsplit(o[i],"\\|")[[1]]
names(out) <- out[1,]
out <- out[-1,]
rownames(out) <- NULL
return(out)
}


###

cross_tabs <- function(df,row_vars,col_vars,cell_fun,file="",title="Table Title",label="tbl",column_head,cell_h_align="Z",table_width=1, footnote=" "){
  
  #Change col_vars, row_vars to numeric vector if input is variables names
  if(!is.numeric(col_vars)){tf <- as.data.frame(matrix(1:ncol(df),nrow=1)); colnames(tf) <- colnames(df); col_vars <- as.numeric(tf[1,col_vars]) }
  if(!is.numeric(row_vars)){tf <- as.data.frame(matrix(1:ncol(df),nrow=1)); colnames(tf) <- colnames(df); row_vars <- as.numeric(tf[1,row_vars]) }
  
  #Change _ to \_ in variable names and factor levels
  colnames(df) <- gsub("_","\\\\_",colnames(df))
  
  
  #Metadata about table
  md <- list()
  md[['n_cols']] <- rep(0,length(col_vars))
  md[['col_reverse_lookup']] <- function(x){ which(col_vars == x) }
  md[['row_reverse_lookup']] <- function(x){ which(row_vars == x) }
  
  for(i in col_vars){
    if(!is.factor(df[,i])) df[,i] <- factor(df[,i])
    md[['n_cols']][ md[['col_reverse_lookup']](i) ] <- length(levels(df[,i]))
    md[['col_levels']][[ md[['col_reverse_lookup']](i) ]] <- levels(df[,i])
  }
  
  if(is.list(md[['col_levels']])){  md[['col_combinations']] <- do.call(expand.grid,md[['col_levels']])}
  else{ md[['col_combinations']] <- data.frame(a=md[['col_levels']]) }
    
  for(j in 1:ncol( md[['col_combinations']] )){  md[['col_combinations']][,j] <-  as.character(md[['col_combinations']][,j]) }
  md[['ncol_df']] <- ncol(df)
  
  
  df[ 1:prod(md[['n_cols']]) + md[['ncol_df']] ] <- " "
  
  df$idx <- apply(df[,row_vars,drop=FALSE],1,paste,collapse="|")
  df$idy <- apply(df[,col_vars,drop=FALSE],1,paste,collapse="|")
  
  df2 <-by(df,df$idx,function(x){
    #browser()
    oo <- x[1,unique(c(row_vars,md[['ncol_df']]+1:prod(md[['n_cols']])))]
    for(j in 1:nrow(md[['col_combinations']])){
      idz <- x$idy==paste(as.character(md[['col_combinations']][j,]),collapse="|")
      if(sum(idz)==0){next}
      z <- x[idz,]
      oo[1,j+length(unique(c(row_vars)))] <- cell_fun(z)
    }
    return(oo)
  })
  df3 <- do.call(rbind,df2)
  out <- df3
  
  ###
  ### Format table
  ###
  
  # Numeric to Character
  for(j in row_vars){ df3[,md[['row_reverse_lookup']](j)] <- as.character(df3[,md[['row_reverse_lookup']](j)]); 
                      df3[is.na(df3[,md[['row_reverse_lookup']](j)]),md[['row_reverse_lookup']](j)] <- "Missing"
  }
  
  # Remove repeated row labels, add horizontal lines, add optional column heading
  #df3 <- rbind(colnames(df3),df3)
  df3$lines <- Inf

  for(i in nrow(df3):2){
    if(i == 1){break}
    for(j in row_vars){
      if( df3[i,md[['row_reverse_lookup']](j)] == df3[i-1,md[['row_reverse_lookup']](j)] ){ df3[i,md[['row_reverse_lookup']](j)] <- " "}
      else{
        df3[i-1,'lines'] <- min(df3[i-1,'lines'],md[['row_reverse_lookup']](j))
      }
    }
  }
  
  df3$lines2 <- paste("\\cmidrule(l){",df3$lines,"-",ncol(df3)-1,"}",sep="")
  df3$lines2[df3$lines == Inf] <- "\\bottomrule"
  df3$lines2[df3$lines == md[['row_reverse_lookup']](rev(row_vars)[1])] <- " " 
  
  #Body of Table
  df3$lines <- NULL
  # Optional Column Heading
  if(!is.null(column_head)){df3 <- rbind(c(rep("\\rule[-.95em]{0em}{1em} ",length(row_vars)),rep(column_head,prod(md[['n_cols']])),"\\rule{0pt}{0pt}%\\\\\n "),df3) }
  df3[,ncol(df3)-1] <- paste(df3[,ncol(df3)-1],"\\\\\n",df3$lines2,sep="")
  df3$lines2 <- NULL
  
  
  
  
  lt1 <- paste(apply(df3,1,paste,collapse="&"),collapse="")
  
  #Header of Table
  ht <- df3[rep(1,length(col_vars)),]
  ht$lines <- " "
  
  header <- list()
  lines <- list()
  for(i in 1:length(col_vars)){
    if(i == length(col_vars)){   
      header[[i]] <- paste("\\multicolumn{1}{c}{", levels(df[,col_vars[i]]) ,"}",sep="")
      lines[[i]] <- 1
    }
    else{
      header_col_width <- prod(md[['n_cols']][(i+1):length(col_vars)])
      header[[i]] <- paste("\\multicolumn{",header_col_width,"}{c}{", levels(df[,col_vars[i]]) ,"}",sep="")
      lines[[i]] <- header_col_width
    }
  }
  
  
  header_lines <- list()
  for(i in length(col_vars):1){
    header_col_repeats <- ifelse(i==1,1,prod(md[['n_cols']][1:(i-1)]))
    n_levels <- md[['n_cols']][i]
    length <- lines[[i]]
    left_cmidrules <- rep(1,header_col_repeats*n_levels)
    right_cmidrules <- left_cmidrules
    r <- 0
    for(j in 1:header_col_repeats){
      for(k in 1:n_levels){
        r <- r+1
        if(r>1){ left_cmidrules[r] <- right_cmidrules[r-1] + 1 }
        right_cmidrules[r] <- left_cmidrules[r] + length - 1
      }
    }
    midrule <- paste("\\cmidrule(lr){",length(row_vars)+left_cmidrules,"-",length(row_vars)+right_cmidrules,"}",sep="")
    
    #PUT IN ROW VAR HEADER ON BOTTOM LINE, A SPACE IN ALL OTHER HEADER LINES
    if(i == length(col_vars)){
      row_var_headers <- colnames(df)[row_vars]
      #PUT MIDRULE LINES ON ROW VAR HEADERS
      midrule <- paste(c(paste("\\cmidrule(lr){",1:length(row_vars),"-",1:length(row_vars),"}",sep=""),midrule),collapse="")
    }
    else{
      row_var_headers <- rep(" ",length(row_vars))
    }
    
    
    header_lines[[i]] <- paste(paste(c(row_var_headers,rep(header[[i]],header_col_repeats)),collapse="&"),"\\\\",paste(midrule,collapse=""),"\n",sep="")
    
  }
  
  
  cat(
    "\\refstepcounter{table}\\label{",label,"}\\begin{tabularx}{",table_width ,"\\linewidth}{",
    paste(rep("l",length(row_vars)),collapse=""),
    paste(rep(cell_h_align,prod(md[['n_cols']])),collapse=""),
    '}\\multicolumn{',length(row_vars) + prod(md[['n_cols']]),'}{@{}l@{}}{\\bfseries Table \\ref{',label,'}. \\parbox[t]{',table_width,'\\linewidth -\\widthof{Table XXXX.}}{',title,'\\rule[-1ex]{0ex}{1ex}}}\\\\\\toprule\n',
    do.call(paste,header_lines),
    lt1,
  '\\multicolumn{',length(row_vars) + prod(md[['n_cols']]),'}{@{}l@{}}{\\parbox[t]{',table_width,'\\linewidth}{ \\rule{0pt}{1.5em}\\raggedright ',footnote, '}}\\\\',
    "\\end{tabularx}",
    sep="",
    file=file
  )
  
  return(out)
}



###########################################################################################################

rate_table <- function(data,binary,variables){
  
  out_table <- NULL
  for(i in 1:length(variables)){
    df <- data.frame(Y=data[,binary],X=data[,variables[i]],V=gsub("\\_"," ",variables[i]),level=" ",stringsAsFactors=FALSE,rate=0,lb=0,ub=1)
    rtl <- by(df,df$X,function(x){bt <- binom.test(sum(x$Y),nrow(x)); x[1,'rate'] <- bt$estimate; x[1,c("lb","ub")] <- c(bt$conf.int); x[1,'level'] <- as.character(x$X[1]); return(x[1,])  })
    rtl <- do.call(rbind,rtl)

    out_table <- rbind(out_table,rtl[rtl$level!="REMOVE",])
  }
  
  rownames(out_table) <- NULL;
  out_table$level <- gsub("^[0-9]+\\. ","",out_table$level)
  return(out_table[,c("V","level","rate","lb","ub")])
}



###########################################################################################################
plot_rate_table <- function(rt,plot_width,table_width,table_title,label,file="D:/erase/rate_table.txt"){
  
 out<- rt
 out$plot <- paste("\\begin{tikzpicture} 
                    \\draw [draw=black!50] (0in,0ex) -- (",plot_width,"in,0ex); 
                    \\draw [fill=black,draw=black] (",plot_width*rt$rate,"in,0ex) circle (.5ex); 
                    \\draw [draw=black,ultra thick] (",plot_width*rt$lb,"in,0ex) -- (",plot_width*rt$ub,"in,0ex);
                    \\end{tikzpicture}",sep="" )
 out$prate <- sprintf("%3.0f",100*out$rate)
 out$ci <- paste("(",round(100*out$lb,0),", ",round(100*out$ub,0),")",sep="")
 
 for(j in nrow(out):2){
   if(out$V[j] == out$V[j-1]){ out$V[j] <- " "}
   else{ out$plot[j-1] <- paste(out$plot[j-1],"\n\\rule{0pt}{0pt} \\\\",sep="") }
 
 }

axis <- paste(
"\\begin{tikzpicture}",
"\\draw [draw=black] (0in,.5ex) -- (0in,0ex) -- (",plot_width,"in,0ex) -- (",plot_width,"in,.5ex);",
"\\node (zero) at (0in,0ex) [anchor=north west] {\\scriptsize 0\\%};",
"\\node (zero) at (",plot_width,"in,0ex) [anchor=north east] {\\scriptsize 100\\%};",
"\\end{tikzpicture}",
sep="")

out <- rbind(out,out[1:2,])
out[nrow(out),"plot"] <- axis
out[nrow(out)-1,"plot"] <- " "
out[nrow(out)-c(0:1),c("V","level","prate","ci")] <- " "
 
cat(
"\\refstepcounter{table}\\label{",label,"}\n",

#"\\begin{tabularx}{", table_width,"\\linewidth}",
"\\begin{tabular}",

"{llrrl}",
"\\multicolumn{5}{@{}l@{}}{Figure \\ref{",label,"}. ",table_title,"}\\\\ \\toprule\n",
"Group & Level &Rate &CI \\\\ \\midrule\n",
 paste(apply(out[,c("V","level","prate","ci","plot")],1,paste,collapse="&"),
       collapse="\\\\")
,

#"\\end{tabularx}",
"\\end{tabular}",

sep="",
file=file
)
}


######

mypar2 <- function(){ 
 par(mar=c(1,1,.1,.1)*2,
  mgp=c(1,0,0),
  tcl=.2,
  cex.axis=.75,
  col.axis="black",
  pch=16)
}

bg.color<-function(color){
  ttt<-par()$usr
  rect(ttt[1],ttt[3],ttt[2],ttt[4],col=color)
  xxx<-axTicks(1)
  yyy<-axTicks(2)
  xxx<-c(xxx[1]-diff(xxx)[1]/2,xxx,xxx+diff(xxx)[1]/2)
  yyy<-c(yyy[1]-diff(yyy)[1]/2,yyy,yyy+diff(yyy)[1]/2)
  grid(col="white",lwd=1,lty=1)
  abline(col="white",
   lwd=.5,
   lty=1,
   v=xxx,
   h=yyy)
  }

#bg.color<-function(color){
#  ttt<-par()$usr
#  rect(ttt[1],ttt[3],ttt[2],ttt[4],col=color)
#  grid(col="white",lwd=1.2,lty=1)
#  }

sigma<-function(lm) summary(lm)$sigma 

TexLegend<-function(x="topleft",y=NULL,legend,env=e1,...){

	w<-cbind(lapply(strsplit(legend,character(0)),function(x){sum(x=="$")}))
	w<-as.numeric(w)
	e<-sapply(nchar(legend)-w-1,function(x){paste(rep("x",x,sep=""),collapse="")})
	psfragLegend<-paste(3:(length(legend)+2),e,sep="")
	legend(x,y,legend=psfragLegend,...)
	show(env)
	env$psfrag.labels<-c(env$psfrag.labels,psfragLegend)
	env$psfrag.real<-c(env$psfrag.real,legend)
}

TexPlot<-function(x,y,xlab="x",ylab="y",width=4,height=3,file="TexPlot",caption="",main="",MorePlots=NULL,...){

   postscript(file=paste(file,".eps",sep=""),width=width,height=height,horizontal=FALSE,family="Courier")

	par(mar=c(1,1,1,1))

#	Xnchar<-nchar(xlab)-sum(strsplit(xlab,split=character(0))[[1]]=="$")
#	Ynchar<-nchar(ylab)-sum(strsplit(ylab,split=character(0))[[1]]=="$")

	psfrag.labels<-NA

#	psfrag[1]<-paste(rep("x",Xnchar),sep="",collapse="")
#	psfrag[2]<-paste(rep("y",Ynchar),sep="",collapse="")

	psfrag.labels[1]<-"1"
	psfrag.labels[2]<-"2"

	plot(x,y,xaxt="n",yaxt="n",xlab="",ylab="",...)
	xaxis<-axis(1,lwd=0,line=5)
	yaxis<-axis(2,lwd=0,line=5)

#	xlabels<-paste("x",letters[1:length(xaxis)],sep="")
#	ylabels<-paste("y",letters[1:length(yaxis)],sep="")

	xlabels<-letters[1:length(xaxis)]
	ylabels<-LETTERS[1:length(yaxis)]

		xaddchar<-nchar(xaxis)-1
		xadd<-sapply(xaddchar,function(x){paste(rep("0",x),collapse="")})
		yaddchar<-nchar(yaxis)-1
		yadd<-sapply(yaddchar,function(x){paste(rep("0",x),collapse="")})
		xlabels<-paste(xlabels,xadd,sep="")
		ylabels<-paste(ylabels,yadd,sep="")

	mtext(psfrag.labels[1],line=0,side=1,adj=.02)
	mtext(psfrag.labels[2],line=0,side=2,adj=.02)
	xticks<-axis(3,tck=0.01)
	axis(1,tck=0.01)
	yticks<-axis(4,tck=0.01)
	axis(2,tck=0.01)
	axis(3,at=xaxis,labels=xlabels,line=-1,lty=0,cex.axis=1)
	axis(4,at=yaxis,labels=ylabels,line=-1,lty=0,cex.axis=1)

	psfrag.labels<-c(psfrag.labels,xlabels,ylabels)
	psfrag.real<-c(xlab,ylab,xaxis,yaxis)

	if(!is.null(MorePlots)){
		e1 <- environment()
		show(e1)
		show(psfrag.labels)
		write.table(MorePlots,file="temp.R",quote=FALSE,row.names=FALSE,col.names=FALSE,append=FALSE)
		sys.source("temp.R", envir=e1)
		show(psfrag.labels)
	}

	dev.off()

	psfrag<-paste("\\psfrag{",psfrag.labels,"}{",psfrag.real,"}",sep="")

	write.table(paste("\\begin{myfig}{",main,"}",sep=""),file=paste(file,".tex",sep=""),
		quote=FALSE,row.names=FALSE,col.names=FALSE,append=FALSE)
	write.table("\\begin{center}",file=paste(file,".tex",sep=""),
		quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
	write.table(psfrag,file=paste(file,".tex",sep=""),
		quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
	write.table( paste(c("\\includegraphics{",file,"}"),sep="",collapse=""),
		file=paste(file,".tex",sep=""),
		quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
	label<-paste(c("\\label{",file,"}"),sep="",collapse="")
	write.table(label,file=paste(file,".tex",sep=""),
		quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
	write.table("\\end{center}",file=paste(file,".tex",sep=""),
		quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
	if(caption!="") caption<-paste(caption,"\\\\",sep=" ")
	write.table(caption,file=paste(file,".tex",sep=""),
		quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
	write.table("\\end{myfig}",file=paste(file,".tex",sep=""),
		quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
}


TexTable<-function(tbl,
  file="TexTable",
  main="",
  hline.after=0,
  include.rownames=TRUE,
  digits=NULL,
  width=.8,
  column=FALSE,
  font.size="",
  ...){

	label<-paste(c("\\label{",file,"}"),sep="",collapse="")
	rule<-paste(c("\\rule{", width, "\\textwidth}{1pt} \\\\"),sep="",collapse="")
	if(column==TRUE) rule<-paste(c("\\rule{", width, "\\columnwidth}{1pt} \\\\"),sep="",collapse="")

	write.table(c("\\begin{center} \\refstepcounter{table}", label , rule,font.size),file=paste(file,".tex",sep=""),
		quote=FALSE,row.names=FALSE,col.names=FALSE,append=FALSE)
	
	xtbl<-xtable(tbl,digits=digits)
	ncols<-length(align(xtbl))
	if(include.rownames==FALSE){
			align(xtbl)[3:ncols]<-rep(">{\\raggedleft\\arraybackslash}X",ncols-2)
			align(xtbl)[2]<-paste(width,"\\textwidth}{l",sep="")
			if(column==TRUE) align(xtbl)[2]<-paste(width,"\\columnwidth}{l",sep="")
			align(xtbl)[ncols]<-paste(c(align(xtbl)[ncols],"} \n  \\multicolumn{",ncols-1,"}{l}{\\bf\\large Table \\arabic{table}: ",main ,"} \\\\ \\ \\\\ \n % {"),sep="",collapse="")  
	}else{
	align(xtbl)[2:ncols]<-rep(">{\\raggedleft\\arraybackslash}X",ncols-1)
	align(xtbl)[1]<-paste(width,"\\textwidth}{l",sep="")
	if(column==TRUE) align(xtbl)[1]<-paste(width,"\\columnwidth}{l",sep="")

	align(xtbl)[ncols]<-paste(c(align(xtbl)[ncols],"} \n  \\multicolumn{",ncols,"}{l}{\\bf\\large Table \\arabic{table}: ",main ,"} \\\\ \\ \\\\ \n % {"),sep="",collapse="") 
	} 
	print(xtbl,floating=FALSE,tabular.environment="tabularx",append=TRUE, file=paste(file,".tex",sep=""),hline.after=hline.after,include.rownames=include.rownames,...)

	write.table(c( rule, "\\end{center}"),file=paste(file,".tex",sep=""),
		quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
}

SummaryTable <- function(data,
  units=NULL,
  file="SummaryTable",
  main="Variable Summary",
  show.factors=TRUE,...){
        
        if(!is.null(units)) units<-paste("(",units,")",sep="")
        colnames(data)<-paste(colnames(data),units)
        
        Numeric.Flag<-sapply(data,is.numeric)
        Factor.Flag<-sapply(data,is.factor)

        out<-NULL
        SumTable1<-NULL

        if(sum(Numeric.Flag)>0){
                data1<-as.matrix(data[,Numeric.Flag])

                SumTable1<-cbind(
                        apply(data1,2,function(x){NA}),
                        apply(data1,2,function(x){r<-sum(is.na(x)) ; return(100*r/length(x))}),
                        apply(data1,2,min,na.rm=T),
                        apply(data1,2,mean,na.rm=T),
                        apply(data1,2,max,na.rm=T))
        }

        if(sum(Numeric.Flag)==1) rownames(SumTable1)<-colnames(data)[Numeric.Flag]

        if(sum(Factor.Flag)==1){
                data2<-data[,Factor.Flag]
                numlevels<-nlevels(data2)
                out<-matrix(NA,ncol=5,nrow=1+numlevels)
                rownames(out)<-rep("",nrow(out))
                rownames(out)[1]<-colnames(data)[Factor.Flag]
                r<-sum(is.na(data2))/length(data2)
                out[1,2]<-r

                for(j in 1:numlevels){
                        out[1+j,1]<-100*(table(data2)/length(data2))[j]
                        rownames(out)[1+j]<-paste("\\rule{2ex}{0pt}",levels(data2)[j])
                }
        }

        if(sum(Factor.Flag)>1){
                data2<-data[,Factor.Flag]
                SumTableF1<-apply(data.matrix(data2),2,table)
                if(class(SumTableF1)=="matrix") SumTableF1 <- as.list(as.data.frame(SumTableF1))
                numlevels<-sapply(data2,function(x){nlevels(x)})
                nfac<-ncol(data2)
                tot.numlevels<-sum(numlevels)

                out<-matrix(NA,nrow=nfac+tot.numlevels,ncol=5)
                rownames(out)<-rep("",nrow(out))

                k<-1
                for(i in 1:nfac){
                        r<-sum(is.na(data2[[i]]))
                        out[k,2]<- 100*r/length(data2[[i]])
                        rownames(out)[k]<-colnames(data2)[i]
                        k<-k+1
                        for(j in 1:numlevels[i]){
                                out[k,1]<-100*as.vector((SumTableF1)[[i]]/length(data2[[i]]))[j]
                                rownames(out)[k]<-paste("\\rule{2ex}{0pt}",levels(data2[[i]])[j])
                                k<-k+1
                        }
                        
                }
        }

        out2<-matrix(NA,nrow=2,ncol=5)
        Total.Label<-paste("Total (N=",nrow(data),")",sep="")
        rownames(out2)<-c("",Total.Label)
        out2[2,2]<-100*sum(!complete.cases(data))/nrow(data)
        
        
        SumTable2<-rbind(SumTable1,out,out2)
        colnames(SumTable2)<-c("PPP","NACheat","Min","Mean","Max")
        SumTable3<-SumTable2
        rownames(SumTable3)<-NULL
        tbl<-data.frame(Variables=rownames(SumTable2),SumTable3,stringsAsFactors=FALSE)

        if(show.factors==FALSE) tbl<-tbl[!substr(tbl[,1],1,5)=="\\rule",]
        if(tbl[nrow(tbl),3]==0){ tbl<-tbl[,-3] ; tbl[nrow(tbl),1]<-paste("Total (N=",nrow(data)," all complete)",sep="")}
        if(sum(!is.na(tbl[,2]))==0) tbl<-tbl[,-2]

        col.sanitize.cheat<-function(x){
                PPP.flag<-substr(x,1,3)=="PPP"
                x[PPP.flag]<-"\\%"
                NACheat.flag<-substr(x,1,7)=="NACheat"
                x[NACheat.flag]<-"\\%NA"
                return(x)
        }

        TexTable(tbl,file=file,main=main,include.rownames=FALSE,
                sanitize.text.function=function(x){x},
                sanitize.colnames.function=col.sanitize.cheat,...)
        return(tbl)
}

