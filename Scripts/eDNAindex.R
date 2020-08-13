# Usage eDNAindex(x, sample, taxa, nReads) all elements without quotes


eDNAindex<- function(df, Sample_column, OTU_column, Counts_column, Biological.replicate ){ # 
  
  require(dplyr)
  require(rlang)
  
  
  Sample_column <- rlang::enquo(Sample_column)
  OTU_column    <- rlang::enquo(OTU_column)
  Counts_column <- rlang::enquo(Counts_column)
  Biological.replicate <- rlang::enquo(Biological.replicate)

  
  if (quo_name(Biological.replicate) %in% colnames(df)){
    
    print ("Averaging ratios between Biological replicates")
    
    df %>% 
      
      group_by(!!Sample_column, !!OTU_column, !!Biological.replicate) %>%
      
      summarise ( sumreads = sum(!!Counts_column))  %>%  # This sums technical replicates

      group_by(!!Sample_column,!!Biological.replicate) %>%

      mutate (Tot = sum(sumreads),
              Row.prop = sumreads / Tot)  %>%                      # This creates the proporttion on each biological replicate

      group_by(!!Sample_column) %>%

      mutate (nreps = length(unique(!!Biological.replicate))) %>%

      group_by(!!Sample_column, !!OTU_column) %>%

      summarise (mean.prop = sum (Row.prop) / max(nreps))   %>%

      group_by (!!OTU_column) %>%

      mutate (Colmax = max (mean.prop),
              Normalized.reads = mean.prop / Colmax) %>%

      dplyr::select( -Colmax, -mean.prop) -> output
    return(output)
  }


  print("Calculating eDNAindex directly")

  # IF THERE ARE TECHNICAL REPLICATES, WE NEED TO SUM THOSE VALUES FIRST
  #return(!!Technical_replicate)


  df %>%

    group_by(!!Sample_column, !!OTU_column) %>%

    summarise (sumreads = sum(!!Counts_column)) %>% # In case the sample column is of a higher group than the nrows

    group_by(!!Sample_column) %>%

    mutate (Tot = sum(sumreads),
            Row.prop = sumreads / Tot) %>%

    group_by (!!OTU_column) %>%

    mutate (Colmax = max (Row.prop),
            Normalized.reads = Row.prop / Colmax) %>%
    dplyr::select(-Tot, -Row.prop, -Colmax, -sumreads)  #note, specifying dplyr::select to avoid conflict w MASS package


  
}
