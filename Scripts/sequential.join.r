library(tidyverse)
pkgconfig::set_config("dplyr::na_matches" = "never")
sequential_join <- function(df1, df2,
                            var.to.join.1,
                            var.to.join.2="",
                            var.to.join.3="",
                            confidence,
                            ... ){
  interesting.vars <- rlang::enquos(...)
  join.col <- rlang::ensym(var.to.join.1) # To use it with select
  join.quo <- rlang::enquo(var.to.join.1) # To use it as right hand side of mutate
  varname <- rlang::ensym(confidence)     # To use it as left hand side of mutate
  
  anti_join(df1, df2, by = var.to.join.1) -> left_over # what does not have a match
  
  
  df2 %>% 
    select(..., !!join.col) %>% # Keep only: key variable, interesting variables
    inner_join(df1,., by = var.to.join.1) %>% 
    mutate(!!varname := (!!join.quo)) -> step.1 # The confidence level of the join 
  
  ### Now we do the second step join.   
  
  if(var.to.join.2 !=""){                      # So it only works if there are more than 1 key variables
    join.col <- rlang::ensym(var.to.join.2)
    
    join.quo <- rlang::enquo(var.to.join.2)
    
    anti_join(left_over, df2, by = var.to.join.2) -> left_over2 # Whatever doesn't have a match on either 1 or 2 
    
    
    df2 %>% 
      select(..., !!join.col) %>% 
      
      distinct()  %>% # To stop duplication of column information - maybe change it later on to a group_by - mutate - case_when - both
      add_count(!!join.col) %>% 
      mutate_at(.vars = vars(...), .funs = function(x) ifelse (.$n == 1, x, "undetermined")) %>% 
      distinct() %>% 
      inner_join(left_over, .) %>%
      mutate(!!varname := (!!join.quo)) -> step.2
    
    
    left_over <- left_over2
    step.1 <- bind_rows(step.1, step.2)
    
    
  } 
  #}
  if(var.to.join.3!=""){  
    
    join.col <- rlang::ensym(var.to.join.3)
    
    join.quo <- rlang::enquo(var.to.join.3)
    
    anti_join(left_over, df2, by = var.to.join.3) -> left_over2
    
    df2 %>% 
      select(..., !!join.col) %>% 
      
      distinct()  %>% # To stop duplication of column information - maybe change it later on to a group_by - mutate - case_when - both
      add_count(!!join.col) %>% 
      mutate_at(.vars = vars(...), .funs = function(x) ifelse (.$n == 1, x, "undetermined")) %>% 
      distinct() %>% 
      inner_join(left_over, .) %>%
      mutate(!!varname := (!!join.quo)) -> step.2
    
    left_over <- left_over2
    step.1 <- bind_rows(step.1, step.2)
  }
  
  return(bind_rows(step.1, left_over))
}
