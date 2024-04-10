library(reshape2)
#' clean and combine VAX, AE and pat data from VAERS
#' @return cleaned dataframe
processData <- function(startyear, endyear,
                         case, control,
                         case_cnt=1, control_cnt=1){
  
  year=startyear:endyear
  
  print("loading data...")
  #### data_PAT
  dd_pat=NULL
  for (k in 1:length(year)){
    tmp=read.csv(paste('data/',year[k],"VAERSDATA",".csv",sep=""),na = "")
    dd_pat[[k]]=tmp[,c("VAERS_ID","RECVDATE", "AGE_YRS","CAGE_YR","CAGE_MO","VAX_DATE","NUMDAYS","SEX","STATE","DIED","DATEDIED","ALLERGIES")]
    
  }
  
  data_PAT=Reduce("rbind",dd_pat)
  data_PAT=data_PAT[order(data_PAT$VAERS_ID),]
  
  #-- clean the age variable (CAGE= CAGE_YR+CAGE_MO by vax_date-birthdate; NUMDAYS:ONSET_DATE-VAX_DATE)
  CAGE_MO=ifelse(is.na(data_PAT$CAGE_MO),0,data_PAT$CAGE_MO)
  CAGE=data_PAT$CAGE_YR+CAGE_MO
  data_PAT$AGE=ifelse(is.na(CAGE), data_PAT$AGE_YRS, CAGE)
  
  
  #### data_VAC
  dd_vac=NULL
  for (k in 1:length(year)){
    tmp=read.csv(paste('data/',year[k],"VAERSVAX",".csv",sep=""), na = "")
    dd_vac[[k]]=tmp[,c("VAERS_ID","VAX_NAME","VAX_TYPE")]
  }
  
  data_VAC=Reduce("rbind",dd_vac)
  data_VAC=data_VAC[order(data_VAC$VAERS_ID),]
  
  
  #### data_AE
  dd_ae=NULL
  for (k in 1:length(year)){
    tmp=read.csv(paste('data/', year[k],"VAERSSYMPTOMS",".csv",sep=""), na = "")
    tmp=tmp[,c("VAERS_ID","SYMPTOM1","SYMPTOM2","SYMPTOM3","SYMPTOM4","SYMPTOM5")]
    tmp[tmp==""] = NA
    tmp.long=melt(tmp, id.vars = c('VAERS_ID'), measure.vars = 2:6, na.rm = T, value.name = "AE_NAME")
    dd_ae[[k]]=tmp.long
    dd_ae[[k]]$year=year[k]
  }
  data_ae=Reduce("rbind",dd_ae)
  data_ae=data_ae[order(data_ae$VAERS_ID),]
  
  
  print("preprocess...")
  meddra_list = read.csv("data/final_MedDRA.csv",header=T)
  # merge with PT ID
  meddra_ptlist = meddra_list %>% distinct(PT, PT_ID) %>% transmute(AE_NAME = PT, MEDDRA_ID = PT_ID)
  # If cannot mapped PT, then merge with LLT and use the first PT ID
  meddra_lltlist = meddra_list %>% distinct(LLT, PT_ID) %>% transmute(AE_NAME = LLT, MEDDRA_ID = PT_ID)
  data_AE = merge(data_ae, meddra_ptlist, by = 'AE_NAME', all.x = T) %>% as_tibble()
  data_AE[is.na(data_AE$MEDDRA_ID), ] = merge(data_AE[is.na(data_AE$MEDDRA_ID), 1:4],
                                              meddra_ptlist, by = "AE_NAME", all.x = T) %>% as_tibble()
  #length(unique(data_AE$AE_NAME))  > 7788
  
  #### combine VAX-AE-PAT files
  dd_VAC_AE=merge(data_VAC,data_AE, by=c("VAERS_ID"))
  dd_VAC_AE_PAT=merge(dd_VAC_AE,data_PAT, by=c("VAERS_ID"),all.x = TRUE)
  dds=dd_VAC_AE_PAT[order(dd_VAC_AE_PAT$VAERS_ID),]
  
  dds = dds %>% mutate(VAX_YEAR = substr(VAX_DATE, nchar(VAX_DATE)-3, nchar(VAX_DATE)) %>% as.numeric(),
                       AE_YEAR = year, AE_ONSET_DAYS = NUMDAYS) %>%
    dplyr::select(VAERS_ID, AE_NAME, MEDDRA_ID, RECVDATE, 
                  AGE, SEX,
                  VAX_NAME, VAX_TYPE,  VAX_DATE
    )
  
  
  dds_sub = dds %>% filter( grepl( case, VAX_TYPE, fixed = TRUE) | grepl( control, VAX_TYPE, fixed = TRUE))%>%
    mutate(VAX_LABEL = ifelse(grepl( case, VAX_TYPE, fixed = TRUE), case, control )) 
   
  dd_case <- dds_sub %>% 
    filter(VAX_LABEL==case) %>% 
    group_by(AE_NAME) %>%
    filter(n() > case_cnt)
  
  dd_control <- dds_sub %>% 
    filter(VAX_LABEL==control) %>% 
    group_by(AE_NAME) %>%
    filter(n() > control_cnt)
  
  dds_sub = dds_sub %>%
    filter((AE_NAME %in% dd_case$AE_NAME) & (AE_NAME %in% dd_control$AE_NAME))

  return(dds_sub)

}