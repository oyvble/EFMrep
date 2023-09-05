#' @title efm2
#' @author Oyvind Bleka
#' @description efm2 is a GUI wrapper for the EFMrep package
#' @param envirfile A Rdata file including a saved environment of a project
#' @param envir A environment object
#' @export

#library(EFMrep);envir=NULL;envirfile=NULL
#efm2(envirfile)
efm2 = function(envirfile=NULL,envir=NULL) {
  require(euroformix)
  #size of main window
  mwH <- 400
  mwW <- 1000
  
  #type of gwidgets-kit
  #library(gWidgetstcltk)
  options(guiToolkit="tcltk")
  
  pack = "EFMrep" #name of package and program
  #version:
  version = packageVersion(pack) #follows same version as package number
  softname <- paste0(pack," v",version)
  
  #GUI-Restriction on maximum number of contributors
  maxKsetup <- 4 
  
  #Spacing between widgets
  spc <- 10
  .sep <- .Platform$file.sep # Platform dependent path separator. 
  pgkPath <- path.package(pack, quiet = FALSE) # Get package path.
  
  #Option files are storer in system (opt-settings): These are default values
  optList = list(
    optSetup =  list(AT=50,pC=0.05,lambda=0.01,fst=0), #default global settings
    optFreq =list(freqsize=0,minF=0,normalize=1), #option when new frequencies are found (size of imported database,minFreq), and missmatch options
    optMLE = list(nDone=3,delta=1,difftol=0.01,maxThreads=0,seed=0,steptol=1e-3,equaltol=0.01), #options when optimizing,validation (nDone,delta)
    optDC = list(alphaprob=0.99,maxlist=20), #options when doing deconvolution
    optFreqfile =  NULL #default frequency file 
  )  
  
  configList = list() #list of file name for each element

  #TRAVERSE EACH config FILE (configPrefix) AND CHECK IF CONFIG FILE XISTS
  for(optName in names(optList) ) {
    tmp <- paste0("config",gsub("opt","",optName))
    configFile <- configList[[optName]] <- paste0(pgkPath,.sep,tmp)
    if( file.exists(configFile) ) {  #use default values if not existing
      optF <- scan(file=configFile,what=character(),quiet=TRUE)
      for(i in 1:length(optF)) {
        suppressWarnings({
          tmp = as.numeric(optF[i])
        })
        if(is.na(tmp)) tmp = optF[i] #keep original
        if(length(tmp)==0) tmp = NULL #no value
        optList[[optName]][[i]] = tmp #convert direclty from file
      } 
      if(length(optF)==1) optList[[optName]] = optList[[optName]][[1]]
    }
      
  }

  
  #####################
  #create environment #
  #####################

  if(!is.null(envir)) {
    mmTK = envir #use environment direclty
  } else if(!is.null(envirfile)) {
    load(envirfile) #loading environment
    
    #NEED TO ADJUST loglik tolerance here
    opt = get("optMLE",envir=mmTK) #not used in program (yet)
    if(!is.null(opt$maxIter)) {
      idx = which(names(opt)=="maxIter")
      opt[idx] = 0.01
      names(opt)[idx] = "difftol"
    }
    assign("optMLE",opt,envir=mmTK) #store again
  } else {
    mmTK = new.env( parent = globalenv() ) #create new envornment object (must be empty)
  
    #Toolbar options: can be changed any time by using toolbar
    for(optName in names(optList) ) assign(optName,optList[[optName]],envir=mmTK)  #store to envir
    assign("optRepMarkerSetup",NULL,envir=mmTK) #not used in program (yet)
    
    #initializing environment variables
    assign("workdir",NULL,envir=mmTK) #store work dir
    assign("version",NULL,envir=mmTK) #Store version
    
    #Check frequency file and try import it.:    
    popFreq = NULL #try load frequency data from system file
    if(!is.null( optList$optFreqfile ) ) {
      tryCatch({
        popFreq = euroformix::freqImport( optList$optFreqfile ,url=FALSE,xml=FALSE)[[1]] 
      }, error=function(e) print("System frequency file not valid."))
    }
    assign("popFreq",popFreq,envir=mmTK) 
    
    #imported data:
    assign("mixData",NULL,envir=mmTK) 
    assign("refData",NULL,envir=mmTK) 
    assign("repSettings",NULL,envir=mmTK)  #this is settings for each sample (AT,dropin etc)
    assign("modelSettings",NULL,envir=mmTK)  #this is model settings (how params are relataed)
    assign("relSettings",NULL,envir=mmTK)  #this is settings for relatedness (list of hp and hd separately). Object is modfied in function "kinshipSelectorGUI"
        
    #models: stored setups for model specification
    assign("calcList",NULL,envir=mmTK) #store model + calculated results  

    #results: stored results after calculations
    #assign("resEVID",NULL,envir=mmTK) #assign evidence weighting results (i.e. calculated LR with MLE estimates)
    assign("resDC",NULL,envir=mmTK) #assign deconvolved results (i.e. ranked tables of results)
    assign("resCompare",NULL,envir=mmTK) #assign information when doing model selection searchs (stored tables)
  }

  ####################################
  #auxiliary functions and variables:#
  ####################################
  
  #helpfunction to get small number from log-value
  getSmallNumber = function(logval,sig0=2,scientific="e") {
    log10base = logval/log(10) #convert to 10 base
    power = floor(log10base) #get power number
    remainder = log10base - power
    return( paste0( round(10^remainder,sig0),scientific,power)) #representation of very small numbers (avoid underflow)
  }
  
   NAtoSign <- function(x) {
    x[is.na(x)] <- "-" #NB: New version of gtable does not accept NA values
    return(x)
   }
  # helptext = function(obj,txt) { gWidgets2::addHandlerRightclick(obj,handler=function(h,...) { gWidgets2::gmessage(txt,title="Detailed information") }) }
  helptext = function(obj,txt) { gWidgets2::tooltip(obj) = txt } 
  
  #helpfunction to return minimum frequency (used for new alleles)
  getminFreq = function() {
    popFreq <- get("popFreq",envir=mmTK) #get selected popFreq
    freqsize <- get("optFreq",envir=mmTK)$freqsize #get selected size of frequence-database
    minfreq <- get("optFreq",envir=mmTK)$minF #get selected size of frequence-database
    if(freqsize>0) {
     return(5/(2*freqsize))
    } else if(is.null(minfreq) || minfreq==0 ) {
     return(min(unlist(popFreq))) #minumum observed frequency was used 
    } else {
     return(minfreq) #specified observed frequency was used 
    }
  }
  
  #Function to get data from environment
  #sel used to select a specific datasubset
  getData = function(type,sel=NULL) {
   Data <- NULL
   if(type=="mix") Data <- get("mixData",envir=mmTK) #assign kit to mmTK-environment
   if(type=="ref") Data <- get("refData",envir=mmTK) #assign kit to mmTK-environment 
   if(!is.null(sel)) return(Data[sel]) #returns only selected datasubset
   return(Data)
  }
  
  #function for inserting sample/ref/db-names into existing gWidgets2::gcheckboxgroup
  getDataNames_type = function(type) {
    subD <- getData(type)
    if(!is.null(subD)) { return( names(subD))
    } else { return("") }
  }
  
  #Function which takes rownames and adds to first column
  addRownameTable = function(tab) {
    tmp <- colnames(tab)
    tab <- cbind(rownames(tab),tab)
    colnames(tab) <- c("X",tmp)
    return(tab)
  }
  
  #save result table to file:
  saveTable = function(tab,sep="txt") {
    tabfile  = mygfile(text="Save table",type="save") #csv is correct format!
    if(length(tabfile)==0) return()
     if(length(unlist(strsplit(tabfile,"\\.")))==1) tabfile = paste0(tabfile,".",sep)
     if(sep=="txt" | sep=="tab") write.table(tab,file=tabfile,quote=FALSE,sep="\t",row.names=FALSE) 
     if(sep=="csv") write.table(tab,file=tabfile,quote=FALSE,sep=";",row.names=FALSE) 
     print(paste("Table saved in ",tabfile,sep=""))
  } #end file
  
  strsplit2 <- function(x,spl) {
    if(nchar(x)==0) return("")
    txt <- x
    for(j in 1:length(spl)) {
     txt <- unlist(strsplit(txt,split=spl[j]))
    }
    return(txt)
  }
  
  #Helpfunction to print tables (R v3.5.0) has problem showing tables in the GUI.
  printTable = function(x) {
    print(cbind(rownames(x),x),row.names=FALSE)
  }

  
###################################################################
#############GUI HELPFUNCTIONS#####################################
###################################################################

 #TAKEN FROM CASESOLVER: This function is written since the encoding in  gWidgets2::gfile is fixed to UTF-8 which doesn't handle special letters
 mygfile <- function(text,type,filter=list(),initf=NULL) { #Possible bug: text ignored when type="selectdir"
   file <- gWidgets2::gfile(text=text,type=type,filter=filter,initial.filename=initf)
   Encoding(file) <- options()$encoding #Set to local encoder: Handle special cases.
   return(file)
 }
 
 #Helpfunction to get focus 
 getFocus = function() {
   gWidgets2::visible(mainwin) <- TRUE
   gWidgets2::focus(mainwin) <- TRUE
 }
 
 #Menu bar file-lists:
 f_setwd = function(h,...) {
  dirfile = mygfile(text="Select folder",type="selectdir")
  if(length(dirfile)==0) return()
  setwd(dirfile)
  assign("workdir",dirfile,envir=mmTK) #assign working directory
 }
  
 f_openproj = function(h,...) {
  projfile = mygfile(text="Open project",type="open", filter=list("Project"=list(patterns=list("*.Rdata")) , "All files"=list(patterns=list("*"))))
  if(length(projfile)==0) return()
  gWidgets2::dispose(mainwin)
  efm2(projfile) #send environment into program
 }
 
 f_saveproj = function(h,...) {
  projfile = mygfile(text="Save project",type="save")
  if(length(projfile)==0) return()
  
   if(length(unlist(strsplit(projfile,"\\.")))==1) projfile = paste0(projfile,".Rdata") #add extension if missing
   tmp = sort(sapply(mmTK,object.size)/1e6,decreasing=TRUE)
   #print(paste0("Size of stored objects (in MB): ",round(sum(tmp),2))) #prints size of each stored object
   #print(tmp) #prints size of each stored object
   selectDataToModel(h=list(action="SAVE")) #store data before saving project
   #save(mmTK,file=projfile,compress="xz") #save environment, #Dont compress!
   save(mmTK,file=projfile,compress="xz",eval.promises=FALSE,precheck=FALSE,compression_level=2)
   print(paste("Project saved in ",projfile,sep=""))
  
 }

 
 f_settings = function(h,...) { #creates a table, with different settings
   isGLOBAL = h$action=="GLOBAL" #whether settings are aglobal or not
   
   title = "Settings"
   opt <- get("optSetup",envir=mmTK)  #get global (or rep) param settings (always as default)
   vars = names(opt) #c("AT","pC","lambda") 
   opt2 = opt #make a copy (use first element as default)

   if(!isGLOBAL) { #if evid specific settings
     evidName = h$action #obtain evid name
     title = paste0(title," for ",evidName)
     for(var in vars) {
       if(!is.na(opt2[[var]][evidName])) opt2[[var]][1] = opt[[var]][evidName][1] #set first
     }
   }
   
   setwin <- gWidgets2::gwindow(title,visible=FALSE)
   tabval = gWidgets2::glayout(spacing=0,container=(setwin)) 
   w0 <- 15 #width of text box
   
   #Model parameters: 
   tabval[1,1] <- gWidgets2::glabel(text="Analytical threshold (AT)",container=tabval)
   tabval[1,2] <- gWidgets2::gedit(text=opt2$AT[1],width=w0,container=tabval)
   tabval[2,1] <- gWidgets2::glabel(text="Probability of drop-in (PrC)",container=tabval)
   tabval[2,2] <- gWidgets2::gedit(text=opt2$pC[1],width=w0,container=tabval)
   tabval[3,1] <- gWidgets2::glabel(text="Drop-in hyperparam (lambda)",container=tabval)
   tabval[3,2] <- gWidgets2::gedit(text=opt2$lambda[1],width=w0,container=tabval)
   
   if(isGLOBAL) { #set fst if global
     tabval[4,1] <- gWidgets2::glabel(text="Fst-correction (theta)",container=tabval)
     tabval[4,2] <- gWidgets2::gedit(text=opt2$fst[1],width=w0,container=tabval)
   }
   
   tabval[5,1] <- gWidgets2::gbutton("Save", container=tabval,handler = function(h, ...) { 
     opt2$AT[1] <- as.numeric(gWidgets2::svalue(tabval[1,2]))
     opt2$pC[1] <- as.numeric(gWidgets2::svalue(tabval[2,2])) 
     opt2$lambda[1] <- as.numeric(gWidgets2::svalue(tabval[3,2]))
     if(isGLOBAL) opt2$fst[1] <- as.numeric(gWidgets2::svalue(tabval[4,2])) 

     if( any( sapply(opt2,function(x) is.na(x[1])) ) ) { 
       gWidgets2::gmessage("Invalid input in settings, please set another value!",title="Wrong input",icon="error")
       return()
     }
     
     if(isGLOBAL) {     
       saveOpt =  c(opt2$AT[1],opt2$pC[1],opt2$lambda[1],opt2$fst[1]) #obtain relevant vars to save
       write(saveOpt,file=configList$optSetup)    #Save to file only for global settings
       opt = opt2 #must make a copy to keep in environment
     } else { #replicate specific
       for(var in vars) {
         opt[[var]][evidName] = opt2[[var]][1] #refresh with selected value 
       }
     }
     assign("optSetup",opt,envir=mmTK)  #assign user-value to opt-list
     
     gWidgets2::dispose(setwin)
   })
   gWidgets2::visible(setwin) <- TRUE
   
 }

 f_quitproj = function(h,...) {
  ubool <- gWidgets2::gconfirm("Do you want to save project?",title="Quit Program",icon="info")
  if(ubool) {
    f_saveproj(h)
  } else { 
   print("Program terminated without saving")
  }
  gWidgets2::dispose(mainwin) #remove window!
 }

 #helpfunction to get value in from user and store
 setValueUser <- function(what1,what2,txt,allowNULL=FALSE,allowText=FALSE) {
   listopt <- get(what1,envir=mmTK) #get object what 1.
   val <- listopt[[what2]]
   if(is.null(val)) val ="" #gwidgets2 does not handle NULL, must use empty string
   sw <- gWidgets2::gwindow(title="User input",visible=FALSE, width=300,height=50)
   grid <- gWidgets2::glayout(spacing=0,container=sw )
   grid[1,1] <- gWidgets2::glabel(txt, container=grid)
   grid[1,2] <- gWidgets2::gedit(text=val,container=grid,width=30)
   grid[2,1] <- gWidgets2::gbutton("OK", container=grid,handler = function(h, ...) { 
    GUIval = gWidgets2::svalue(grid[1,2]) #obtain GUI value
    if(allowNULL && GUIval=="") { #if accepting empty string
      tmp = NULL #Insert NULL
    } else {
      tmp <- GUIval
      if(!allowText) {
        tmp <- as.numeric(GUIval) #insert new value
        if(is.na(tmp)) {
          NAerror(what2)
          return()
        }
      }
    }
    listopt[[what2]] <- tmp
    assign(what1,listopt,envir=mmTK) #assign user-value to opt-list
    write(unlist(listopt),file=configList[[what1]]) #store selected values to file   
    gWidgets2::dispose(sw)
   } )
   grid[2,2] <- gWidgets2::gbutton("Cancel", container=grid,handler = function(h, ...) { gWidgets2::dispose(sw) } )
   gWidgets2::visible(sw) <- TRUE
 }

 #helpfunction to get value from user
 getValueUser <- function(txt="",val=0) {
  val2 <- gWidgets2::ginput(txt, text=val, title="User input",icon="question")
  return(val2)   
 }

 NAerror <- function(what) {
  gWidgets2::gmessage(paste0(what," must be specified as a valid value"),title="Wrong input",icon="error")
  #stop("Wrong user-input")
 }

 #helpfunction which checks that at value is in interval of [0,1]
 checkProb = function(x,what) {
  if(is.na(x)) NAerror(what)
  if(x < 0 || x>1) {
   gWidgets2::gmessage(paste0(what," must be specified in interval [0,1] "),title="Wrong input",icon="error")
   stop("Wrong user-input")
  }
 }
 checkPositive = function(x,what,strict=FALSE) {
  if(is.na(x)) NAerror(what)
  if(x < 0 ) {
   gWidgets2::gmessage(paste0(what," cannot be a negative number"),title="Wrong input",icon="error")
   stop("Wrong user-input")
  }
  if(strict && x==0) {
   gWidgets2::gmessage(paste0(what," cannot be zero"),title="Wrong input",icon="error")
   stop("Wrong user-input")
  }
 }
 checkPosInteger = function(x,what) {
  if(is.na(x)) NAerror(what)
  if(x < 1 || round(x)!=x) {
   gWidgets2::gmessage(paste0(what," must be a positive integer"),title="Wrong input",icon="error")
   stop("Wrong user-input")
  }
 }

 #helpfunction for printing evidence sample to terminal
 printEvid = function(subD) {
  locs <- names(subD) #get unique loci
  mixtab <- matrix(ncol=2,nrow=length(locs))
  for(loc in  locs) { #for each locus
        mixA <- subD[[loc]]$adata
        mixH <- subD[[loc]]$hdata
        if(!is.null(mixA)) mixtab[which(loc==locs),1] <- paste0(mixA ,collapse="/")
        if(!is.null(mixH)) mixtab[which(loc==locs),2] <- paste0(mixH ,collapse="/")
  }
  rownames(mixtab) <- locs
  colnames(mixtab) <- c("Allele","Height")
  printTable(mixtab)
 }  
 
 #helpfunction for printing reference sample to terminal
 printRefs = function(refD,refSel=NULL) {
   nR <- length(refSel) #number of selected references
   locs <- unique(unlist(lapply(refD,names))) #get unique loci
   reftab <- matrix(ncol=nR,nrow=length(locs)) #last row is RMP
   for(rsel in refSel) {
    for(loc in  locs) { #for each locus
      refA <-refD[[rsel]][[loc]]$adata
      if(!is.null(refA)) {
       reftab[which(loc==locs),which(rsel==refSel)] <- paste0(refA ,collapse="/")
      }
    }
   }
   rownames(reftab) <- locs
   colnames(reftab) <- refSel 
   printTable(reftab)
  }

  showTable = function(df,header="", dec=2, printTable=TRUE) {
    df = cbind(X=rownames(df),round(df,dec))
    if(printTable) print(df)
    dfwin <- gWidgets2::gwindow(header, visible=FALSE)#,height=mwH)
    tab <- gWidgets2::gdf(df,container=dfwin) #create table
    gWidgets2::visible(dfwin) <- TRUE
    gWidgets2::focus(dfwin)
    
  }
 
########################################################################################################################################## 
  
##################################### 
###########GUI WINDOW STARTS#########
##################################### 
 
 ##########
 #Menu bar#
 ##########
 mblst = list( #project saving and so on
  File=list(  
    gWidgets2::gaction('Set directory',handler=f_setwd),
    gWidgets2::gaction('Open project',handler=f_openproj),
    gWidgets2::gaction('Save project',handler=f_saveproj),
    gWidgets2::gaction('Settings',handler=f_settings, action="GLOBAL"),
    gWidgets2::gaction('Quit',handler=f_quitproj,icon="close")
  ),
  Frequencies=list(
    gWidgets2::gaction('Set size of frequency database',handler=function(h,...) {  
      setValueUser(what1="optFreq",what2="freqsize",txt="Set size of imported freq database \n(Min observed used if not spesified):") 
    }),
    gWidgets2::gaction('Set minimum frequency',handler=function(h,...) {  
      setValueUser(what1="optFreq",what2="minF",txt="Set minimum freq for new alleles\n(Min observed used if not spesified):") 
    }),
    gWidgets2::gaction('Set whether to normalize frequencies',handler=function(h,...) {  
      setValueUser(what1="optFreq",what2="normalize",txt="Should frequencies always add up to one\nafter including rare alleles? (1=YES,0=NO)") 
    })
  ),
  Optimization=list(
    gWidgets2::gaction('Set number of optimizations',handler=function(h,...) {  
      setValueUser(what1="optMLE",what2="nDone",txt="Set required number of (identical) optimizations:") 
    }),
    gWidgets2::gaction('Set variation of randomizer',handler=function(h,...) {  
      setValueUser(what1="optMLE",what2="delta",txt="Set variance of start point randomizer:") 
    }),
    gWidgets2::gaction('Set logLik tolerance',handler=function(h,...) {  
      setValueUser(what1="optMLE",what2="difftol",txt="Set number for tolerance:") 
    }),
    gWidgets2::gaction('Set maximum threads for computation',handler=function(h,...) {  
      setValueUser(what1="optMLE",what2="maxThreads",txt="Set max number of threads to be used in parallelisation:") 
    }),
    gWidgets2::gaction('Set seed of randomizer',handler=function(h,...) { 
      setValueUser(what1="optMLE",what2="seed",txt="Set seed of randomizer\n(Not used if zero):",allowNULL=TRUE) 
   }),
   gWidgets2::gaction('Set accuracy of optimization',handler=function(h,...) { 
     setValueUser(what1="optMLE",what2="steptol",txt="Set accuracy of optimization (steptol, see ?nlm):") 
   })
  ),
  Deconvolution=list(
    gWidgets2::gaction('Set required summed probability',handler=function(h,...) {  
      setValueUser(what1="optDC",what2="alphaprob",txt="Set required summed posterior genotype-probability of list:") 
    }),
    gWidgets2::gaction('Set max listsize',handler=function(h,...) {  
      setValueUser(what1="optDC",what2="maxlist",txt="Set size of maximum elements in deconvoluted list:") 
    })
  )
 )

##################################################################################################
########### Program starts #######################################################################
##################################################################################################

 #change working directory to the one stored in mmTK-environment
 wd=get("workdir",envir=mmTK) #assign working directory to mmTK-environment
 if(!is.null(wd) && dir.exists(wd)) setwd(wd)
 
 #Main window:
 mainwin <- gWidgets2::gwindow(softname, visible=FALSE, width=mwW,height=mwH)
 gWidgets2::gmenu(mblst,container=mainwin)
 nb = gWidgets2::gnotebook(container=mainwin)
 tabimport0 = gWidgets2::ggroup(horizontal=TRUE,spacing=10,container=nb,label="Data") #tab2: (imports all files)
 tabmodel = gWidgets2::glayout(spacing=spc,container=nb,label="Model") #tab3: specify model used in weight-of-evidence (INT/MLE) or in a Database search 
 tabMLE = gWidgets2::glayout(spacing=spc,container=nb,label="Results")#,expand=T,fill=T) #results from MLE
 tabDC = gWidgets2::ggroup(horizontal=FALSE,spacing=spc,container=nb,label="Deconvolution") #results from a deconvolution


####################################################
###############Tab 1: Import Data:##################
####################################################

 #When program starts, import assumed model for EVID.

  editboxsize = 3 #length of boxes 
 #b) load/save profiles/database: Supports any filetype
  f_importprof = function(h,...) {
    type=h$action #get type of profile
# type = "mix"
    #  proffile = mygfile(text=paste("Open ",type,"-file",sep=""),type="open",filter=list("text"=list(patterns=list("*.txt","*.csv","*.tab"))))
    proffile = mygfile(text=paste("Open ",type,"-file",sep=""),type="open",filter=list("text"=list(patterns=list("*.txt","*.csv","*.tab")),"all"=list(patterns=list("*"))))
    if(length(proffile)==0) return() 
    
    Data = euroformix::tableReader(proffile) #load profile
    Data = euroformix::sample_tableToList(Data) #convert from table to list 
  
    if(type=="ref") { #check that homozygote alleles are given twice
      txt2 <- "The number of alleles in a genotype must be 2."
      txt1 <- paste0("Only one allele was given for a genotype in a reference profile. ",txt2, " Hence the second allele was automatically set as the first allele.") 
      txt3 <- paste0("Too many alleles where given for a genotype in a reference profile. ",txt2) 
      miss <- FALSE #indicator whether hom. are missing
     for(kn in names(Data)) { #for each profile
      for(loc in names(Data[[kn]])) { #for each profile
         if( length(Data[[kn]][[loc]]$adata)>2 )  {
           gWidgets2::gmessage(txt3,"Wrong file-input",icon="error")
           break #breaking loop if wrong input
         }
         if( length(Data[[kn]][[loc]]$adata)==1) {
          Data[[kn]][[loc]]$adata <- rep(Data[[kn]][[loc]]$adata,2) #duplicated
          miss <- TRUE
         }
          Data[[kn]][[loc]]$hdata = NULL #Updated v2.1.0: Remove peak heights if these have been added for references.
      }
     }
     if(miss) gWidgets2::gmessage(txt1,"Warning",icon="info")
    }
  
    #get already stored data:
    if(type=="mix") Data2 <- getData("mix") #get data from mmTK-environment
    if(type=="ref") Data2 <- getData("ref") #get data from mmTK-environment
  
    oldNames = names(Data2) #old samples
    if(is.null(Data2)) { #if no previous already there
      Data2 <- Data
    } else {
      for(kn in names(Data)) Data2[[kn]] <- Data[[kn]] #insert dataframe for each profile
    }
    if(type=="mix")  assign("mixData",Data2,envir=mmTK) #assign data to mmTK-environment
    if(type=="ref")  assign("refData",Data2,envir=mmTK) #assign data to mmTK-environment
    sampleNames = names(Data2) #Updated samplenames
    
    #Update table
    if(type=="ref")  tabimportB[1,1][] <- cbind(sampleNames)
    if(type=="mix") {
      indIns = which(!sampleNames%in%oldNames) #indices to extend (not previously stored)
      for(e in indIns) {
        evidName = sampleNames[e] #obtain evidence name
        tabimportC[e,1] = gWidgets2::gcombobox(items=(0:e),container=tabimportC, selected=e+1, editable = TRUE)
        gWidgets2::size(tabimportC[e,1]) <- editboxsize[1]
        tabimportC[e,2] = gWidgets2::gcheckbox(evidName,container=tabimportC,checked = TRUE)
        gWidgets2::size(tabimportC[e,2]) <- max(nchar(sampleNames)) #editboxsize[2]
        tabimportC[e,3] <- gWidgets2::gcombobox(items=initKits, width= max(nchar(initKits)), selected = 0, editable = FALSE, container = tabimportC)
        tabimportC[e,4] <- gWidgets2::gbutton(text="Settings",container=tabimportC, handler=f_settings, action=evidName)
      }
    } 
  }


  #Compare references against each evidence
 f_compare = function(h,...) {
   evidD = getData("mix") #get selected references
   refD <- getData("ref")
   
   nEvids = length(evidD)
   nRefs = length(refD)
   if(nEvids==0 || nRefs==0) print("Missing data to compare")
   
   nLocs <- MAC <- matrix(0,nrow=nEvids,nRefs,dimnames = list(names(evidD),names(refD)))
   for(evid in names(evidD) ) { #for each selected evidence 
     subD <- evidD[[evid]] #selected samples
     locs <- names(subD) #obtain locs

     for(loc in locs) { #for each locus
        evidA = subD[[loc]]$adata #allele vector of evid
        if( length(evidA)==0) next #don't count missing markers
        
        for(ref in names(refD) ) { #for each selected evidence 
          refA <- refD[[ref]][[loc]]$adata
          if(is.null(refA) || length(refA)==0) next #skip if no data
          
          nLocs[evid,ref] = nLocs[evid,ref] + 1 #count locus
          MAC[evid,ref] <- MAC[evid,ref] + sum(refA%in%evidA) #count number of alleles
       }
     } 
   }
   matchrate <- MAC/(2*nLocs)
   nMissMatches = 2*nLocs - MAC
   print("---NUMBER OF MISSMATCHES FOR EACH REFERENCES/REPLICATES:")     
   printTable(nMissMatches)
 }
     
 #prints evidence, references, EPG, databases and population frequencies
 f_viewdata = function(h,...) {
  sep = "--------------------------------------------"
  help_gmessage = function(what) gWidgets2::gmessage(paste0("Please import and select ",what),icon="info")
  
  #Frequencies:
  print(sep)
  print("Frequencies") #get frequencies
  print(sep)
  print(get("popFreq",envir=mmTK)) #get frequencies
  
  #View evidence:
  evidD = getData("mix") #get selected references
  grpID = getRepGrpID( names(evidD) )
  evidGrp = grpID$evid
  kitGrp = grpID$kit
  
  for(evidName in names(evidGrp) ) {
      subD <- evidD[[evidName]] #selected samples
      print("------------------------------------")
      print(paste("Samplename: ",evidName,sep=""))
      printEvid(subD)
  }
  
  #View references:
  refD <- getData("ref")
  refSel <- numeric()
  if(!is.null(refD))  refSel <- gWidgets2::svalue(tabimportB[1,1])  #get selected references
  print("------------------------------------")
  print("References:")
  if(length(refSel)>0) printRefs(refD,refSel)
  
  #CREATE replicate EPG plots WITH SELECTED REPLICATES:
  if(!is.null(evidGrp))  {
    refData2 = NULL
    if(length(refSel)>0) refData2 = refD[refSel] #obtain selected reference data
    grps = unique(evidGrp) #obtain group indices
    for(grp in grps) {
      repEvids = names(evidGrp)[grp==evidGrp] #obtain replicates
      kit0 = unique(names(kitGrp)[grp==kitGrp])
      if(length(kit0)>1) {
        print("Several kit was selected for same group. Plot not provided.")
        next
      }
      if(length(kit0)==0 || is.na(kit0) || kit0=="NONE") { #if no kit provided
        euroformix::plotMPS2(evidD[repEvids],refData = refData2)
      } else { 
        euroformix::plotEPG2(evidD[repEvids],kit = kit0, refData = refData2)
      }
    }
  }
 }  #end viewdata
 
 f_fillsettings = function(h,...) { #helpfunction to fill out settings (kit or params) based on numbering
   # print("FILLING")
   evidNames = getDataNames_type("mix")
   if(evidNames[1]=="") return() #no data
   
   #Obtain grouped replicates (overview)
   grpID = getRepGrpID( evidNames )
   evidGrp <- grpID$evid
   kitGrp = grpID$kit
   
   #Fill in kit or settings info:
   opt = get("optSetup",envir=mmTK)  #Get rep settings from envir (includes default values)
   grps = unique(evidGrp) #obtain group IDs
   for(grp in grps) { #for each unique group
# grp=1
     evids =  names(evidGrp)[evidGrp==grp]
     indEvid = which(evidNames%in%evids) #get index of evidence
     
     if(h$action=="KITS") { #fill in kit if found at least one
       kit0 = na.omit(names(kitGrp)[kitGrp==grp]) #extract selected kit(s)
       if(length(kit0)>0) { #should not be empty
         for(e in indEvid)   gWidgets2::svalue(tabimportC[e,3]) = kit0[1] #insert
       } 
     }
     if(h$action=="SETTINGS") {
       for(var in names(opt) ) {
         optEvids = names(opt[[var]]) #evid stored in settings
         indVec = which(optEvids%in%evids) #obtain already registered
         if(length(indVec)>0) {
           insVal = opt[[var]][optEvids[indVec[1]]] #value to insert
           opt[[var]][evids] = insVal #refresh part of vector
         }
       }
     }
   } #end for each grpps
   if(h$action=="SETTINGS") {
     print("Settings were filled out!")
     assign("optSetup",opt,envir=mmTK) #assign values
   }
 } #end helpfunction
 
 getRepGrpID = function(evidNames) { #helpfunction to obtain GROUP ID for each sample (including kit)
   evidGrp <- kitGrp <- NULL #obtain grouped info about samples
   for(e in seq_along(evidNames)) { #traverse each rows and obtain group id
     if(gWidgets2::svalue(tabimportC[e,2])) { #if selected
       evidName = evidNames[e] #tabimportC[e,2][] #get name (splits name with ws separator). Must be avoided
       grp0 = as.integer( gWidgets2::svalue(tabimportC[e,1]) ) #obtain group the sample belongs to
       if(grp0==0) next #skip if selected as group 0 (ensure similar length) 
       evidGrp = append(evidGrp, setNames(grp0, evidName)  )
       kit0 = gWidgets2::svalue(tabimportC[e,3]) #obtain selected kit
       if(length(kit0)==0) kit0=NA
       kitGrp = append(kitGrp, setNames(grp0, kit0)  )
     }
   }
   return(list(evid=evidGrp,kit=kitGrp))
 }

 ###############
 #start layout:#
 ###############
 tabimport = gWidgets2::ggroup(container=tabimport0,horizontal = FALSE)
 tabimport2 = gWidgets2::ggroup(container=tabimport0)
 
 tabimportA = gWidgets2::glayout(spacing=3,container=gWidgets2::gframe("Step 1) Import Data",container=tabimport)) #kit and population selecter
 tabimportC = gWidgets2::glayout(spacing=1,container=gWidgets2::gframe("Step 2) Select Evidence(s) and Kit",container=tabimport2,expand=T,fill=T),expand=T,fill=T)
 tabimportB = gWidgets2::glayout(spacing=3,container=gWidgets2::gframe("Step 3) Select Reference(s)",container=tabimport,expand=T,fill=T),expand=T,fill=T) #evidence,ref dataframe
 tabimportD = gWidgets2::glayout(spacing=3,container=gWidgets2::gframe("Step 4) Select Interpretation",container=tabimport)) #Tasks button

 #Choose box and import button
 setAsImported = function(fn) { #set as imported
  gWidgets2::svalue(tabimportA[1,2]) <- "imported"
  gWidgets2::font(tabimportA[1,2]) <- list(weight="bold",size=11,color="green")
  helptext(tabimportA[1,2],fn) #set file name in helptext
 }
 
 tabimportA[1,1] = gWidgets2::gbutton(text="Import frequencies",container=tabimportA,handler=
  function(h,...) {
    fn = mygfile(text="Select file",type="open",filter = list(`All files` = list(patterns = c("*"))))
    if(length(fn)==0) return()
    popFreq = euroformix::freqImport(fn,url=FALSE,xml=FALSE)[[1]]
    assign("popFreq",popFreq,envir=mmTK) #assign popFreq
    
    #STORE FILE NAME
    setAsImported(fn) #set label as imported
    write(fn,file= configList$optFreqfile )   #STORE TO SYSTEM FILE
    optList$optFreqfile <<- fn
    assign("optFreqfile",fn,envir=mmTK) #STORE TO PROJECT object
  })
 helptext(tabimportA[1,1],paste0("Choose a frequency file (LRmix/EFM format)."))
 tabimportA[1,2] <- gWidgets2::glabel("none", container=tabimportA) #Create label of whether freq data is imported:
 if(!is.null(get("popFreq",envir=mmTK))) setAsImported( optList$optFreqfile  )
 
 #Set other import buttons
 tabimportA[2,1] = gWidgets2::gbutton(text="Import reference",container=tabimportA,handler=f_importprof,action="ref")
 helptext(tabimportA[2,1],"Imports reference profile(s) from a selected file. \n\nThe column names must contain 'sample..', 'marker', 'allele..'.")

 tabimportA[2,2] = gWidgets2::gbutton(text="Import evidence",container=tabimportA,handler=f_importprof,action="mix")
 helptext(tabimportA[2,2],"Imports evidence  profile(s) from a selected file. \n\nThe column names must contain 'sample..', 'marker', 'allele..', 'height..'.")
 
 tabimportA[3,1] = gWidgets2::gbutton(text="Fill out kits",container=tabimportA,handler=f_fillsettings,action="KITS")
 helptext(tabimportA[3,1],"Filling out kitnames based on which replicates belongs together (id-column).")
 
  
 tabimportA[3,2] = gWidgets2::gbutton(text="Fill out settings",container=tabimportA,handler=f_fillsettings,action="SETTINGS")
 helptext(tabimportA[3,2],"Filling out settings based on which replicates belongs together (id-column).")
 
 tabimportB[1,1] = gWidgets2::gcheckboxgroup(items= getDataNames_type("ref"), container = tabimportB, use.table=TRUE)
 gWidgets2::size(tabimportB[1,1]) <- c(470,150) #c(600,150)

  
 #IMPORT EVIDENCE:
 initKits <- c("NONE",euroformix::getKit()) #get kit names (shortname)
 evids = getDataNames_type("mix")
 if(evids[1]!="") {
   repSettings = get("repSettings",envir=mmTK) #get stored settings about replicates
   evidGrps = repSettings$evidGrp #obtain group info
   selKits = repSettings$kit  #get selected kit (earlier)
  
#   stop()
   nEvids = length(evids) #get number of evidence
   for( e in seq_len(nEvids) ) {
     evidName = evids[e]
     evidGrp <- 1:e #inset stored group ID
     evidSel = evidGrps[e] #set selected
     tabimportC[e,1] = gWidgets2::gcombobox(items=evidGrp, container = tabimportC, editable=TRUE, selected = evidSel)
     gWidgets2::size(tabimportC[e,1]) <- editboxsize[1]
          
     tabimportC[e,2] = gWidgets2::gcheckbox(evidName, container = tabimportC, checked = TRUE)
     gWidgets2::size(tabimportC[e,2]) <- max(nchar(evids)) #editboxsize[2]
     
     evidKit = 0
     if(!is.null(selKits) && evidName%in%names(selKits) )  {
       evidKit = which(initKits==selKits[evidName]) #extract kit
       if(length(evidKit)==0) evidKit = 0 #kit did not exist in dropdown (reseting)
     }
     tabimportC[e,3] <- gWidgets2::gcombobox(items=initKits, width=max(nchar(initKits)), selected = evidKit, editable = FALSE, container = tabimportC)
     tabimportC[e,4] <- gWidgets2::gbutton(text="Settings",container=tabimportC, handler=f_settings, action=evidName)
   }
 }

 #helpfunction used to extract selected importdata-elements to further model-setup
 selectDataToModel <- function(h,....) {
   #All: popFreq must be imported!
   #EVID: must have both mixture and reference profiles
   #DC: Deconvolution requires only mixtures. Reference profiles is optional
   popFreq <- get("popFreq",envir=mmTK)
   mixSel <- refSel <-  numeric()
   opt = get("optSetup",envir=mmTK)  #Get rep settings from envir (includes default values)
   
   if(length(tabimportB[1,1][])>0) refSel <- gWidgets2::svalue(tabimportB[1,1])  #get selected references

   evidD = getData("mix") #get selected mixtures
   evidNames = names(evidD)
   grpID = getRepGrpID( evidNames )
   evidGrp <- grpID$evid
   kitGrp = grpID$kit
   
   evidNames = names(evidD) #obtain evidence names
   repSettings = NULL
   if(length(evidGrp)>0) {
     mixSel = names(evidGrp) #obtain selected evidence 

     #Set repSettings (set values from )
     def = function(x)  setNames(rep(x,length(mixSel)),mixSel)  #set defaul values
     repSettings = list(AT=def(opt$AT[1]),pC=def(opt$pC[1]),lambda=def(opt$lambda[1]),fst=opt$fst[1]) #init if not yet found
     repSettings$kit = def(NA) #default kit defined for each rep
     repSettings$evidGrp=evidGrp #store grouped evidence (last element)
     
     #obtain kit settings from GUI:
     for(evid in mixSel) { #traverse each selected and obtain
    #  evid=mixSel[1]
        grp0 = evidGrp[evid] #obtain group from evidence
        kit0 = na.omit(names(kitGrp)[kitGrp==grp0]) #obtain kits 
        kit0 = setdiff(kit0,"NONE")
        if(length(kit0)>0) repSettings$kit[ evid ] = setNames(kit0[1],evid) #insert selected kit
        
        #Set other settings:
        for(var in names(repSettings)[1:3]) {
          if(!is.na( opt[[var]][evid] )) repSettings[[var]][evid] = opt[[var]][evid] #insert if value was given
        }
      }
    } #end if any evid existing
    assign("repSettings",repSettings,envir=mmTK) #store replicate settings
   
    #indicate which samples that belong together (suggest model settings)
    mod = list()
    parNames = c("mixProp","PHexp","PHvar", "DEG","BWS","FWS") #this are parameter names
    for(par in parNames[1:3]) mod[[par]] = setNames( 1:length(mixSel),mixSel)
    for(par in parNames[4:6]) mod[[par]] = setNames( rep(0,length(mixSel)),mixSel)
    
    #Insert grouped info to default parameters:
    grps = unique(evidGrp)
    nGrps = length(grps) #obtain number of unique grps
    for(g in seq_len(nGrps)) { #for each group
      #g=1
      grp = grps[g] #obtain group
      reps = names(evidGrp)[ evidGrp==grp ] #obtain replicates for specific groups 
      for(par in parNames[1:3]) mod[[par]][reps] = g #insert group info (incremented index)

      kits0 = na.omit(unique(names(kitGrp)[kitGrp==grp])) #obtain kit
      if(length(kits0)==1) mod[[parNames[4]]][reps] = g #insert group info (incremented index)
    }
    assign("modelSettings",mod,envir=mmTK)  #Get rep settings from envir (includes default values)
    
   if(h$action=="SAVE") return() #stop here if only saving
    
   if(is.null(popFreq)) {
     gWidgets2::gmessage("No frequencies was specified!\n Please import table with population frequencies.")
   } else if(length(mixSel)==0) {
     gWidgets2::gmessage("Please import and select evidence-profile!")
   } else if(h$action=="EVID" && length(refSel)==0) {
     gWidgets2::gmessage("Please import and select reference-profiles for weight of evidence!")
   } else {
     refreshTabModel(mixSel,refSel,h$action) #refresh table with selected data
     gWidgets2::svalue(nb) <- 2 #change tab of notebook
   }
 } #end selectDataToModel
 
 #Button-choices further (INTERPRETATION:
 tabimportD[1,1] = gWidgets2::gbutton(text="View data",container=tabimportD,handler=f_viewdata)
 helptext(tabimportD[1,1],"Print all selected data to console. Selected Evidence profiles are shown as EPGs (replicates for those in same group).")

 tabimportD[1,2] = gWidgets2::gbutton(text="Compare data",container=tabimportD,handler=f_compare)
 helptext(tabimportD[1,2],"Obtain number of missmatches between all references and evidences.")
  
 tabimportD[1,3] = gWidgets2::gbutton(text="Weight-of-Evidence",container=tabimportD,handler=selectDataToModel,action="EVID")
 helptext(tabimportD[1,3],"A module for calculating the Likelihood Ratio for selected sample(s) (treated as replicates). \n\nSelected reference(s) can later be conditioned on in the hypotheses.")
 
 tabimportD[1,4] = gWidgets2::gbutton(text="Deconvolution",container=tabimportD,handler=selectDataToModel,action="DC")
 helptext(tabimportD[1,4],"A module for ranking the most likely profiles of the unknown contributors for selected sample(s) (treated as replicates).\n\nThe user will first need to fit a gamma-model based on maximum likelihood estimation.")
 
 tabimportD[1,5] = gWidgets2::gbutton(text="RESTART",container=tabimportD,handler=function(h,...) {
   gWidgets2::dispose(mainwin) #remove window!
   efm2() #and open EFM again
 })
 
 
####################################################################################################################
#######################################Tab 2: Hypotheses:##########################################################
#####################################################################################################################

 
#  mixSel= names(getData("mix"));refSel=names(getData("ref"));type = "EVID"
  refreshTabModel = function(mixSel,refSel,type) { 

    #helpfunction which takes GUI settings and stores them in "set'type'"
    storeSettings = function(type2="EVID") {
      mixData = get("mixData",envir=mmTK) 
      refData = get("refData",envir=mmTK) 
      popFreq = get("popFreq",envir=mmTK)
      repSettings = get("repSettings",envir=mmTK)  #this is settings for each sample (AT,dropin etc)
      optFreq = get("optFreq",envir=mmTK)  #Get default settings from envir
      relSettingsList = get("relSettings",envir=mmTK)  #this is settings for each sample (AT,dropin etc)
      
      #Obtain selected reps:
      mixSel2Ind = numeric() #Keep index
      for(e in 1:nEvids) { #check which reps to continue with
        if(gWidgets2::svalue(tabmodelB[e+2,1])) mixSel2Ind = append(mixSel2Ind, e)
      }
      mixSel2 = mixSel[mixSel2Ind] #get name of selected reps
      
      #Obtain model specification for selected refs:
      parList = list()
      for(i in 1:length(parNames)) parList[[parNames[i]]] = 1:length(mixSel2Ind) #defining parameter index
      
      #going through GUI and re-specify list
      for(c in 1:length(parIndCols)) { #for each param type
        col = parIndCols[c] #obtain column
        for(e in seq_len(length(mixSel2Ind))) { #check which reps to continue with
          row = mixSel2Ind[e] + 2 #obtain row index
          val = as.integer(gWidgets2::svalue(tabmodelB[row,col])) #obtain selected value
          parList[[c]][e] = val
        }
      }
      #Re-structure par-list: Be sure to make proper structure on it
      parList2 = parList #create copy for updated list
      for(i in 1:length(parList)) {
        uniq = setdiff(unique(parList[[i]]),0) #get unique vals (but not zero, those are unaffected)
        for(j in 1:length(uniq)) {
          parList2[[i]][parList[[i]]==uniq[j]] = j #insert increasly ordered index
        }
      }
      
      #Select data to continue with:
      samples=mixData[mixSel2]
      refData=refData[refSel]
      
      #Update par-settings (keep only settings of those selected)
      for(elem in names(repSettings)) {
        tmp = repSettings[[elem]][mixSel2]  #to keep
        if(!any(is.na(tmp))) repSettings[[elem]] = tmp #insert back
      }

      #prepare data for showing output: 
      if(type2%in%c("VIEW","COMPARE")) {
        dat = prepareData2(samples, refData, popFreq, minF=getminFreq(),normalize=as.logical(optFreq$normalize), AT=repSettings$AT )
        if(type2=="VIEW") {
          print(dat) #print data to screen
          print(repSettings)
          #grps = unique(repSettings$evidGrp) #obtain groups
		  
		  #VISUALIZE EPG?
        }
        return(dat) #return with data
      }
      
      #EXTRACT HYPOTHESIS SPECIFICATION (tabmodelA1, tabmodelA2)
      
      #get specified preposition 
      NOC = as.integer(gWidgets2::svalue(tabmodelA1[1,2])) #number of contributors is fixed! 
      checkPosInteger(NOC,"Number of contributors under Hp/Hd")
      nC_hp <- nC_hd <- NOC #number of contributors in model:
      if(type=="DC") nC_hp = NULL #SET TO NULL TO INDICATE THAT Hp is not to be calculated

      condOrder_hp <- condOrder_hd <- rep(0,nRefs)
      knownref_hp <- knownref_hd <- NULL #Typed profiles (known non-contributors)
      for(rsel in refSel) { #traverse each reference
        row = which(rsel==refSel)
        
        #Hd: EVID + DC
        valhd <- as.integer(gWidgets2::svalue(tabmodelA2$hd[row,2]))
        condOrder_hd[row] <- valhd + valhd*max(condOrder_hd)
        
        #Hp: ONLY FOR EVID (NOT DC)
        if(type=="EVID") {
          valhp <- as.integer(gWidgets2::svalue(tabmodelA2$hp[row,2])) 
          condOrder_hp[row] <- valhp +  valhp*max(condOrder_hp)
        }
      }
	  
      if(type=="EVID") { #only for Evidence
        knownref_hp <- which(condOrder_hp==0) #those not conditioned on under Hp
        if(length(knownref_hp)==0) knownref_hp <- NULL
      }
      #THIS IS DONE FOR BOTH DC AND EVID
      knownref_hd <- which(condOrder_hd==0) #those not conditioned on under Hd
      if(length(knownref_hd)==0) knownref_hd <- NULL
      
      #Obtain specified kinship 
      getRelData = function(relSettings) { #small helpfunction to structure data for contLikMLE
        retList = list(knownRel=NULL,ibd=NULL)
        if(!is.null(relSettings)) {
          nrows = nrow(relSettings) #number of definitions
          retList$knownRel = match(relSettings[,2],refSel) #obtain indices
          names(retList$knownRel) = relSettings[,3] #insert type of kinship as header
          ibdList = list()
          for(row in seq_len(nrows)) ibdList[[row]] = as.numeric(relSettings[row,4:6])
          retList$ibd = ibdList
        }
        return(retList)
      }

      #get input to list: note: "fit_hp" and "fit_hd" are list-object from fitted model
      hyp = list(hp=getRelData(relSettingsList$hp),hd=getRelData(relSettingsList$hd))
      if(repSettings$fst>0 && (!is.null(hyp$hp$knownRel) || (!is.null(hyp$hd$knownRel)))) {
        warning("NOT IMPLEMENTED CORRECTLY FOR FST>0")
      }
      hyp$hp = append(hyp$hp ,list(nC=nC_hp,condOrder=condOrder_hp,knownRef=knownref_hp))
      hyp$hd = append(hyp$hd ,list(nC=nC_hd,condOrder=condOrder_hd,knownRef=knownref_hd))
      set <- list(samples=samples,refData=refData,popFreq=popFreq,model=parList2,hyp=hyp,param=repSettings)     

      #get specified preposition 
      nHp =  sum(hyp$hp$condOrder>0) +  length(hyp$hp$knownRel) #number of conditional + related unknowns
      nHd =  sum(hyp$hd$condOrder>0) +  length(hyp$hd$knownRel) #number of conditional + related unknowns
      if( nHp > NOC || nHd > NOC) {
        gWidgets2::gmessage("The number of contributors was speified too low!",title="Wrong setting",icon="error")
        return(0) #return with error
      }
      
      #Store to a listed object
      calcList = get("calcList",envir=mmTK) #obtain already calculated objects
      if(is.null(calcList)) calcList = list()
      calcList[[length(calcList)+1]] = set #insert object
      assign("calcList",calcList,envir=mmTK) #store data to envir
      
      return(1) #success
    } #end store settings from GUI to environment
    
    
    f_compare = function(h,...) { #helpfunction to compare replicates (is comparable or not)
      dat <- storeSettings("COMPARE") #Obtain DATA to compare (from stored settings)

      locs = names(dat) #obtain loci to evaluate
      allSampleNames = unique(unlist(lapply(dat,function(x) names(x$samples)))) #get all samples names 
      sumPHmatrix = NULL
      for(loc in locs) {
        sumPH = sapply( dat[[loc]]$samples, function(x) sum(x)/2)
        sumPHmatrix = rbind(sumPHmatrix, sumPH[allSampleNames])
      }
      colnames(sumPHmatrix) = allSampleNames
      rownames(sumPHmatrix) = locs
      sumPHmatrix[is.na(sumPHmatrix)] <- 0 #insert zero as proxy
      VAR = rowMeans(sumPHmatrix^2)-rowMeans(sumPHmatrix)^2 #obtain variances
      sumPHmatrix = sumPHmatrix[VAR>0,,drop=FALSE] #drop markers with zero variance
      
      repSettings = get("repSettings",envir=mmTK)  #this is settings for each sample (AT,dropin etc)
      kitNames = repSettings$kit[allSampleNames] #obtain selected kit names for each sample

# X=t(sumPHmatrix);grps = kitNames;getClusters=TRUE
      txt = "Similarity of samples (based on summed PHs)"
      clusters = getPCAclusters(X=t(sumPHmatrix),grps = kitNames, txt=txt)
      #print(clusters)      
    }
    
    #common variables
    boxgrp1 = c("Unequal","Equal") #for required params always ON
    boxgrp2 = c("OFF",boxgrp1) #Include also OFF
    boxwidth = max(nchar(boxgrp2))+1 #set fixed width of dropdown
    
    #helpfunction to set box values
    f_setBoxValues = function(h,...) {
      col = h$action #first element is which column
      boxSel = gWidgets2::svalue(tabmodelB[2,col]) #obtain selection typ
      
      idSel = 1:nEvids #get parameter id selections (unequal is default)
      if(boxSel=="OFF")  idSel=0*idSel
      if(boxSel=="Equal")  idSel= 1 + 0*idSel
      for(e in 1:nEvids) { #for each selected replicate (put a dropdown box)
        gWidgets2::svalue(tabmodelB[e+2,col]) = idSel[e]
      }    
    }
    
    #type={EVID",DC"}
    gWidgets2::visible(mainwin) <- FALSE
    # dispose(tabmodel) 
    tabmodeltmp <- gWidgets2::glayout(spacing=spc,container= tabmodel[1,1] <- gWidgets2::ggroup(container=tabmodel) ) 
    tabmodelCC = gWidgets2::glayout(spacing=10,container=(tabmodeltmp[1,1] <-gWidgets2::gframe(spacing=10,container=tabmodeltmp)))  
    tabmodelA = gWidgets2::glayout(spacing=5,container=(tabmodelCC[1,1] <-gWidgets2::gframe("Hypothesis specification",container=tabmodelCC))) 
    tabmodelC = gWidgets2::glayout(spacing=0,container=(tabmodelCC[2,1] <-gWidgets2::gframe("Further:",container=tabmodelCC)))  
    tabmodelB = gWidgets2::glayout(spacing=0,container=(tabmodeltmp[1,2] <-gWidgets2::gframe("Model specification",container=tabmodeltmp))) 
    
    #Obtain evid settings
    repSettings <- get("repSettings",envir=mmTK) #get evidence specific settings from envir
    kits = repSettings$kit #obtain defined kit for eachsample
      
    edwith = 6 #edit width
    nRefs = length(refSel)
    nEvids = length(mixSel)
    NOCsel = nRefs + 1
	
    #Hypothesis selection: subframe of A
    Krange <- 1:maxKsetup #default Contr range
    txt <- "Contributor(s) under H"
    tabmodelA1 = gWidgets2::glayout(spacing=0,container=(tabmodelA[1,1] <-gWidgets2::gframe("Number of contributors",container=tabmodelA))) 
    tabmodelA1[1,1] <- gWidgets2::glabel("NOC:",container=tabmodelA1)
    tabmodelA1[1,2] <- gWidgets2::gcombobox(items=Krange,selected=NOCsel,editable=TRUE,container=tabmodelA1)
    gWidgets2::size(tabmodelA1[1,2]) = 3 #set with
    

#    tabmodelTmp= tabmodelA2$hp
    createHypSetup = function(hypsel,checked=TRUE) { #Helpfunction
      tabmodelTmp = tabmodelA2[[hypsel]]
      for(rsel in refSel) { #indicate all selected references
        rowind = which(rsel==refSel)
        tabmodelTmp[rowind,1]  <- gWidgets2::glabel(paste0("C",rowind),container=tabmodelTmp)
        tabmodelTmp[rowind,2]  <- gWidgets2::gcheckbox(rsel,container=tabmodelTmp,checked=checked)
      }
      
      tabmodelTmp[length(refSel) + 1,2]  <- gWidgets2::gbutton("Add kinship",container=tabmodelTmp,handler=function(h,...){
        NOC = as.integer(gWidgets2::svalue(tabmodelA1[1,2])) #obtain defined NOC
        nCond = 0 #count number of conditionals
        for(rsel in refSel) nCond = nCond + gWidgets2::svalue(tabmodelTmp[which(rsel==refSel),2]) #count number of selected
        nCadd = NOC - nCond #number of contributors extra
        if(nCadd > 0) kinshipSelectorGUI(refSel,Crange=nCond + 1:nCadd ,hypsel,mmTK) #Crange is contribution range of the unknowns
      })
    }
    
    tabmodelA2 = list()
    if(type%in%c("EVID")) {
      tabmodelA2$hp = gWidgets2::glayout(spacing=0,container=(tabmodelA[2,1] <-gWidgets2::gframe(paste0(txt,"p:"),container=tabmodelA))) 
      createHypSetup("hp",TRUE)
    }
    tabmodelA2$hd = gWidgets2::glayout(spacing=0,container=(tabmodelA[3,1] <-gWidgets2::gframe( paste0(txt,"d:"),container=tabmodelA)))
    createHypSetup("hd",TRUE)
    
    #Model specification (replicates): 
    #tabmodelB[1,1] <- gWidgets2::glabel(text="",container=tabmodelB)
    tabmodelB[2,1] <- gWidgets2::glabel(text="Replicate",container=tabmodelB)
    tabmodelB[2,2] <- gWidgets2::glabel(text="Kit",container=tabmodelB)

    #Specify parameter names
    mod = get("modelSettings",envir=mmTK)  #Get rep settings from envir (includes default values)
    parNames = names(mod) #c("mixProp","PHexp","PHvar", "DEG","BWS","FWS") #this are parameter names
    
    parIndCols = 3:8 #column index for parameter option choices
    DEGcol = 6 #column index of degradation
    #provide fast selection assignation (drop-down meny)
    for(col in parIndCols) {
      param = parNames[which(parIndCols==col)] #obtain parameter name
      tabmodelB[1,col] <- gWidgets2::glabel(text=param,container=tabmodelB) #insert param name
      
      boxgrp = boxgrp1
      if(col>=DEGcol) boxgrp = boxgrp2
      tabmodelB[2,col] <- gWidgets2::gcombobox(items=boxgrp,selected=0,editable=FALSE,container=tabmodelB, handler=f_setBoxValues,action=col)
      gWidgets2::size(tabmodelB[2,col]) = boxwidth
    } 
    
    Srange1 = 1:nEvids #must always be selected
    Srange2 = c(0,Srange1) #include also turn off
    for(e in 1:nEvids) { #for each selected replicate (put a dropdown box)
      #tabmodelB[e+2,1] <- gWidgets2::glabel(text=mixSel[e],container=tabmodelB) #insert name
      tabmodelB[e+2,1] <- gWidgets2::gcheckbox(text=mixSel[e],container=tabmodelB,checked=TRUE) #insert name w/checkbox
      for(col in parIndCols) {
        Srange = Srange1
        sel = mod[[ which(col==parIndCols) ]][e] #selecting replicate index from modelSettings object (decided when storingSettings)
        if(col >= DEGcol)  {
          Srange = Srange2
          sel = sel+1 #set as zero if turned off
        }
        tabmodelB[e+2,col] <- gWidgets2::gcombobox(items=Srange,selected=sel,editable=FALSE,container=tabmodelB)
        gWidgets2::size(tabmodelB[e+2,col]) = 4
      }        
      #CHECK AFTER (use defined kit info)
      kit0 = kits[e]
      if(!is.na(kit0)) {
        tabmodelB[e+2,2] <- gWidgets2::glabel(text=kit0,container=tabmodelB) #insert name
      } else { #kit not found
        gWidgets2::svalue( tabmodelB[e+2,DEGcol] ) <- 0 #set selection to zero
        gWidgets2::enabled( tabmodelB[2,DEGcol] ) <- FALSE
        gWidgets2::enabled( tabmodelB[e+2,DEGcol] ) <- FALSE
      }
    } #
    

  #View evaluating evidence/databases 
  tabmodelC[1,1] = gWidgets2::gbutton(text="Print data",container=tabmodelC,handler= function(h,...) { 
    storeSettings("VIEW") #SHOW DATA
   })
  
  tabmodelC[2,1] = gWidgets2::gbutton(text="Compare reps",container=tabmodelC,handler=f_compare)
  helptext(tabmodelC[2,1],"Giving comparison plots of replicates. Which ones that are replicates or not.")
  
 #Calculation button:  
  tabmodelC[3,1] = gWidgets2::gbutton(text="CALCULATE",container=tabmodelC,handler=
  function(h,...) {
    ret = storeSettings(type) #store settings
    if(ret==0) return() #don't continue if value is 0
    refreshTabMLE() #refresh MLE fit tab (i.e. it fits the specified model)
    gWidgets2::svalue(nb) <- 3 #go to mle-fit window (for all cases) when finished
  }) #end cont. calculation button

  
   gWidgets2::visible(mainwin) <- TRUE
   gWidgets2::focus(mainwin) <- TRUE
} #end refresh setup tab-frame


###################################################################################################
############################Tab 3: MLE estimation (RESULTS:########################################
###################################################################################################

  #helpfunction ran when call deconvolution
  doDC <- function(mleobj) {
     dcopt <- get("optDC",envir=mmTK) #options when Deconvolution
     dcobj <- deconvolve2(mlefit=mleobj,alpha=dcopt$alphaprob,maxlist=dcopt$maxlist) 
     DCtable1<-addRownameTable(dcobj$table2)
     colnames(DCtable1)[1] <- "Locus"
     DCtable2<-dcobj$table1
     DCtable3<-dcobj$table3
     DCtable4<-dcobj$table4
     assign("resDC",list(DCtable1,DCtable2,DCtable3,DCtable4),envir=mmTK) #assign deconvolved result to environment
     refreshTabDC() #update table with deconvolved results
     gWidgets2::svalue(nb) <- 4 #go to deconvolution results window (for all cases) when finished     
     gWidgets2::focus(mainwin) <- TRUE
   }

  #helpfunction to translate fitted MLE objects, get and store LR values
  getLRvalues = function(set) { 
    fithp = set$mlefit_hp
    fithd = set$mlefit_hd
    logLRmle <- fithp$fit$loglik - fithd$fit$loglik
    logLikHp = fithp$logLiki #obtain logLik per marker 
    logLikHd = fithd$logLiki #obtain logLik per marker 
    #sum(logLikHp-logLikHd)==logLRmle
    
    LRi <- exp(logLikHp-logLikHd)
    log10LRi <- (logLikHp-logLikHd)/log(10)
    LRmle <- exp(logLRmle)
    log10LRmle <- logLRmle/log(10)

    #Get maximum attainable LR based on random match probability of POI for each markeres (conditional on refs and fst under Hd)
    POIind = fithd$knownRef #obtain known reference under HD (not conditional) 
    if(length(POIind)==1) {
      #hdcond = which(fithd$condOrder>0) #obtain condtional under Hd
      locs = names(set$refData[[POIind]])
      rmp = setNames(rep(NA,length(locs)),locs)
      for(loc in locs) {
        freq = set$popFreq[[loc]]
        if(length(freq)==0) next
        poi = set$refData[[POIind]][[loc]]$adata
        tmp = euroformix::calcGjoint(freq=freq/sum(freq),nU=1) #get genotype prob conditioning on POI and other refs
        ind = which( (tmp$G[,1]==poi[1] & tmp$G[,2]==poi[2]) | (tmp$G[,1]==poi[2] & tmp$G[,2]==poi[1])) #obtain index to use
        if(length(ind)==0) next #skip if none found
        rmp[loc] =  tmp$Gprob[ind][1] #insert geno
      }
      LRupper = 1/prod(rmp, na.rm=TRUE) #prod(scale_fst/rmp) #obtain upper boundary of LR (adjusted with scale)
    } else {
      LRupper = NA
    }
    resEVID <- list(LRmle=LRmle,LRi=LRi,LRupper=LRupper) 
    return(resEVID)
  }

  f_showPerMarkerLR = function(h,...) {
    LRi = h$action
    tab = cbind(LRi,log10(LRi))
    colnames(tab) = c("LR","log10LR")
    showTable(tab,"LR perMarker",2)
  }
  
  #helpfunction to print msg to screen
  doValidMLEModel = function(mlefit,txt="") { #function to get significance level in Validation plot
    alpha <- 0.01 #as.numeric(getValueUser("Set significance level \nin model validation:",0.01))
    #checkPositive(alpha,"The significance level",strict=TRUE)
    valid = validMLEmodel2(mlefit,txt,alpha=alpha,createplot = TRUE, verbose = TRUE) #return table
    return(valid)
  }  
  
#BEGIN MAIN FUNCTION (wrapper:  
  refreshTabMLE = function(resID=NULL) {  #resID is index to show (elemnent in calcList)
    dec <- 4 #number of significant numbers to have in MLE print
    gWidgets2::visible(mainwin) <- FALSE
    
    calcList <- get("calcList",envir=mmTK)  #obtain list with calculations:    
    nCalcs = length(calcList)  #number of calculations done
    if(is.null(resID)) resID = nCalcs #uselast element if not selected
    
    set = calcList[[resID]] #use selected element
    if(is.null(set)) return() #NO RESULTS ARE GIVEN (return...)
   
    f_showParam = function(h,...) {
      if(h$action=="hp") mlefit = set$mlefit_hp
      if(h$action=="hd") mlefit = set$mlefit_hd
      df = mlefit$fit$par2
      showTable(df,paste0("PARAMS (",h$action,")"),3)
    }
    
    f_showPerMarkerLogLik = function(h,...) {
      if(h$action=="hp") mlefit = set$mlefit_hp
      if(h$action=="hd") mlefit = set$mlefit_hd
      df = cbind(logLik=mlefit$logLiki)
      showTable(df,paste0("LogLik per marker (",h$action,")"),2)
    }
    
    #PREPARE CALCULATIONS: DONE IF mlefit_hd is not obtained#
    #take out relevant parameters from stored list
    mod <- set$model #obtain param settings of model
    hyp <- set$hyp      #obtain hypotheses
    par <- set$param #obtain parameters
    calcHp = !is.null( hyp$hp$nC) #whether Hp is considered

    if( is.null(set$mlefit_hd ) ) {
        
       calcMLE = function(hy) { #helpfunction for optimizing likelihood function (hypothesis setup is changed)
         #Obtain optimizing options
         opt <- get("optMLE",envir=mmTK) #options when optimizing (nDone,delta)
         seed0 = opt$seed
         if(seed0==0) seed0=NULL #convert seed to NULL (0 means none)
         minFreq = getminFreq() #obtain minium freq
         print(paste0(opt$nDone," random startpoints with variation ",opt$delta," are applied in the optimizer.")) 
         
         contLikMLE2(hy$nC,set$samples,set$popFreq,set$refData,hy$condOrder, hy$knownRef, 
           kit=par$kit,AT=par$AT,pC=par$pC,lambda=par$lambda,fst=par$fst, 
           mixProp=mod$mixProp,PHexp=mod$PHexp, PHvar=mod$PHvar, DEG=mod$DEG, stuttBW=mod$BWS, stuttFW=mod$FWS,
           minF=minFreq,normalize = as.logical(get("optFreq",envir=mmTK)$normalize),
           difftol=opt$difftol, steptol=opt$steptol, nDone=opt$nDone, maxThreads=opt$maxThreads, delta=opt$delta, seed=seed0,
           knownRel=hy$knownRel, ibd=hy$ibd)
       } 
       
       #fit under hp: (only for evidence)
       mlefit_hp <- NULL #not used otherwise
       if( calcHp  ) { #considering HP
          #nUhp <- mod$nC_hp-sum(mod$condOrder_hp>0) #number of unknowns
          print("Calculating under Hp...")
          time <- system.time({ mlefit_hp <- calcMLE(hyp$hp) })[3]
          print(paste0("Optimizing under Hp took ",format(time,digits=5),"s"))
          if(!is.null(set$mlefit_hp) && set$mlefit_hp$fit$loglik>mlefit_hp$fit$loglik )  mlefit_hp <- set$mlefit_hp #the old model was better
       }
     
       #fit under hd: (does it for all methods)
       print("Calculating under Hd...")
       time <- system.time({    mlefit_hd <- calcMLE(hyp$hd) })[3]
       print(paste0("Optimizing under Hd took ",format(time,digits=5),"s"))
       if(!is.null(set$mlefit_hd) && set$mlefit_hd$fit$loglik>mlefit_hd$fit$loglik )  mlefit_hd <- set$mlefit_hd #the old model was better
  
       fixmsg <- "The specified model could not explain the data.\nPlease re-specify the model."
       if(is.infinite(mlefit_hd$fit$loglik)) {
         gWidgets2::gmessage(fixmsg,title="Wrong model specification (Hd)",icon="error")
       } else if(calcHp && !is.infinite(set$mlefit_hd$fit$loglik) && is.infinite(mlefit_hp$fit$loglik)) {
         gWidgets2::gmessage(fixmsg,title="Wrong model specification (Hp)",icon="error")
       }

       #store MLE result: also store best mle-values once again (possible with re-running optim)
       set$mlefit_hp=mlefit_hp #store fitted mle-fit
       set$mlefit_hd=mlefit_hd #store fitted mle-fit
       
       #Store by-products of calculation:
       if( calcHp  ) set$resLR <- getLRvalues(set) #obtain get LR results directly (but store to object)
       
       #Obtain text to show in GUI
       refNames = unique(names(set$refData)) #extract reference names
       refHp = refNames[ set$hyp$hp$condOrder>0 ]
       refHd = refNames[ set$hyp$hd$condOrder>0 ]
       nUhp = set$hyp$hp$nC - length(refHp)
       nUhd = set$hyp$hd$nC - length(refHd)
       relHp = setNames(refNames[set$hyp$hp$knownRel],names(set$hyp$hp$knownRel))
       relHd = setNames(refNames[set$hyp$hd$knownRel],names(set$hyp$hd$knownRel))
              
       getHypTxt = function(ref,nU,rel) {
         txt=""
         if(length(ref)>0) txt = paste0(txt, paste0(ref,collapse="/"))
         if(nU>0) {
           tmp = paste0(nU," unknown")
           if(nU>1) tmp = paste0(tmp,"s") #plural
           if(length(ref)>0) tmp = paste0(" + ",tmp) #additional to refs
           txt = paste0(txt, tmp)
           
           if(length(rel)>0) {
             reltxt = paste0("C",length(ref) + 1:length(rel),":",names(rel)," to ",rel)
             txt = paste0(txt, "\n(", paste0(reltxt,collapse="/") ,")")
           }
         }
         return(txt)
       }
       set$hypTxt =  c(NA,getHypTxt(refHd,nUhd,relHd)) #create hyp-text
       if( calcHp  ) set$hypTxt[1] = getHypTxt(refHp,nUhp,relHp)
       names(set$hypTxt) = paste0("H",c("p","d"),": ",set$hypTxt) #insert

       #Show model parameter info:
       repNames = names(set$samples) #obtain replicate names
       set$modelTable = matrix( unlist(set$model),nrow = length(repNames),dimnames = list(repNames,names(set$model)) )
    } #ENDING CALCULATIONS
    
    #######################
    #GUI (common under Hd)#
   ###################################
   #helpfunction used to show MLE fit#
   tableMLE <- function(hyp,tabmleX,sig0=2) {
     if(hyp=="hp") mlefit = set$mlefit_hp
     if(hyp=="hd") mlefit = set$mlefit_hd
     
     if(is.null(mlefit)) return()
     tabmleX1 = gWidgets2::glayout(spacing=0,container=(tabmleX[1,1] <-gWidgets2::gframe("Maximum Likelihood value",container=tabmleX))) 
     nparam1 = (mlefit$nC-1)*length(unique(mlefit$mixProp)) + length(unique(mlefit$PHexp)) + length(unique(mlefit$PHvar))
     nparam2 = length(setdiff(unique(mlefit$DEG),0)) + length(setdiff(unique(mlefit$stuttBW),0)) + length(setdiff(unique(mlefit$stuttFW),0))
     nparam = nparam1 + nparam2
     
     tabmleX1[1,1] =  gWidgets2::glabel(text="nContr=",container=tabmleX1)
     tabmleX1[1,2] =  gWidgets2::glabel(text=mlefit$nC,container=tabmleX1)
     
     tabmleX1[2,1] =  gWidgets2::glabel(text="nParam=",container=tabmleX1)
     tabmleX1[2,2] =  gWidgets2::glabel(text=nparam,container=tabmleX1)
     #tabmleX1[3,1] =  gWidgets2::glabel(text="logLik=",container=tabmleX1)
     #tabmleX1[3,2] =  gWidgets2::glabel(text= round(mlefit$fit$loglik,sig0),container=tabmleX1)
     tabmleX1[3,1] =  gWidgets2::glabel(text="adj.loglik=",container=tabmleX1) #show adj.loglik=-AIC/2, where AIC= -2*logLik + 2*nparam -AIC/2
     tabmleX1[3,2] =  gWidgets2::glabel(text=round((mlefit$fit$loglik - nparam),sig0),container=tabmleX1)
     #tabmleX1[5,1] =  gWidgets2::glabel(text="Lik=",container=tabmleX1)
     #tabmleX1[5,2] =  gWidgets2::glabel(text=getSmallNumber(mlefit$fit$loglik,sig0),container=tabmleX1)
     
     tabmleX2 = gWidgets2::glayout(spacing=0,container=(tabmleX[2,1] <-gWidgets2::gframe("Show more",container=tabmleX))) 
     
     tabmleX2[1,1] =  gWidgets2::gbutton(text="Parameters",container=tabmleX2 ,handler=f_showParam,action=hyp)
     tabmleX2[2,1] =  gWidgets2::gbutton(text="logLik per marker",container=tabmleX2,handler=f_showPerMarkerLogLik,action=hyp)
     
     tabmleX3 = gWidgets2::glayout(spacing=0,container=(tabmleX[3,1] <-gWidgets2::gframe("Further Action",container=tabmleX))) 
     tabmleX3[1,1] <- gWidgets2::gbutton(text="Model validation",container=tabmleX3,handler=function(h,...) {
       valid = doValidMLEModel(mlefit,paste0("PP-plot under ",hyp))  #obtain validation table
       if(hyp=="hp") set$mlefit_hp$nFailed <<- sum(valid$Significant)
       if(hyp=="hd") set$mlefit_hd$nFailed <<- sum(valid$Significant)
       calcList[[resID]] <<- set #update set
       assign("calcList",calcList,envir=mmTK) #store validation results to table
      } )
     tabmleX3[2,1] <- gWidgets2::gbutton(text="Deconvolution",container=tabmleX3,handler=function(h,...) { doDC(mlefit) }  )
     #tabmleA3[2,1] <- gWidgets2::gbutton(text="MCMC simulation",container=tabmleA3,handler=function(h,...) { doMCMC(mlefit_hd) } )
   }#end function
     
    #CREATING GUI FOR RESULTS
    
    #This panel is for selecting different calculated results and showing comparison of calculated models
    tabMLEinfo <- gWidgets2::glayout(spacing=spc,container=(tabMLE[1,1] <- gWidgets2::ggroup(container=tabMLE)))
    tabMLEinfo[1,1] <- gWidgets2::glabel("Go to calculation (#id):",container=tabMLEinfo)
    tabMLEinfo[1,2] <- gWidgets2::gcombobox(items=1:nCalcs,container=tabMLEinfo, selected=resID, editable = FALSE, 
      handler = function(h,...) {
        refreshTabMLE( as.integer(gWidgets2::svalue(tabMLEinfo[1,2])) ) #recurse call when selecting her (updating calc)
      })
    gWidgets2::size(tabMLEinfo[1,2]) <- editboxsize[1]
    
    #SHOW MODEL COMPARISON IN RIGHT PANEL:
    tabMLEinfo[1,3] <- gWidgets2::gbutton("Model comparison",container=tabMLEinfo,handler = function(h,...) { 
      #print("COMPARING CALCULATED MODELS:")
      calcLists <- get("calcList",envir=mmTK)  #obtain list with calculations:    
      nCalcs = length(calcLists)  #number of calculations done
      sig0 = 2
      #RUN THROUGH EACH MODEL AND EXTRACT ESSENITAL INFO:
      outtab = NULL
      for(cc in seq_len(nCalcs)) {
        set = calcLists[[cc]] #obtain list
        nparam = length(set$mlefit_hd$fit$phihat) #nparam
        repNames = paste0(rownames(set$modelTable),collapse="/") #get sample name
        modvec = apply(set$modelTable,2,function(x) paste0(x,collapse="/")) #obtain model param overview
        
        nFailedHp <- nFailedHd <- NA #defualt is no calcs
        if(!is.null(set$mlefit_hp$nFailed))  nFailedHp = set$mlefit_hp$nFailed
        if(!is.null(set$mlefit_hd$nFailed))  nFailedHd = set$mlefit_hd$nFailed
        nFailedTxt = paste0(nFailedHp,"/",nFailedHd)
        
        LRval = NA
        if(!is.null(set$resLR$LRmle)) LRval = round(log10(set$resLR$LRmle),sig0)
        newrow = setNames(c(paste0("#",cc),LRval,round(set$mlefit_hd$fit$loglik-nparam,sig0),nparam,modvec,repNames,set$hypTxt,nFailedTxt)
                          ,c("Model","log10LR","adjLogLik","nParam",names(modvec),"Evid","Hp","Hd","nFailed(hp/hd)"))
        outtab = rbind(outtab, newrow)
      }
      ord = order( as.numeric(outtab[,3]),decreasing = TRUE) #order to show models
      compTable = outtab#[ord,,drop=FALSE]
      compTable[is.na(compTable)] = "" #remove NA
      assign("resCompare",compTable,envir=mmTK) #store table
      
      #Show as table in GUI
      width = c(40,70,80,40, rep(70,ncol(compTable)-7), rep(120,3))
      tabSearch_win = gWidgets2::gwindow("Model comparison results",visible = FALSE,width=sum(width),height=min(nrow(compTable)*100,1000))
      tabSearch_GUI = gWidgets2::gtable(compTable ,container=tabSearch_win) #show table to user
      gWidgets2::size(tabSearch_GUI) = list(column.widths= width)
      gWidgets2::svalue(tabSearch_GUI,index=1) = ord[1] #optimInd #highlight best model
      gWidgets2::visible(tabSearch_win) = TRUE
      
    })

    tabMLEinfo[1,4] =  gWidgets2::gbutton(text="Create report",container=tabMLEinfo,handler = function(h,...) { 
      saveTable(getReportText(set, mmTK),"txt")
    })
    
    
    #THis panel will come after showing results for particular calculuation (with buttons etc)
    tabMLEtmp <- gWidgets2::glayout(spacing=spc,container=(tabMLE[2,1] <- gWidgets2::ggroup(container=tabMLE))) 
    tabmleD = gWidgets2::glayout(spacing=0,container=(tabMLEtmp[1,1] <- gWidgets2::gframe("Hypotheses + Model",container=tabMLEtmp,expand=T,fill=T)),expand=T,fill=T) 
    tabmleD[1,1] = gWidgets2::glabel( names(set$hypTxt)[1]  ,container=tabmleD) #Hp
    tabmleD[2,1] = gWidgets2::glabel( names(set$hypTxt)[2]  ,container=tabmleD) #Hd
    
    tabmleD[3,1] = gWidgets2::gbutton("Parameter info",container=tabmleD, handler = function(h,...) {
        showTable(set$modelTable,"Parameter sharing info") })
    helptext(tabmleD[3,1],"Show which replicates sharing parameters")
    
    tabmleA = gWidgets2::glayout(spacing=0,container=(tabMLEtmp[2,1] <- gWidgets2::gframe("Estimates under Hd",container=tabMLEtmp,expand=T,fill=T)),expand=T,fill=T) 
    tableMLE(hyp="hd",tabmleA) #hyp="hd";tabmleX=tabmleA;sig0=2
    

    #INSERT RESULTS BASED ON LR IF EVID:
    if( calcHp ) { #used only for weight-of-evidence
      tabmleB = gWidgets2::glayout(spacing=0,container=(tabMLEtmp[2,2] <-gWidgets2::gframe("Estimates under Hp",container=tabMLEtmp,expand=T,fill=T)),expand=T,fill=T) 
      tableMLE(hyp="hp",tabmleB)
      resLR = set$resLR #obtain possibly already calculated object

      tabmleC1 = gWidgets2::glayout(spacing=0,container=(tabMLEtmp[1,2] <-gWidgets2::gframe("Joint LR",container=tabMLEtmp))) 
      #tabmleC1 = gWidgets2::glayout(spacing=0,container=(tabmleC[1,1] <-gWidgets2::gframe("Joint LR",container=tabmleC))) 
      tabmleC1[1,1] =  gWidgets2::glabel(text="LR=",container=tabmleC1)
      tabmleC1[1,2] =  gWidgets2::glabel(text=format(resLR$LRmle,digits=dec),container=tabmleC1)
      tabmleC1[2,1] =  gWidgets2::glabel(text="log10LR=",container=tabmleC1)
      tabmleC1[2,2] =  gWidgets2::glabel(text=format(log10(resLR$LRmle),digits=dec),container=tabmleC1)
      
      txt = paste0("Upper boundary LR: log10 1/RMP=",format(log10(resLR$LRupper),digits=dec)) #highlight upper boundary
      helptext(tabmleC1[2,2],txt)
      #tabmleC1[3,1] =  gWidgets2::glabel(text="Upper boundary=",container=tabmleC1)
      #tabmleC1[3,2] =  gWidgets2::glabel(text=format(log10(resLR$LRupper),digits=dec),container=tabmleC1)
      tabmleC1[3,1] =  gWidgets2::gbutton(text="LR per-marker",container=tabmleC1,handler = f_showPerMarkerLR,action=resLR$LRi)
    } #end if LR results
    
    gWidgets2::visible(mainwin) <- TRUE
    gWidgets2::focus(mainwin) <- TRUE #focus window after calculations are done
    
    #LAST: STORE CALCULATED RESULTS
    calcList[[resID]] = set #put calculated results back to object
    assign("calcList",calcList,envir=mmTK) #store  
    
  } #end refresh tab-frame of MLE-fit

  refreshTabMLE() #Show already calculted evidence-results when program starts


##############################################################
###############Tab 5: Deconvolution results:##################
##############################################################

 f_savetableDC = function(h,...) {
   if(is.null(DCtable)) {
    gWidgets2::gmessage("There is no deconvolution results available!")
   } else {
    saveTable(DCtable[], "txt") #save deconvolution results
   }
 }
  
  f_saveDCasRef= function(h,...) {
    topRanked <- get("resDC",envir=mmTK)[[1]] #get deconvolved results (top ranked table)
    ncol = ncol(topRanked)
    numRefs = (ncol-1)/3 #number of references to export
    outtab = NULL #get outtabe
    for(r in seq_len(numRefs)) {
      genos = topRanked[,3*(r-1) + 2] #obtain genotypes
      genos =  t(matrix(unlist(strsplit(genos,"/")),nrow=2)) #obtain alleles per genotyes
      new = cbind(paste0("C",r),topRanked[,1],genos)
      outtab = rbind(outtab,new)
    }
    colnames(outtab) = c("SampleName","Marker","Allele 1","Allele 2")
    saveTable(outtab,"csv")
  }
  
 refreshTabDC = function(dctype=1) { #1=table1 (top marginal results),2=table2 (joint results), 3=table3 (all marginal results per genotype), 4=table4 (all marginal results per alleles)
   DCtables <- get("resDC",envir=mmTK) #get deconvolved results
   if(!is.null(DCtables)) {
     DCtable[] = NAtoSign(DCtables[[dctype]])  #update Table
    }
 }

 #CREATE DECONV-GUI 
 tabDCa = gWidgets2::glayout(spacing=1,container=tabDC) #table layout
 tabDCb = gWidgets2::ggroup(spacing=1,container=tabDC,expand=T,fill=T)
 itemvecDC = c("Top Marginal","All Joint","All Marginal (G)","All Marginal (A)")
 tabDCa[1,1] <- gWidgets2::glabel("Select layout:",container=tabDCa)
 tabDCa[1,2] <-  gWidgets2::gradio(items=itemvecDC,selected=1,horizontal=TRUE,container=tabDCa,handler=function(x) {
   refreshTabDC( which(itemvecDC==gWidgets2::svalue(tabDCa[1,2])) )
 })
 tabDCa[2,1] <- gWidgets2::gbutton(text="Save table",container=tabDCa,handler=f_savetableDC)  
 tabDCa[2,2] <- gWidgets2::gbutton(text="Save profiles",container=tabDCa,handler=f_saveDCasRef)  
 
 #ADD DC TABLE
 DCtable = gWidgets2::gtable(items="",multiple = TRUE,container = tabDCb,expand=T,fill=T)
 gWidgets2::add(tabDCb,DCtable,expand=T,fill=T)
 refreshTabDC() #open results when program starts


 for(nbvisit in 4:1) gWidgets2::svalue(nb) <- nbvisit #visit tables 
 getFocus()

} #end funcions
