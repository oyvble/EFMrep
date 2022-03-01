#' @title kinshipSelectorGUI
#' @author Oyvind Bleka
#' @description GUI window to specify kinship 
#' @param refs Name of references
#' @param Crange range of contributors to allow
#' @param hypsel name of hypothesis to consider (obtain object from env)
#' @param env an environment object from efm2
#' @export 

#Helpfunction to (popup window)
kinshipSelectorGUI = function(refs,Crange,hypsel="hd",env=NULL) {

  relSettingsList = get("relSettings",envir=env) #obtain specified kinship object 
  relSettings = relSettingsList[[hypsel]] #obtain object of specific hypothesis
  
  #CREATE TABLE WITH RELATEDNESS (versus IBD coefficients)
  kinshipIBD = rbind( c(1,0,0), t(replicate(2,c(0,1,0))) , c(1/4,1/2,1/4),  t(replicate(5,c(1/2,1/2,0))) ,c(3/4,1/4,0), c(0,0,1) )  #Defined at https://strbase.nist.gov/pub_pres/Lewis-Towson-kinship-Apr2010.pdf
  rownames(kinshipIBD) <- c("Unrelated" , "Parent","Child", "Sibling" , "Uncle","Nephew","Grandparent","Grandchild","Half-sibling" , "Cousin" , "Twin (ident.)")
  refs2 <- c("",refs) #selection list of references
  
  inwin = gWidgets2::gwindow("Specify kinship for unknown contributors",visible=FALSE,expand=TRUE,fill=TRUE)
  gg = gWidgets2::ggroup(container=inwin,horizontal = FALSE) #set up window layout
  butlay <- gWidgets2::glayout(container=gWidgets2::gframe( paste0("Save selection") ,container=gg)) 

  grid <- gWidgets2::glayout(container=gWidgets2::gframe( paste0("Specify kinship") ,container=gg,expand=TRUE,fill=TRUE),expand=TRUE,fill=TRUE,spacing=10) #evidence,ref dataframe
  for(row in seq_len(length(Crange))) {
    grid[row,1] = paste0("C",Crange[row],":")
    grid[row,2] = gWidgets2::gcombobox(items=rownames(kinshipIBD),container=grid,editable=FALSE) 
    grid[row,3] = gWidgets2::glabel(text="to",container=grid)
    grid[row,4] = gWidgets2::gcombobox(items=refs2,container=grid,editable=FALSE)
    
    if(!is.null(relSettings)) { #if previously specifed
      indrow = which(as.integer(relSettings[,1])==Crange[row]) #check if already set
      if( length(indrow)>0) {
        gWidgets2::svalue(grid[row,2]) = relSettings[indrow,3] #update Ref
        gWidgets2::svalue(grid[row,4]) = relSettings[indrow,2] #update type of rel
      }
    }
  }
  
  f_quit = function(h,...) {
    if(h$action=="save")  {
      relSettings = NULL #store settings from GUI into list (refreshed)
      for(row in seq_len(length(Crange))) { #traverse rows and obtain info
        relType = gWidgets2::svalue(grid[row,2])
        relRef = gWidgets2::svalue(grid[row,4])
        relInd = which(rownames(kinshipIBD)==relType)
        if(length(relRef)==0) relInd = 1 #no reference selected
        
        if(relInd==1) next # skip to next if none selected
        relSettings = rbind(relSettings, c(Crange[row],relRef,relType, kinshipIBD[relInd,]) )
      }
      if(length(relSettings)>0)  colnames(relSettings) = c("C","Ref","Rel","K0","K1","K2") #insert column name if a tale
      #STORE SPECIFICATIONS FROM TABLE          
      relSettingsList[[hypsel]] = relSettings #update list
      assign("relSettings",relSettingsList,envir=env)  #this is settings for kinship (list of hp and hd separately). Object is modfied in function "kinshipSelectorGUI"
    }
    
    bool = gWidgets2::gconfirm("Are you sure you want to exit?")
    if(!bool) return() #
    gWidgets2::dispose(inwin) #close window
  }
  butlay[1,1] = gWidgets2::gbutton("Save and exit",container=butlay,handler=f_quit,action="save")
  butlay[1,2] = gWidgets2::gbutton("Exit without saving",container=butlay,handler=f_quit,action="quit")

  gWidgets2::visible(inwin) = TRUE
}