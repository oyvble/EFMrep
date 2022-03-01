#rm(list=ls());library(EFMrep);library(euroformix);library(testthat)
library(testthat)
library(euroformix)

#Path must set manually if using test_check direclty:
#setwd("C:\\Users\\oyvbl\\Dropbox\\Forensic\\MixtureProj\\myDev\\EFMrep\\tests") 

source("helpfunctions\\getLogLiki.R")
source("helpfunctions\\calcGjoint.R")
source("helpfunctions\\testFuns.R")

test_check("EFMrep")
