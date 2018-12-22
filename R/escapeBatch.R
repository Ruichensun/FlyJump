source("escapeAnalysis.R")
source("readAVI8.R")
library(methods)
rdir <- "/oasis/tscc/scratch/tkatsuki/data/Images"
bgstart <- 1
bgend <- 10000
bgskip <- 10
start <- 1
end <- 0
large <- 3000
interval <- 8000
maxdist <- 200
size <- 200
unit <- 0.0875
fps <- 160
maskmovie <- F
speedmovie <- F
objectmovie <- F
moviejp <- T
maskmoviejp <- F
DLO <- F
DLOonly <- F
ram <- 0
gender <- "N"
spthresh <- 3000
thresh <- 50
useres <- F  
options(bitmapType = "cairo")
dir <- commandArgs(TRUE)[2]
filelist <- list.files(dir, pattern=".mp4$")
args <- as.numeric(commandArgs(TRUE)[1])
for(i in args){
  ptm <- proc.time()
  filename <- substr(filelist[i], 1, nchar(filelist[i])-4)
  if(file.exists(paste0(dir, "/", filename, ".avi"))==F){
    com <- paste0("ffmpeg -i ", dir, "/", filelist[i], " -pix_fmt gray -vcodec rawvideo -y ", dir, "/", filename, ".avi")
    system(com)
  }else{
    print("avi file already exists. Skipped conversion.")
  }
  avifile <- paste0(filename, ".avi")
  escapeAnalysis(rdir, dir, avifile, bgstart=bgstart, bgend=bgend, start=start, end=end, interval=interval, maxdist=maxdist, 
                 size=size, unit=unit, fps=fps, large=large, maskmovie=maskmovie, speedmovie=speedmovie, objectmovie=objectmovie,
                 moviejp=moviejp, maskmoviejp=maskmoviejp, DLO=DLO, DLOonly=DLOonly, ram=ram, gender=gender, spthresh=spthresh, thresh=thresh, useres=useres)
  unlink(paste0(dir, "/", avifile))
  print(proc.time() - ptm)
  detach("package:ggplot2", unload=TRUE)
  unloadNamespace("data.table")
  unloadNamespace("reshape2")
  unloadNamespace("scales")
  unloadNamespace("plyr")
}