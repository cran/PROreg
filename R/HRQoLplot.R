HRQoLplot <- function(data,legend=FALSE,title="Short Form-36 Health Survey",dimlabel.cex=NULL,legend.cex=1,linewidth=3,title.cex=1,lty=1){
  
  number <- 0
  number.cex <- rep(1,8)
  
  
maxmin <- data.frame(
  PF=c(20,0),
  RP=c(4,0),
  BP=c(9,0),
  GH=c(20,0),
  VT=c(20,0),
  SF=c(8,0),  
  RE=c(3,0),
  MH=c(13,0))
#colnames(maxmin) <- c("PF \n(20)","RP \n(4)","BP \n(9)","GH \n(20)","VT \n(20)","SF \n(8)","RE \n(3)","MH \n(13)")
colnames(maxmin) <- c("PF (20)","RP (4)","BP (9)","GH (20)","VT (20)","SF (8)","RE (3)","MH (13)")


 
#names(data) <- names(maxmin)
#colnames(data) <- c("PF \n(20)","RP \n(4)","BP \n(9)","GH \n(20)","VT \n(20)","SF \n(8)","RE \n(3)","MH \n(13)")
colnames(data) <- c("PF (20)","RP (4)","BP (9)","GH (20)","VT (20)","SF (8)","RE (3)","MH (13)")


dat <- rbind(maxmin,data)

#radarchart(dat, axistype=1, seg=8, plty=1, vlabels=c("RP", "RE","PF", "VT", "SF","BP","GH","MH"), 
#           title="Short Form-36 Health Survey",pcol=brewer.pal(8,"Set1"),col=brewer.pal(8,"Set1"))

#radarchart(dat, axistype=2, pcol=brewer.pal(8,"Set1"), plty=1, pdensity=30, pfcol=brewer.pal(8,"Set1"),
#           title="Short Form-36 Health Survey")

#radarchart(dat, axistype=3, pty=32, plty=1, axislabcol="grey", na.itp=FALSE,
#           title="Short Form-36 Health Survey")

#radarchart(dat, axistype=1, plwd=1:8, pcol=1, centerzero=TRUE, 
#           seg=4, caxislabels=c("worst", "", "", "", "best"),
#           title="Short Form-36 Health Survey")

radarchart(dat, axistype=number, seg=4, pty=32, plty=lty,plwd=linewidth, na.itp=FALSE,cglcol="black",
           title=title,pcol=brewer.pal(8,"Set1"),vlcex=dimlabel.cex,cex.main=title.cex,axislabcol="black",
           calcex=number.cex)

#Legend
if (legend==TRUE){
  legend("topright", legend=c(rownames(dat[-c(1,2),])),text.col=brewer.pal(8,"Set1"),bty="n",cex=legend.cex,
         lty=lty,lwd=legend.cex,col=brewer.pal(8,"Set1"))
}

}