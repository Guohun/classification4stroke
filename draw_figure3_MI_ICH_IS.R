library(lattice)
new.palette=colorRampPalette(c("black","red","yellow","blue","green","white"),space="rgb")

MyPath="./"
bleedingdata<- read.delim(file=paste(MyPath,"M_p105014_bleeding6.csv",sep=""), header=TRUE, sep=",")
clotdata<- read.delim(file=paste(MyPath,"M_p105014_clot6.csv",sep=""), header=TRUE, sep=",")

B_outPut<-as.matrix(bleedingdata[,1:256])
C_outPut<-as.matrix(clotdata[,1:256])
NewB_outPut<-array(0,c(128,128))
NewC_outPut<-array(0,c(128,128))
for (i in 1:128) {
  {
    for (j in 1:128) {
          {
        NewB_outPut[i,j]=B_outPut[i,j]
        NewC_outPut[i,j]=C_outPut[i,j]
          }
    }    
  }
}

levelplot(outPut[1:49,1:49],col.regions=new.palette(20),main="ICH",xlab="degree sequence: S(1,1)-S(8,16)", ylab="degree sequence: S(1,1)-S(8,16)")
levelplot(NewB_outPut[1:128,1:128],col.regions=new.palette(20),main="GDMI of an ICH (ID:p105014_b6)",xlab="degree sequence: S(1,1)-S(8,16)", ylab="degree sequence: S(1,1)-S(8,16)")
levelplot(NewC_outPut[1:128,1:128],col.regions=new.palette(20),main="GDMI of an IS (ID:p105014_c6)",xlab="degree sequence: S(1,1)-S(8,16)", ylab="degree sequence: S(1,1)-S(8,16)")


B2_group=array(0,c(45,45))
C2_group=array(0,c(45,45))
k=1
l=1
for (i in 1:16){
  for (j in 1:16)  {
    if (i!=j){
      diff=abs(i-j)
      for (s in 1:16){
        for (t in 1:16){
          diff1=abs(s-t)
          diff2=abs(s-t-8)
          if (diff1!=diff&diff1!=diff2&&k<=ncol(B2_group)) {  #&B_outPut[(i-1)*16+j,(s-1)*16+t]<1  
            x<-(i-1)*16+j
            y<-(s-1)*16+t
            if (x!=y&(j!=s-1)&(t!=i-1)){
              B2_group[k,l]<-B_outPut[x,y]
              C2_group[k,l]<-C_outPut[x,y]
              l=l+1
              if (l>ncol(B2_group)){
                k=k+1
                l=1
              }
            }
          } 
        }
      }
    }
  }
}

data_fram <- data.frame(D2=c(B2_group, C2_group),group=c(rep('B',45*45),rep('C',45*45)))
boxplot(D2~group,data=data_fram)
