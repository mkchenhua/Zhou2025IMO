#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
my %opts;
my $VERSION = "v2.1-20140214"; 

GetOptions( \%opts,"i=s","m=s","e=i","o=s","pc=s","g=s","d=s","dim=i","w=f","h=f","lp=s","l=s");
my $usage = <<"USAGE";
       Program : $0   
       Discription:   
       Version : $VERSION
       Contact : guo.yu\@majorbio.com
       Usage :perl $0 [options]		
			-i	 * input distance matrix 
			-o * output dir
			-m	input mapping file if you want set points's color and pch by groups. default:none
			-g	group name in mapping file . default:none
			-e	draw ecclipse number,default :0
			-dim	Number of dimensions of NMDS space,default: 2
                        -w      default:6
                        -h      default:6
			-lp     legend place;default:topright
                        -l      T/F; display samplename or not when -m setted;default:T
	   Example:$0 

USAGE

die $usage if(!($opts{i}&&$opts{o}));
die $usage if($opts{m}&& !$opts{g});
die $usage if(!$opts{m}&& $opts{g});

$opts{e}=defined $opts{e}?$opts{e}:0;
$opts{m}=defined $opts{m}?$opts{m}:"none";
$opts{g}=defined $opts{g}?$opts{g}:"none";
$opts{dim}=defined $opts{dim}?$opts{dim}:2;
$opts{w}=defined $opts{w}?$opts{w}:6;
$opts{h}=defined $opts{h}?$opts{h}:6;
$opts{lp}=defined $opts{lp}?$opts{lp}:"topright";
$opts{l}=defined $opts{l}?$opts{l}:"T";
$opts{i}=~/([^\/]+)$/;
my $bn=$1;

if(! -e $opts{o}){
		`mkdir $opts{o}`;
}
		

open CMD,">cmd.r";
print CMD "

library(vegan)
basename=paste(\"nmds\",\"$bn\",sep=\"_\")
mycol <-c(\"#CD0000\",\"#3A89CC\",\"#769C30\",\"#D99536\",\"#7B0078\",\"#BFBC3B\",\"#6E8B3D\",\"#00688B\",\"#C10077\",\"#CAAA76\",\"#EEEE00\",\"#458B00\",\"#8B4513\",\"#008B8B\",\"#6E8B3D\",\"#8B7D6B\",\"#7FFF00\",\"#CDBA96\",\"#ADFF2F\")
#mypch <-c(21,22,24,23,25,11,13,8)
mypch <-c(21:25,3,4,7,9,8,10,15:18,0:14)
pch=21
col=\"#1E90FF\"

# read otubased metics 
da <-read.table(\"$opts{i}\",sep=\"\\t\",head=T,check.names=F)
rownames(da) <-da[,1]
da <-da[,-1]

# read sample design file
map=\"$opts{m}\"
if(map !=\"none\"){
		sd <-read.table(\"$opts{m}\",head=T,sep=\"\\t\",comment.char = \"\",check.names = FALSE)
		rownames(sd) <- as.character(sd[,1])
		sd[,1] <-as.character(sd[,1])
		sd\$$opts{g} <-as.character(sd\$$opts{g} )
		legend <- as.matrix(unique(sd\$$opts{g})) 

		da <-da[as.character(sd[,1]),]
		da <-da[,as.character(sd[,1])]	
		basename =paste(basename,\"$opts{g}\",sep=\"_\")
}		

# if nmds analysis

pc_x =1
pc_y =2

da <-as.dist(da)
nmds <-metaMDS(da,k=$opts{dim})
pc12 <- nmds\$points
sites=paste(\"$opts{o}/\",basename,\"_sites.xls\",sep=\"\")
#rotat=paste(\"$opts{o}/\",basename,\"_rotation.xls\",sep=\"\")

write.table(nmds\$points,sites,sep=\"\\t\")

if(map !=\"none\"){
		#   set class color and pch by default mycol & mypch
		class_count <-as.matrix(table(sd\$$opts{g}))
		class_color  <-mycol[1:(length(class_count))]
		class_pch <- mypch[1:(length(class_count))]
		class <-data.frame(count=class_count,color=as.character(class_color),pch=class_pch)
		col=as.character(class[sd[rownames(pc12),]\$$opts{g},]\$color)
		pch=class[sd[rownames(pc12),]\$$opts{g},]\$pch  
		lcol <-as.vector(class[legend,]\$color)
		lpch <-as.vector(class[legend,]\$pch)		
}		


# plot nmds
pdf=paste(\"$opts{o}/\",basename,\".pdf\",sep=\"\")
pdf(pdf,width=$opts{w},height=$opts{h})   # save graph to pdf
par(mar=c(6,6,2,2))
mex<-0.2*abs(max(pc12[,1])-min(pc12[,1])) 
mey<-0.2*abs(max(pc12[,2])-min(pc12[,2]))
tmpp <- nmds\$stress
#head1 <- paste(\"NMDS\",\"\(\",\"stress:\",tmpp,\"\)\",sep=\"\")
plot(pc12,xlim=c(min(pc12[,1])-mex,max(pc12[,1])+mex),ylim=c(min(pc12[,2])-mey,max(pc12[,2])+mey),xlab=\"nmds1\",ylab=\"nmds2\",main=paste(\"NMDS\",\"\(\",\"stress:\",signif(tmpp,3),\"\)\",sep=\"\"),cex=1,las=1,pch=pch,col=col,bg=paste(col,\"FF\",sep=\"\"))

# point label
library(\"maptools\")
label=\"$opts{l}\"
  if(label==\"T\"){
          pointLabel(x=pc12[,1],y=pc12[,2],labels=paste(\"\\n    \",rownames(pc12),\"    \\n\",sep=\"\"),cex=0.7,col=col)
                  }else{}       

if(map !=\"none\"){
	if(length(legend)>1){
		legend(\"$opts{lp}\",legend=legend,col=lcol,pch=lpch,pt.bg=paste(lcol,\"FF\",sep=\"\"))
	}
}

ell =$opts{e}
if(ell!=0){
		# draw ellipse by cluster analysis 
		library(\"vegan\")

		da.dist <-vegdist(da,method=\"bray\")

		hc <-hclust(da.dist,method=\"complete\")
		library(\"cluster\")
		n=$opts{e}  # set cluster number 
		cn <-cutree(hc,k=n)
		cu <- sapply(1:n,function(x) names(cn[cn==x]))

		# plot ellipses
		drawelli <-function(xy){
			if(length(xy)>4){
				eh <-ellipsoidhull(xy)
				lines(ellipsoidPoints(eh\$cov,eh\$d2,eh\$loc),col=\"#BEBEBE\")
			}else{ print(\"less than  two point,ignore.\")}
		}
		sapply(1:n,function(i) drawelli(pc12[unlist(cu[i]),]))
}

dev.off()   # plot over

";

`R --restore --no-save < cmd.r`;
