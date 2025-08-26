#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
my %opts;
use FindBin qw($Bin);

my $VERSION = "2024-03-19";
GetOptions (\%opts,"d=s","type=s","i=s","gr=s","t=s","filter=s","top=i","tv=f","color=s","n=i","so=s","m=i","g=s","b_ncol=i","b_ncol_w=i","b_blas=i","b_w=f","b_h=f","p_w=i","b_scex=f","p_h=f","b_hi=s","b_wi=s","b_lcex=f","b_l1=f","b_lx=s","p=f","s=i","p_ltp=s","od=s","bar=s","pie=s","p_cex=f","p_xlim=f","p_ylim=f","clockwise=s","angle=i","horiz=s","oth=s","space=f","width=f","b_w_w=f","b_h_w=f","b_bo=f");


my $usage = <<"USAGE";
        Program : $0
        Version : $VERSION
        Contact : huifang.an\@majorbio.com
        Discription: plot bar and pie chart
        Usage:perl $0 [options] 
             base options:
                -d *	work directory
                -i *	inputfile,a frame.
		-filter filter by p/top,p:num/all<=p push to Low_abundance,top:top abundance ,default:p
		-top     must > 0 integer ,invalid when -filter p default:20
                -color  color list file
                -gr*       map file ;
                -type     group/mean,default group
                -n	collect the max n in each sample
                -m	max m in all samples' sum
                -p	if num/all<=p push to   Low_abundance invalid when -filter top,  defalt :0.01
                -s	if num<=s  push to Low_abundance     defalt :0  
                -od	sort table with samples sum, [d/i/n],d:decreasing; i:increasing; n:not do anything; default:d
                -so	sort table by sample names. default:T
                -g	calculte sum of groups for samples eg. [   a1 A
                                              				a2 A
                                              				a3 A
                                              				b1 B
                                              				b2 B ...
                                           			    ]
                -bar  [T/F]plot bar or not . defult :T  
                -pie [T/F]plot pie for each sample or not . defult :T   
                
             bar options :                        			       
                -b_ncol  colums of legend defalt=3
		-b_ncol_w	colums of legend2 defalt=1
                -b_lcex  legend cex defalt:0.8
                -b_scex  sample name cex defalt:1
                -b_blas  sample name direction.[1,2,3,4]
                -b_l1    make the width of legend colume1 large defalt:1.3
                -b_w  width  defalt=9
                -b_h  height    defalt=8
                -b_hi split of height 2:1
		-b_wi split of widths 1:1
		-b_w_w	width  defalt=8
		-b_h_w	height    defalt=5
                -b_lx [T/F] make legend looks good  defalt :F 
		-space space between bars,default:1.2
		-width	bar width:default:1
		-b_bo   bar bottom space;default:3
             pie options :
                -p_w  	width      default=7
                -p_h        height    default=5
                -p_cex     label cex .defalt:0.8
                -p_ltp          label type [t/l],t:text around pie ;l :legend on the pie right.defalt :t 
                -p_xlim		xlim of pie;default:0.7
                -p_ylim	    ylim of pie:defualt:0.7
                -clockwise	[T/F] default: F
                -angle		[0-360] default :0
		-horiz          [T/F]  defult :F 
		-oth            [T/F]  plot other or not.       defult :T
USAGE
die $usage if ( !(defined $opts{d}&&$opts{i}&&$opts{gr}));

#define defalts
#$opts{gr}=defined $opts{gr}?$opts{gr}:"none";###map
$opts{t}=defined $opts{t}?$opts{t}:"F";###map
$opts{tv}=defined $opts{tv}?$opts{tv}:"0.5";###map
$opts{filter}=defined $opts{filter}?$opts{filter}:"p";###map
$opts{top}=defined $opts{top}?$opts{top}:20;###map



$opts{n}=defined $opts{n}?$opts{n}:-1;
$opts{m}=defined $opts{m}?$opts{m}:-1;
$opts{p}=defined $opts{p}?$opts{p}:0.01;
$opts{s}=defined $opts{s}?$opts{s}:-1;
$opts{b_ncol}=defined $opts{b_ncol}?$opts{b_ncol}:3;
$opts{b_ncol_w}=defined $opts{b_ncol_w}?$opts{b_ncol_w}:1;
$opts{b_lcex}=defined $opts{b_lcex}?$opts{b_lcex}:0.8;
$opts{b_scex}=defined $opts{b_scex}?$opts{b_scex}:1;
$opts{p_cex}=defined $opts{p_cex}?$opts{p_cex}:0.8;
$opts{b_blas}=defined $opts{b_blas}?$opts{b_blas}:1;
$opts{b_l1}=defined $opts{b_l1}?$opts{b_l1}:1.3;
$opts{b_w}=defined $opts{b_w}?$opts{b_w}:9;
$opts{b_h}=defined $opts{b_h}?$opts{b_h}:8;
$opts{b_w_w}=defined $opts{b_w_w}?$opts{b_w_w}:9;
$opts{b_h_w}=defined $opts{b_h_w}?$opts{b_h_w}:7;
$opts{p_w}=defined $opts{p_w}?$opts{p_w}:8;
$opts{p_h}=defined $opts{p_h}?$opts{p_h}:9;
$opts{b_hi}=defined $opts{b_hi}?$opts{b_hi}:"2:1";
$opts{b_wi}=defined $opts{b_wi}?$opts{b_wi}:"1:1";
$opts{space}=defined $opts{space}?$opts{space}:1.2;
$opts{width}=defined $opts{width}?$opts{width}:1;
$opts{b_bo}=defined $opts{b_bo}?$opts{b_bo}:3;
$opts{g}=defined $opts{g}?$opts{g}:"ALL";
$opts{b_lx}=defined $opts{b_lx}?$opts{b_lx}:"F";
$opts{od}=defined $opts{od}?$opts{od}:"d";
$opts{so}=defined $opts{so}?$opts{so}:"F";
$opts{p_ltp}=$opts{p_ltp}?$opts{p_ltp}:"t";
$opts{bar}=$opts{bar}?$opts{bar}:"T";
$opts{pie}=$opts{pie}?$opts{pie}:"T";
$opts{p_xlim}=defined $opts{p_xlim}?$opts{p_xlim}:0.7;
$opts{p_ylim}=defined $opts{p_ylim}?$opts{p_ylim}:0.7;

$opts{clockwise}=defined $opts{clockwise}?$opts{clockwise}:"F";
$opts{angle}=defined $opts{angle}?$opts{angle}:0;
$opts{horiz}=$opts{horiz}?$opts{horiz}:"F";
$opts{oth}=$opts{oth}?$opts{oth}:"T";
$opts{color}=$opts{color}?$opts{color}:"FALSE";
$opts{type}=$opts{type}?$opts{type}:"group";

my $work_dir = $opts{d};
chdir($work_dir);
print "work_dir=$work_dir\n";

if($opts{type} eq "group") {
my $labtext;
my $legend;
if($opts{p_ltp}=~/^t$/){
          $labtext="TRUE";$legend="FALSE";
}elsif($opts{p_ltp}=~/^l$/){
          $labtext="FALSE";$legend="TRUE";
}else{
          print "opts -ltp must be 't' or 'l'!\n";exit;
}

$opts{g}=~/\/*([^\/]+)$/;
my $gs=$1;
my $del=0;
open RCMD, ">cmd.r";

print RCMD "
#mycol <-c(119,132,147,454,89,404,123,529,463,104,552,28,54,84,256,100,558,43,652,31,610,477,588,99,81,503,562,76,96,495)
#mycol <- c(34, 51, 142, 23, 50, 27, 31, 75, 525, 62, 119, 46, 475, 554, 622, 483, 657, 545, 402, 477, 503, 40, 115, 5, 376,473,546,482)
mycol <- c(34, 51, 142, 26, 31, 371, 36, 7, 12, 30, 84, 88, 116, 121, 77, 56, 386, 373, 423, 435, 438, 471, 512, 130, 52, 47, 6, 11, 43, 54, 367, 382, 422, 4, 8, 375, 124, 448, 419, 614, 401, 403, 613, 583, 652, 628, 633, 496, 638, 655, 132, 503, 24)
mycol <-colors()[rep(mycol,20)]

otu <-read.table(file=\"$opts{i}\",header=T,quote=\"\",comment.char=\"\",check.names=FALSE,sep=\"\\t\")
rownames(otu) <- otu[,1]
rownames(otu) <-sapply(rownames(otu),function(x) gsub(\"_*{\.+}\",\" \",x,perl = TRUE)) 

otu <-otu[,-1]

#al <- which(rownames(otu) \%in% c(\"All\",\"No_Rank\",\"Trimed\"))
al <- which(rownames(otu) \%in% c(\"All\"))
if(length(al)) otu <-otu[-al,]

#gs <-\"$opts{gr}\"
#if(gs!=\"ALL\"){
#      group <- read.table(\"$opts{g}\",skip=1)
#      glst <- lapply(1:length(unique(group[,2])),function(x)group[which(group[,2] \%in\% unique(group[,2])[x]),1])
#      names(glst) <-as.character(unique(group[,2]))
#      tab <-sapply(1:length(glst),function(x) apply(otu[as.character(as.vector(glst[[x]]))],1,sum))
#      otu <-tab[apply(tab,1,function(x)any(x>0)),]      
#      colnames(otu) <-unique(group[,2])
#	  rm(glst)
#	  rm(tab)
#}

n <-$opts{n}
m <-$opts{m}
p <-$opts{p}
s <-$opts{s}
od <-\"$opts{od}\"
oth <-\"$opts{oth}\"
filter <- \"$opts{filter}\"
top <- $opts{top}
rowsum <-sapply(1:nrow(otu),function(x) sum(otu[x,]))

if(od ==\"d\"){
	otu<-otu[order(rowsum,decreasing=TRUE),]
}else if(od ==\"i\"){
	otu<-otu[order(rowsum,decreasing=FALSE),]
}

so <-\"$opts{so}\"
if(so==\"T\"){
	otu<-otu[,order(colnames(otu),decreasing=FALSE)]
}


if(n>0){
     #otu[order(otu[,1],decreasing=T)[1:10],1]
	
     max_names <- sapply(1:ncol(otu),function(x) rownames(otu)[order(otu[,x],decreasing=T)[1:n]])
     otu <-otu[which(rownames(otu) \%in% unique(as.vector(max_names))),]
     
     #unique(as.vector(max_names))
     #use_list <-which(rownames(otu) \%in% unique(as.vector(max_names)))
     #rownames(otu)[which(rownames(otu) \%in% unique(as.vector(max_names)))]
    
}else if(m>0){
        otu<-otu[order(rowsum,decreasing=TRUE),]
		otu<-otu[1:m,] 
}else if((p >0) & (filter == \"p\")){
	if(oth==\"T\"){
	 otu_pec <- otu
	 otu_pec <- sapply(1:ncol(otu),function(x) otu_pec[,x] <-otu[,x]/sum(otu[,x]))
	 minp <-sapply(1:nrow(otu_pec),function(y) all(otu_pec[y,]<=p)) 	 
	 otu_xp <-otu[minp,]
	 other <- sapply(1:ncol(otu_xp),function(y) sum(otu_xp[,y]))
	 otu <-rbind(otu[!minp,],other)
	 rownames(otu)[nrow(otu)] <-\"Others\"
	 rm(otu_pec)
	}else{
	 otu_pec <- otu
         otu_pec <- sapply(1:ncol(otu),function(x) otu_pec[,x] <-otu[,x]/sum(otu[,x]))
         minp <-sapply(1:nrow(otu_pec),function(y) all(otu_pec[y,]<=p))          
         otu_xp <-otu[minp,]
         other <- sapply(1:ncol(otu_xp),function(y) sum(otu_xp[,y]))
         otu <-rbind(otu[!minp,])
         rm(otu_pec)
	}
}else if(s >0){
	 minp <-sapply(1:nrow(otu),function(y) all(otu[y,]<=s)) 	 
	 otu_xp <-otu[minp,]
	 other <- sapply(1:ncol(otu_xp),function(y) sum(otu_xp[,y]))
	 otu <-rbind(otu[!minp,],other)
	 rownames(otu)[nrow(otu)] <-\"Others\"	 
	 rm(otu_xp)
}else if ((top >0) & (filter == \"top\")){
if ($opts{top} >=  dim(otu)[1]){
top <- dim(otu)[1]}else {
top <- top } 
otu <- otu[names(sort(apply(otu,1,sum),decreasing=T)),]
#otu<-otu[sort(rowsum,decreasing=TRUE),]
otu1<-otu[1:top,]
if (top < dim(otu)[1]) {
otu_xp<-otu[(top+1):nrow(otu),]
other <- sapply(1:ncol(otu_xp),function(y) sum(otu_xp[,y])) 
otu <-rbind(otu[1:top,],other)
rownames(otu)[nrow(otu)] <-\"Others\"}
##rm(otu_pec)
}



del <-unlist(sapply(1:ncol(otu),function(x) if(sum(otu[,x])==0) x))

if(!is.null(del)) {
  mydel <-colnames(otu)[c(del)]
  otu <-otu[,-c(del)]
}

bar <-\"$opts{bar}\"
pie <-\"$opts{pie}\"
horiz <-\"$opts{horiz}\"
################# plot  bar     ########################################
if(bar ==\"T\"){
	dat <-sapply(1:ncol(otu),function(x) otu[,x]/sum(otu[,x]))
	colnames(dat) <-colnames(otu)
	rownames(dat) <-rownames(otu)
	dat1 <-sapply(1:ncol(otu),function(x) otu[,x])
	colnames(dat1) <-colnames(otu)
	rownames(dat1) <-rownames(otu)
        dat1 <- as.data.frame(dat1)
        dat1\$ID <- rownames(dat1)
        dat1 <- dat1[,c(\"ID\",colnames(dat1)[1:(ncol(dat1)-1)])]
        #write.table(dat1,\"$opts{g}.new.$opts{i}\",sep=\"\\t\",eol=\"\\n\",quote=FALSE)
        write.table(dat1,\"all.bar\",sep=\"\\t\",eol=\"\\n\",quote=FALSE,row.names=F)
	}
";


`R --restore --no-save < cmd.r`;
#system ('rm cmd.r');
system ('rm Rplots.pdf') if -e "Rplots.pdf";


#my $temp = "$opts{g}.new.$opts{i}";
my $xlables=0;
if ($opts{b_blas} == 1) {$xlables=0;}
if ($opts{b_blas} == 2) {$xlables=90;}
my $leg = $opts{b_lcex} * 10;

my $yabb = $opts{b_scex} * 10 ; 
#if ($opts{color} ne "FLASE") {
#print "python $Bin/barplot_for_group.py  -i all.bar -g $opts{gr}  -hi $opts{b_h} -wd $opts{b_w} -ang  $xlables -leg 9  -color $opts{color}\n";

{`python $Bin/barplot_for_group.py  -i all.bar -g $opts{gr}  -hi $opts{b_h} -wd $opts{b_w} -ang  $xlables -leg $leg  -labCex $yabb   -color $opts{color} -b_ncol $opts{b_ncol} `; }
#if ($opts{t} eq "T"){`python $Bin/chongji-lianxian_barplot.py -i all.bar -ratio $opts{tv}  -hi $opts{b_h} -wd $opts{b_w} -ang  $xlables -leg 9 -ylabCex 9 -color $opts{color}`;}
}


else {

my $labtext;
my $legend;
if($opts{p_ltp}=~/^t$/){
          $labtext="TRUE";$legend="FALSE";
}elsif($opts{p_ltp}=~/^l$/){
          $labtext="FALSE";$legend="TRUE";
}else{
          print "opts -ltp must be 't' or 'l'!\n";exit;
}

$opts{g}=~/\/*([^\/]+)$/;
my $gs=$1;
my $del=0;
open RCMD, ">cmd.r";

print RCMD "
#mycol <-c(119,132,147,454,89,404,123,529,463,104,552,28,54,84,256,100,558,43,652,31,610,477,588,99,81,503,562,76,96,495)
#mycol <- c(34, 51, 142, 23, 50, 27, 31, 75, 525, 62, 119, 46, 475, 554, 622, 483, 657, 545, 402, 477, 503, 40, 115, 5, 376,473,546,482)
mycol <- c(34, 51, 142, 26, 31, 371, 36, 7, 12, 30, 84, 88, 116, 121, 77, 56, 386, 373, 423, 435, 438, 471, 512, 130, 52, 47, 6, 11, 43, 54, 367, 382, 422, 4, 8, 375, 124, 448, 419, 614, 401, 403, 613, 583, 652, 628, 633, 496, 638, 655, 132, 503, 24)
mycol <-colors()[rep(mycol,20)]

otu <-read.table(file=\"$opts{i}\",header=T,quote=\"\",comment.char=\"\",check.names=FALSE,sep=\"\\t\")
rownames(otu) <- otu[,1]
rownames(otu) <-sapply(rownames(otu),function(x) gsub(\"_*{\.+}\",\" \",x,perl = TRUE)) 

otu <-otu[,-1]

#al <- which(rownames(otu) \%in% c(\"All\",\"No_Rank\",\"Trimed\"))
al <- which(rownames(otu) \%in% c(\"All\"))
if(length(al)) otu <-otu[-al,]




gs <-\"$opts{gr}\"
if(gs!=\"ALL\"){
      group <- read.table(\"$opts{gr}\",skip=1)
      glst <- lapply(1:length(unique(group[,2])),function(x)group[which(group[,2] \%in\% unique(group[,2])[x]),1])
      names(glst) <-as.character(unique(group[,2]))
      tab <-sapply(1:length(glst),function(x) apply(otu[as.character(as.vector(glst[[x]]))],1,mean))
      otu <-tab[apply(tab,1,function(x)any(x>0)),]      
      colnames(otu) <-unique(group[,2])
	  rm(glst)
	  rm(tab)
}

n <-$opts{n}
m <-$opts{m}
p <-$opts{p}
s <-$opts{s}
od <-\"$opts{od}\"
oth <-\"$opts{oth}\"
rowsum <-sapply(1:nrow(otu),function(x) sum(otu[x,]))

if(od ==\"d\"){
	otu<-otu[order(rowsum,decreasing=TRUE),]
}else if(od ==\"i\"){
	otu<-otu[order(rowsum,decreasing=FALSE),]
}

so <-\"$opts{so}\"
if(so==\"T\"){
	otu<-otu[,order(colnames(otu),decreasing=FALSE)]
}


if(n>0){
     #otu[order(otu[,1],decreasing=T)[1:10],1]
	
     max_names <- sapply(1:ncol(otu),function(x) rownames(otu)[order(otu[,x],decreasing=T)[1:n]])
     otu <-otu[which(rownames(otu) \%in% unique(as.vector(max_names))),]
     
     #unique(as.vector(max_names))
     #use_list <-which(rownames(otu) \%in% unique(as.vector(max_names)))
     #rownames(otu)[which(rownames(otu) \%in% unique(as.vector(max_names)))]
    
}else if(m>0){
        otu<-otu[order(rowsum,decreasing=TRUE),]
		otu<-otu[1:m,] 
}else if(p >0){
	if(oth==\"T\"){
	 otu_pec <- otu
	 otu_pec <- sapply(1:ncol(otu),function(x) otu_pec[,x] <-otu[,x]/sum(otu[,x]))
	 minp <-sapply(1:nrow(otu_pec),function(y) all(otu_pec[y,]<=p)) 	 
	 otu_xp <-otu[minp,]
	 other <- sapply(1:ncol(otu_xp),function(y) sum(otu_xp[,y]))
	 otu <-rbind(otu[!minp,],other)
	 rownames(otu)[nrow(otu)] <-\"Others\"
	 rm(otu_pec)
	}else{
	 otu_pec <- otu
         otu_pec <- sapply(1:ncol(otu),function(x) otu_pec[,x] <-otu[,x]/sum(otu[,x]))
         minp <-sapply(1:nrow(otu_pec),function(y) all(otu_pec[y,]<=p))          
         otu_xp <-otu[minp,]
         other <- sapply(1:ncol(otu_xp),function(y) sum(otu_xp[,y]))
         otu <-rbind(otu[!minp,])
         rm(otu_pec)
	}
}else if(s >0){
	 minp <-sapply(1:nrow(otu),function(y) all(otu[y,]<=s)) 	 
	 otu_xp <-otu[minp,]
	 other <- sapply(1:ncol(otu_xp),function(y) sum(otu_xp[,y]))
	 otu <-rbind(otu[!minp,],other)
	 rownames(otu)[nrow(otu)] <-\"Others\"	 
	 rm(otu_xp)
}

del <-unlist(sapply(1:ncol(otu),function(x) if(sum(otu[,x])==0) x))

if(!is.null(del)) {
  mydel <-colnames(otu)[c(del)]
  otu <-otu[,-c(del)]
}

bar <-\"$opts{bar}\"
pie <-\"$opts{pie}\"
horiz <-\"$opts{horiz}\"
################# plot  bar     ########################################
if(bar ==\"T\"){
	dat <-sapply(1:ncol(otu),function(x) otu[,x]/sum(otu[,x]))
	colnames(dat) <-colnames(otu)
	rownames(dat) <-rownames(otu)
	dat1 <-sapply(1:ncol(otu),function(x) otu[,x])
	colnames(dat1) <-colnames(otu)
	rownames(dat1) <-rownames(otu)
        dat1 <- as.data.frame(dat1)
        dat1\$ID <- rownames(dat1)
        dat1 <- dat1[,c(\"ID\",colnames(dat1)[1:(ncol(dat1)-1)])]
        #write.table(dat1,\"$opts{g}.new.$opts{i}\",sep=\"\\t\",eol=\"\\n\",quote=FALSE)
        write.table(dat1,\"all.bar\",sep=\"\\t\",eol=\"\\n\",quote=FALSE,row.names=F)
        mycol<-c(\"#EE3B3B\",\"#458100\",\"#FFD700\",\"#0000FF\",\"#8A2BE2\",\"#8B3A62\",\"#912323\",\"#898380\",\"#458F74\",\"#00008B\",\"#8B008B\",\"#A2CD5A\",\"#FF147F\",\"#00BFFF\",\"#EEAD0E\",\"#8B4F13\",\"#89864E\",\"#FF6A6A\",\"#8B5F65\",\"#8470FF\",\"#B0C4DE\",\"#5D478B\",\"#8B4789\",\"#1C86EE\",\"#D2691E\",\"#7FFF00\",\"#CDB6B0\",\"#66CDAA\",\"#98F5FF\",\"#EE7611\",\"#FF69B4\",\"#F0E68C\",\"#CD8C95\",\"#FFF9DB\",\"#7FFFD4\",\"#CD5555\",\"#009ACD\",\"#32CD32\",\"#FFB6C1\",\"#008B45\",\"#B2DFF8\",\"#68858B\",\"#00CD66\",\"#8B8686\",\"#FFFF00\",\"#CDB5CD\",\"#CD4F39\",\"#9ACD32\",\"#00C5CD\",\"#CDCD00\")
        mycol <- rep(mycol,50)
	if (\"$opts{color}\" == \"FALSE\") {mycol<-mycol[1:(dim(dat1)[1])]}
        else {mycol <- readLines(\"$opts{color}\")}
        
	lab <-rownames(dat)
	str <-strwidth(lab, units = \"inches\")
	cnum <-ceiling(length(str)/$opts{b_ncol})
	nstr <- lapply(1:$opts{b_ncol},function(x) if(cnum*x<=length(str)){str[(cnum*(x-1)+1):(cnum*x)]}else{str[(cnum*(x-1)+1):length(str)]})
	nlab <- lapply(1:$opts{b_ncol},function(x) if(cnum*x<=length(lab)){lab[(cnum*(x-1)+1):(cnum*x)]}else{lab[(cnum*(x-1)+1):length(lab)]})
	mw <-unlist(lapply(nstr,max))
	mw[1] <-mw[1]*$opts{b_l1}

	lx <-\"$opts{b_lx}\"
	tiff(file=\"bar.$gs.$opts{i}.tiff\",width=$opts{b_w},height=$opts{b_h},pointsize=15)
	
	#if(ncol(dat)>3){pdf(file=\"bar.$gs.$opts{i}.pdf\",width=$opts{b_w},height=$opts{b_h})}
	#else{pdf(file=\"bar.$gs.$opts{i}.pdf\",width=$opts{b_w_w},height=$opts{b_h_w})}
	
	if(ncol(dat)>3){pdf(file=\"bar.pdf\",width=$opts{b_w},height=$opts{b_h})}
	else{pdf(file=\"bar.pdf\",width=$opts{b_w_w},height=$opts{b_h_w})}
	if(lx==\"T\"){
	  la <-layout(matrix(c(rep(1,$opts{b_ncol}),2:($opts{b_ncol}+1)),2,$opts{b_ncol},byrow=TRUE),width=mw,heights=c($opts{b_hi}))
	}else{
	  layout(matrix(1:2,2,1),heights=c($opts{b_hi}))
	}
	#par(mar=c(3,5,2,2))
	las <-1
	if(ncol(dat)>15) las <-3

	blas <-$opts{b_blas}
	if(blas >0) las <-blas
	if(horiz ==\"T\"){
	par(mar=c(5,5,2,2))
	barplot(dat*100,horiz=T,width=1,space=1.2,plot=T,las=las,col=mycol[1:nrow(dat)],cex.axis=$opts{b_scex},cex.names=$opts{b_scex},border=NA,xlab=\"Relative abundance (%)\",offset=0,cex.lab=$opts{b_scex})
	}else{
	if(ncol(dat)>3){
        par(mar=c($opts{b_bo},5,2,2))
	barplot(dat*100,width=$opts{width},space=$opts{space},plot=T,las=las,col=mycol[1:nrow(dat)],cex.axis=$opts{b_scex},cex.names=$opts{b_scex},border=NA,ylab=\"Relative abundance (%)\",offset=0,cex.lab=$opts{b_scex})}
	#abline(h=0,lwd=1.2)
	else{
	layout(matrix(1:2,1,2),widths=c($opts{b_wi}))
	par(mar=c(5,5,2,0))
	barplot(dat*100,width=$opts{width},space=$opts{space},plot=T,las=las,col=mycol[1:nrow(dat)],cex.axis=$opts{b_scex},cex.names=$opts{b_scex},border=NA,ylab=\"Relative abundance (%)\",offset=0,cex.lab=$opts{b_scex})}
	
	}

	if(lx==\"T\"){
	 for(i in 1:$opts{b_ncol}){
	   if(i==1){par(mar=c(1,5,1,0))}
	   else {par(mar=c(1,0,1,0))}
	   plot.new()
	   legend(\"topleft\",legend=nlab[[i]],fill=mycol[(cnum*(i-1)+1):(cnum*i)],cex=$opts{b_lcex},bty=\"n\")
	 }
	}else{
	if(ncol(dat)>3){
	 par(mar=c(2,5,1,1))
	 plot.new()
	 legend(\"topleft\",legend=rownames(dat),ncol=$opts{b_ncol},fill=mycol[1:nrow(dat)],cex=$opts{b_lcex},bty=\"n\")
	}
	else{par(mar=c(2,3,2,1))
         plot.new()
	legend(\"topleft\",legend=rownames(dat),ncol=$opts{b_ncol_w},fill=mycol[1:nrow(dat)],cex=$opts{b_lcex},bty=\"n\")
	}
	dev.off()
	 }
	if(!is.null(del)){
		prt <-paste(\"Warning: Samples:\",mydel,\"were none,disregard them in the plot.\")
	}
#	rm(dat)


######## bar finished ####################
###png
	if(ncol(dat)>3){png(file=\"bar.png\",width=$opts{b_w}*100,height=$opts{b_h}*100,res=300)}
	else{png(file=\"bar.png\",width=$opts{b_w_w}*100,height=$opts{b_h_w}*100,res=300)}
	if(lx==\"T\"){
	  la <-layout(matrix(c(rep(1,$opts{b_ncol}),2:($opts{b_ncol}+1)),2,$opts{b_ncol},byrow=TRUE),width=mw,heights=c($opts{b_hi}))
	}else{
	  layout(matrix(1:2,2,1),heights=c($opts{b_hi}))
	}
	#par(mar=c(3,5,2,2))
	las <-1
	if(ncol(dat)>15) las <-3

	blas <-$opts{b_blas}
	if(blas >0) las <-blas
	if(horiz ==\"T\"){
	par(mar=c(5,5,2,2))
	barplot(dat*100,horiz=T,width=1,space=1.2,plot=T,las=las,col=mycol[1:nrow(dat)],cex.axis=$opts{b_scex},cex.names=$opts{b_scex},border=NA,xlab=\"Relative abundance (%)\",offset=0,cex.lab=$opts{b_scex})
	}else{
	if(ncol(dat)>3){
        par(mar=c($opts{b_bo},5,2,2))
	barplot(dat*100,width=$opts{width},space=$opts{space},plot=T,las=las,col=mycol[1:nrow(dat)],cex.axis=$opts{b_scex},cex.names=$opts{b_scex},border=NA,ylab=\"Relative abundance (%)\",offset=0,cex.lab=$opts{b_scex})}
	#abline(h=0,lwd=1.2)
	else{
	layout(matrix(1:2,1,2),widths=c($opts{b_wi}))
	par(mar=c(5,5,2,0))
	barplot(dat*100,width=$opts{width},space=$opts{space},plot=T,las=las,col=mycol[1:nrow(dat)],cex.axis=$opts{b_scex},cex.names=$opts{b_scex},border=NA,ylab=\"Relative abundance (%)\",offset=0,cex.lab=$opts{b_scex})}
	
	}

	if(lx==\"T\"){
	 for(i in 1:$opts{b_ncol}){
	   if(i==1){par(mar=c(1,5,1,0))}
	   else {par(mar=c(1,0,1,0))}
	   plot.new()
	   legend(\"topleft\",legend=nlab[[i]],fill=mycol[(cnum*(i-1)+1):(cnum*i)],cex=$opts{b_lcex},bty=\"n\")
	 }
	}else{
	if(ncol(dat)>3){
	 par(mar=c(2,5,1,1))
	 plot.new()
	 legend(\"topleft\",legend=rownames(dat),ncol=$opts{b_ncol},fill=mycol[1:nrow(dat)],cex=$opts{b_lcex},bty=\"n\")
	}
	else{par(mar=c(2,3,2,1))
         plot.new()
	legend(\"topleft\",legend=rownames(dat),ncol=$opts{b_ncol_w},fill=mycol[1:nrow(dat)],cex=$opts{b_lcex},bty=\"n\")
	}
	dev.off()
	 }
	if(!is.null(del)){
		prt <-paste(\"Warning: Samples:\",mydel,\"were none,disregard them in the plot.\")
	}
#	rm(dat)


######## bar finished ####################
#svg
	if(ncol(dat)>3){svg(file=\"bar.svg\",width=$opts{b_w},height=$opts{b_h})}
	else{svg(file=\"bar.svg\",width=$opts{b_w_w},height=$opts{b_h_w})}
	if(lx==\"T\"){
	  la <-layout(matrix(c(rep(1,$opts{b_ncol}),2:($opts{b_ncol}+1)),2,$opts{b_ncol},byrow=TRUE),width=mw,heights=c($opts{b_hi}))
	}else{
	  layout(matrix(1:2,2,1),heights=c($opts{b_hi}))
	}
	#par(mar=c(3,5,2,2))
	las <-1
	if(ncol(dat)>15) las <-3

	blas <-$opts{b_blas}
	if(blas >0) las <-blas
	if(horiz ==\"T\"){
	par(mar=c(5,5,2,2))
	barplot(dat*100,horiz=T,width=1,space=1.2,plot=T,las=las,col=mycol[1:nrow(dat)],cex.axis=$opts{b_scex},cex.names=$opts{b_scex},border=NA,xlab=\"Relative abundance (%)\",offset=0,cex.lab=$opts{b_scex})
	}else{
	if(ncol(dat)>3){
        par(mar=c($opts{b_bo},5,2,2))
	barplot(dat*100,width=$opts{width},space=$opts{space},plot=T,las=las,col=mycol[1:nrow(dat)],cex.axis=$opts{b_scex},cex.names=$opts{b_scex},border=NA,ylab=\"Relative abundance (%)\",offset=0,cex.lab=$opts{b_scex})}
	#abline(h=0,lwd=1.2)
	else{
	layout(matrix(1:2,1,2),widths=c($opts{b_wi}))
	par(mar=c(5,5,2,0))
	barplot(dat*100,width=$opts{width},space=$opts{space},plot=T,las=las,col=mycol[1:nrow(dat)],cex.axis=$opts{b_scex},cex.names=$opts{b_scex},border=NA,ylab=\"Relative abundance (%)\",offset=0,cex.lab=$opts{b_scex})}
	
	}

	if(lx==\"T\"){
	 for(i in 1:$opts{b_ncol}){
	   if(i==1){par(mar=c(1,5,1,0))}
	   else {par(mar=c(1,0,1,0))}
	   plot.new()
	   legend(\"topleft\",legend=nlab[[i]],fill=mycol[(cnum*(i-1)+1):(cnum*i)],cex=$opts{b_lcex},bty=\"n\")
	 }
	}else{
	if(ncol(dat)>3){
	 par(mar=c(2,5,1,1))
	 plot.new()
	 legend(\"topleft\",legend=rownames(dat),ncol=$opts{b_ncol},fill=mycol[1:nrow(dat)],cex=$opts{b_lcex},bty=\"n\")
	}
	else{par(mar=c(2,3,2,1))
         plot.new()
	legend(\"topleft\",legend=rownames(dat),ncol=$opts{b_ncol_w},fill=mycol[1:nrow(dat)],cex=$opts{b_lcex},bty=\"n\")
	}
	dev.off()
	 }
	if(!is.null(del)){
		prt <-paste(\"Warning: Samples:\",mydel,\"were none,disregard them in the plot.\")
	}
#	rm(dat)
}

######## bar finished ####################
if(pie ==\"T\"){

######## function spie ##########
spie <-function (x, labels = names(x), edges = 200, radius = 0.8, clockwise = FALSE, 
    init.angle = if (clockwise) 90 else 0, density = NULL, angle = 45, 
    col = NULL, border = NULL, lty = NULL, main = NULL, textlab=TRUE, legend = FALSE,...) 
{   
    if (!is.numeric(x) || any(is.na(x) | x < 0)) 
        stop(\"\'x\' values must be positive.\")
    if (is.null(labels)) 
        labels <- as.character(seq_along(x))
    else labels <- as.graphicsAnnot(labels)

    if(legend) { layout(matrix(1:2,1,2) ) 
                 par1 <-par(mar=c(1,1,1,0))
                 par2 <-par(mar=c(1,0,4,1))
    }

    x <- c(0, cumsum(x)/sum(x))
    dx <- diff(x)
    nx <- length(dx)
    plot.new()
    pin <- par(\"pin\")
    xlim <- ylim <- c(-1, 1)
    if (pin[1L] > pin[2L]) 
        xlim <- (pin[1L]/pin[2L]) * xlim
    else ylim <- (pin[2L]/pin[1L]) * ylim
    vx<-$opts{p_xlim}+0.3
	vy<-$opts{p_ylim}+0.3
	ylim[2] <-ylim[2]*vy	
    if(legend) {v <- 0.8;par(par1)}
    plot.window(xlim*vx, ylim, \"\", asp = 1)
    if (is.null(col)) 
        col <- if (is.null(density)) 
            c(\"white\", \"lightblue\", \"mistyrose\", \"lightcyan\", 
                \"lavender\", \"cornsilk\")
        else par(\"fg\")
    col <- rep(col, length.out = nx)
    border <- rep(border, length.out = nx)
    lty <- rep(lty, length.out = nx)
    angle <- rep(angle, length.out = nx)
    density <- rep(density, length.out = nx)
    twopi <- if (clockwise) 
        -2 * pi
    else 2 * pi
    t2xy <- function(t) {
        t2p <- twopi * t + init.angle * pi/180
        list(x = radius * cos(t2p), y = radius * sin(t2p))
    }

    pl1 <- pl2 <- c(0,-radius*1.5)
    lb1 <- lb2 <-array()
    lb1[1] <- lb2[1] <- \"lab0\"
    li1 <- li2 <-1
    for (i in 1L:nx) {
        n <- max(2, floor(edges * dx[i]))
        P <- t2xy(seq.int(x[i], x[i + 1], length.out = n))
        polygon(c(P\$x, 0), c(P\$y, 0), density = density[i], angle = angle[i], 
            border = border[i], col = col[i], lty = lty[i])
        pm <- t2xy(mean(x[i + 0:1]))
        lab <- as.character(labels[i])
        if (!is.na(lab) && nzchar(lab)) {
            if(pm\$x >=0){ 
                 l1 <-c(pm\$x,pm\$y)
                 pl1 <-rbind(pl1,l1)
                 li1 <-li1+1
                 lb1[li1] <-lab
            }else { 
                 l2 <-c(pm\$x,pm\$y)
                 pl2 <-rbind(pl2,l2)
                 li2 <-li2+1
                 lb2[li2] <-lab
            }
        }
    }
    row.names(pl1) <-lb1
    row.names(pl2) <-lb2


############### labelxy ####
    labelxy <- function(pl) {
         pl <-pl[order(pl[,2]),]
         by <-pl[1,2]
         bx <-pl[1,1]
         d1 <-1.3
         for(j in 2:nrow(pl)) {
             ply <- pl[j,2]
             plx <- pl[j,1]
             pmy <- 1.3 * ply
             pmx <- 1.3 * plx
             d2 <-1.3
             if(ply<0 & abs(plx)<0.3*radius){ 
                  if(d1>1.1) {d1 <-d1-0.04}    
					# while(abs(pmx-bx)<0.1 ){                   
                     # d2 <-d2+0.01
                     # pmx <- plx*d2}				  
					if(abs(pmx-bx)<0.1){
							if(pmx-bx>0){
									pmx=0.1+bx
							}else{
									pmx=-0.1+bx
							}
							d2<-pmx/plx
					}
             }else{ d1 <-1.2 }
			 if(pmy-by<0.08) {pmy <-0.08+by}
             #while(pmy-by<0.08){ pmy <- pmy+0.01 }
             d3 <-d2
			 if(pmy*pmy+pmx*pmx <=radius*radius){
				 if(pmx>0){
						pmx <- sqrt(radius*radius-pmy*pmy)
				 }else{
						pmx <- -sqrt(radius*radius-pmy*pmy)
				 }
				 d3<-pmx/plx
			 }
             # while(pmy*pmy+pmx*pmx <=radius*radius){
                  # d3 <-d3+0.01
                  # pmx <- plx*d3
             # } 
             
             lines(c(1, d1) * plx,c(1, d1) * ply) #line1                                                  
             lines(c(d1* plx, pmx) , c(d1*ply, pmy) ) #line2
             lines(c(pmx, pmx*1.03) , c(pmy,pmy ))  #line3
             text(1.036*pmx, pmy, rownames(pl)[j], xpd = TRUE, adj = ifelse(plx < 0, 1, 0), ...)
             by <-pmy
             bx <-pmx
         }
    }
###############

  if(textlab) {       
            labelxy(pl1)
            labelxy(pl2)
  }
  title(main = main, ...)
  if(legend) { plot.new()
               par(par2)
               legend(\"topleft\",legend=labels,fill=col)

  }


  invisible(NULL)


}

#####################################################################

########################################################################################pdf start
##### function ppie : plot a table of samples #########
ppie <-function(dat,label,col,smp){
	#tiff(paste(\"pie.\",smp,\".$opts{i}.tiff\",sep=\"\"),width=$opts{p_w},height=$opts{p_h},pointsize=15)
	#pdf(paste(\"pie.\",smp,\".$opts{i}.pdf\",sep=\"\"),width=$opts{p_w},height=$opts{p_h})
	pdf(paste(\"pie.\",smp,\".pdf\",sep=\"\"),width=$opts{p_w},height=$opts{p_h})
	label <-sapply(1:nrow(dat),function(x) paste(label[x],\" \",round(dat[x,1]/sum(dat)*100,digits=2),\"%\",sep=\"\"))
	spie(dat,border=NULL,labels=label,col=col,main=smp,cex=$opts{p_cex},textlab=$labtext,legend=$legend,clockwise=$opts{clockwise},init.angle =$opts{angle})
	dev.off()
}
pcol=mycol[1:nrow(otu)]
dl <-lapply(1:ncol(otu),function(y) which(otu[,y]>0))
sapply(1:ncol(otu),function(y) ppie(dat=as.matrix(otu[dl[[y]],y]),label=rownames(otu)[dl[[y]]],col=pcol[dl[[y]]],smp=colnames(otu)[y]))




########################################################################################png start 
ppie <-function(dat,label,col,smp){
	#tiff(paste(\"pie.\",smp,\".$opts{i}.tiff\",sep=\"\"),width=$opts{p_w},height=$opts{p_h},pointsize=15)
	#pdf(paste(\"pie.\",smp,\".$opts{i}.pdf\",sep=\"\"),width=$opts{p_w},height=$opts{p_h})
#png
	png(paste(\"pie.\",smp,\".png\",sep=\"\"),width=$opts{p_w}*100*.8,height=$opts{p_h}*100*.8,res=300)
	label <-sapply(1:nrow(dat),function(x) paste(label[x],\" \",round(dat[x,1]/sum(dat)*100,digits=2),\"%\",sep=\"\"))
	spie(dat,border=NULL,labels=label,col=col,main=smp,cex=$opts{p_cex},textlab=$labtext,legend=$legend,clockwise=$opts{clockwise},init.angle =$opts{angle})
	dev.off()
}
pcol=mycol[1:nrow(otu)]
dl <-lapply(1:ncol(otu),function(y) which(otu[,y]>0))
sapply(1:ncol(otu),function(y) ppie(dat=as.matrix(otu[dl[[y]],y]),label=rownames(otu)[dl[[y]]],col=pcol[dl[[y]]],smp=colnames(otu)[y]))


########################################################################################svg start 
ppie <-function(dat,label,col,smp){
	#tiff(paste(\"pie.\",smp,\".$opts{i}.tiff\",sep=\"\"),width=$opts{p_w},height=$opts{p_h},pointsize=15)
	#pdf(paste(\"pie.\",smp,\".$opts{i}.pdf\",sep=\"\"),width=$opts{p_w},height=$opts{p_h})

#svg
	svg(paste(\"pie.\",smp,\".svg\",sep=\"\"),width=$opts{p_w},height=$opts{p_h})
	label <-sapply(1:nrow(dat),function(x) paste(label[x],\" \",round(dat[x,1]/sum(dat)*100,digits=2),\"%\",sep=\"\"))
	spie(dat,border=NULL,labels=label,col=col,main=smp,cex=$opts{p_cex},textlab=$labtext,legend=$legend,clockwise=$opts{clockwise},init.angle =$opts{angle})
	dev.off()	
}
pcol=mycol[1:nrow(otu)]
dl <-lapply(1:ncol(otu),function(y) which(otu[,y]>0))
sapply(1:ncol(otu),function(y) ppie(dat=as.matrix(otu[dl[[y]],y]),label=rownames(otu)[dl[[y]]],col=pcol[dl[[y]]],smp=colnames(otu)[y]))
}






";


`R --restore --no-save < cmd.r`;
#system ('rm cmd.r');
system ('rm Rplots.pdf') if -e "Rplots.pdf";

}









