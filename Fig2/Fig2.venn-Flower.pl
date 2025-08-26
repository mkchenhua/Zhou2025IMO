#!/usr/bin/perl -w
###################################################################
#####Author : GuoqiLiu                                             #
#####Date   : 2024-04-17                                           #
#####Copyright (C) 2019~ Mingke Biotechnology (Hangzhou) Co., Ltd. #
#####Contact: liuguoqi@mingkebio.com                               #
#####Suppose: Auto diversity report  program                       #
#####step :                                                        #
#####1. diver_pipe_V3.pl generate basic report                     #
#####2. modify html and data                                       #
#####3. html2pdf report                                            #
#####4. Auto-send report ? whether or not                          #
#####Platform :                                                    #
##############Ubuntu 18.04 & Windows10 ,Perl v5.26 #################
####################################################################
#Log:
##V1 version 2024-06-12 

##
use strict;
use warnings;
use Getopt::Long;
use POSIX; 

my %opts;
my $v = "v2024-06-24";

my $help;

GetOptions (\%opts,"i=s","o=s","sample=s","l=s","type=i","sz=f","alpha=f","cex=f",,"c=i","method=s","std=s","w=f","h=f","m=s","g=s","color=s","rlc=f","point=s","clc=f","llc=f","lg=i","angle=i","test=s");#,"scale=s");
my $usage = <<"USAGE";
#######################################################################################################################################################
#                 Program : perl $0 
#                 Version :  $v
#######################################################################################################################################################
#                -i*  	   <str>             matrix file  
#                -m        <str>             map file
#                -method   <str>             group calculate method sum,mean default:mean
#                -g        <str>             sample or group sort list
#                -color    <str>             sample or group color 
#                -type     <int>             can choice 1 or 2,1 mean border 2 mean fill
#                -cex     <float>            point label size default:1
#                -sample   <str>             T/F, T display sample number label ,F not display sample label,default:F  
#		 -w       <float>            the width of the figure,default:6
#		 -h	  <float>            the height of the figure,default:6.5
#
#                Example :                Usage:perl $0  -i otu/ASV_table.xls  
#########################################################################################################################################################

USAGE
#die $usage if (!($opts{i}&&$opts{g}));
##                -lw    split plot in width,defalt (four numbers) :0.1:0.2:4:1
#                -lh    split plot in heigth,defalt (three numbers) :0.3:5.5:1.2


#die $usage if ( !(defined $opts{i} && $opts{m} && $opts{rep}) );

die $usage if (!($opts{i}));
#die $usage if (!($opts{m}));
#die $usage if (!($opts{rep}));
#




#####                -llc     <float>            legend size
#die $usage if (!($opts{i2}));
#$opts{o}=defined$opts{o}?$opts{o}:"output";
$opts{w} ||= 6;
$opts{h} ||= 6.5;
$opts{m} ||= "F";
$opts{g} ||= "F";
$opts{sample} ||= "F";
$opts{color} ||= "F" ;
$opts{method} ||= "mean";
$opts{cex} ||= 1;
$opts{type} ||= 1;

# print "$pathqiime2\n"; 
# 
sub get_time {
	 my $hhh=shift;
	 my $gettime=strftime("%Y-%m-%d %H:%M:%S",localtime());
	 print  "\@\@$gettime====>>>>>>>$hhh\n";
}
&get_time("$0 program start running...") ;




#print "$opts{angle} ------------------\n";

########################################################

sub calculate_more5sam_core{
	my ($input,$sort_sam) = @_;

	my @samples = split /\t/,$sort_sam;


open II,$input;
chomp(my $head = <II>) ; 
my @hh = split/\t/,$head; 
my @array = () ;
my $t=0;
my %sam2index=();
foreach my $jj (@hh[1..$#hh]) {
	$t++;
	$sam2index{$jj} = $t; 
}
foreach my $i (@samples) {
	push @array,$sam2index{$i} ; 
}
open TT,">tmp.sort.txt";
print TT "OTU\t",join("\t",@samples),"\n";
while (<II>) {
	chomp;
	my @tmp=split/\t/,$_;
	my @ttt;
	push @ttt,$tmp[0];
	foreach my $iii (@array) {
		push @ttt,$tmp[$iii] ;}
	print TT join("\t",@ttt),"\n";
}

close TT;close II;














open A,"tmp.sort.txt";

my $h=<A>;	chomp $h;	my @h=split /\t+/,$h;
open O,">venn.sets.xls";
print O "Sam\tCore_otu\tspecial_out\tall_otu\n";
my $core=0;	my %rec;

my %core_hash=();
my %uniq_hash=(); 
open OO,">venn2dat.tsv";
my %s2index=();

foreach my $k(1..$#h){
	$rec{$k}[0]=0;	$rec{$k}[1]=0;
	$s2index{$k} = $h[$k];
}
while(<A>){
	chomp;my @a=split /\t+/;
	my $num=0;
	foreach my $k(1..$#a){
		if($a[$k] != 0){
			$num++;	$rec{$k}[0]++;
			

		}
	}
	if($num == $#a){
		$core++;
		push @{$core_hash{"all"}},$a[0]; 
	}
	if($num == 1){
		foreach my $k(1..$#a){
			if($a[$k] != 0){
				$rec{$k}[1]++;
				push @{$uniq_hash{$s2index{$k}}},$a[0];
			}
		}
	}
}
close A;

foreach my $k(1..$#h){
	print O "$h[$k]\t$core\t$rec{$k}[1]\t$rec{$k}[0]\n";
	#print O "$tmp[$k]\t$core\t$rec{$k}[1]\t$rec{$k}[0]\n";
}
#print OO "All_samples_exists\t",scalar @{$core_hash{"all"}},"\t",join("\t",@{$core_hash{"all"}}),"\n"; 

if (exists $core_hash{"all"}) {
print OO "All_samples_exists\t",join("\t",@{$core_hash{"all"}}),"\n";
}
else {
	print OO "All_samples_exists\t","NA","\n";
}
foreach my $k (keys %uniq_hash){
	print OO "$k\t",join("\t",@{$uniq_hash{$k}}),"\n";
        #	print OO "$k\t",scalar @{$uniq_hash{$k}},"\t",join("\t",@{$uniq_hash{$k}}),"\n";
}


`rm tmp.sort.txt`;
}

# &calculate_more5sam_core ($opts{i}) ;


sub calculate_less_and_eq5sam_core{
	my ($input,$sort_sam,$sort_color) = @_; 
open Venn,">cmd1.r";
print Venn "
lst <-list()
table <-read.table(file=\"$input\",sep=\"\\t\",head=T,check.names = FALSE)
rownames(table) <- as.character(table[,1])
table <-table[,-1]
table <-table[,c($sort_sam)]

for(i in 1:length(colnames(table))){
       samp <-colnames(table)[i]
       lst[[samp]] <- rownames(table)[which(table[,i] !=0)]
}



intersect_n <-function(x) {
        n <-length(x)
        x2 <- x[[1]]
        for(i in 2:n){
               x2 <- intersect(x[[i]],x2)
        } 
        return(x2)        
}

getset <- function(lst){   #####lst is a list
    sn <-length(names(lst))
    sets <- list()
    sls <- list()
    maxl <-0
    #### get all intersect ####  
    for(i in sn:2){
             sl <- combn (names(lst),i,simplify=FALSE)
             inter <-lapply(sl,function(x) intersect_n(lst[x]))
             names(inter) <-lapply(sl,function(x) paste(x,collapse = \"__\"))
             sets <- c(sets,inter)
             sls <- c(sls,sl)
                          
    }
    #### get all unique ####
    for(i in 1:sn){
             uniq <- list(setdiff(unlist(lst[names(lst)[i]]),unlist(lst[names(lst)[-i]])))
             names(uniq) <- paste(names(lst)[i],\"__only\",sep=\"\")
             sets <- c(sets,uniq)
             
    } 
    return(sets)

}
sets <- getset(lst)

###### write sets to file ######
maxl <-max(sapply(sets,length))
sets_table <- sapply(sets,function(x) c(x,rep(\"\",maxl-length(x))))

otuset <-\"venn.sets.xls\"

a <- data.frame(sets_table)
write.table(a,otuset,sep = \"\\t\",eol=\"\\n\",row.names=FALSE,quote=FALSE)   


####dram venn diagram

library(\"VennDiagram\")
mycol <- c(\"red\",\"orange\",\"green\",\"blue\",\"yellow\")
	#venn.plot<-venn.diagram(lst,filename = NULL,col=mycol[1:3],fill=mycol[1:3],print.mode = c(\"percent\", \"raw\"),ext.pos=-180,cat.cex=.7,scaled = F)
	#pdf(\"venn5group.pdf\")
        #grid.draw(venn.plot)
        #dev.off()
if ($opts{type} == 1) {
	venn.plot<-venn.diagram(lst,filename = NULL,col=c($sort_color),cex=$opts{cex},cat.cex=$opts{cex})
	}else {venn.plot<-venn.diagram(lst,filename = NULL,fill=c($sort_color),cex=$opts{cex},cat.cex=$opts{cex})}
pdf(file = \"venn.pdf\",width = $opts{w},height=$opts{h})
#print(p2)
grid.draw(venn.plot)
dev.off()
png(file = \"venn.png\",width = $opts{w}*240,height=$opts{h}*200,res=300)
#print(p2)
grid.draw(venn.plot)
dev.off()
svg(file = \"venn.svg\",width = $opts{w},height=$opts{h})
#print(p2)
grid.draw(venn.plot)
dev.off()




x <- lst
inter <- get.venn.partitions(x)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = '\t')
inter\$sample<-\"a\"
for (i in 1:nrow(inter)) inter\$sample[i] <- paste(names(unlist(inter[,1:(ncol(inter)-5)][i,]))[unlist(inter[,1:(ncol(inter)-5)][i,])],collapse=\"__\") 
#colnames(head(inter[,c(ncol(inter),ncol(inter)-2,ncol(inter)-1)]))
write.table(inter[,c(ncol(inter),ncol(inter)-1)], 'venn2dat.tsv', row.names = FALSE, col.names = FALSE,sep =  \"\\t\", quote = FALSE)




";
`R --no-save<cmd1.r`;
}

# &calculate_more5sam_core ($opts{i}) ;
# &calculate_less_and_eq5sam_core($opts{i}) ;

sub calculate_group{
	        my ($input,$map,$method) = @_;
		open Group_stat,">cmd2.r";
print Group_stat "
otu <-read.table(file=\"$input\",sep=\"\\t\",head=T,check.names = FALSE,row.names=1)
group <- read.table(\"$map\",skip=1,check.names = FALSE)
glst <- lapply(1:length(unique(group[,2])),function(x)group[which(group[,2] %in% unique(group[,2])[x]),1])
names(glst) <-as.character(unique(group[,2]))
tab <-sapply(1:length(glst),function(x) apply(otu[as.character(as.vector(glst[[x]]))],1,$method))
otu <-tab[apply(tab,1,function(x)any(x>0)),]      
colnames(otu) <-unique(group[,2])
data <- as.data.frame(otu)
data\$OTU <- rownames(data)
data <- data[,c(\"OTU\",colnames(data)[1:(ncol(data)-1)])]
write.table(data,file=paste(\"tmp.\",\"$method\",\".txt\",sep=\"\"),sep=\"\\t\",quote=F,row.names=F)
";
`R --no-save<cmd2.r`;
}
#&calculate_group($opts{i},$opts{m},$opts{method})

sub More5sam {
	my ($input) = @_;

        open R5venn,">cmd3.r";
	print R5venn "
	library(\"VennDiagram\")
	mycol <- c(\"red\",\"orange\",\"green\",\"blue\",\"yellow\")
	venn.plot<-venn.diagram(lst,filename = NULL,col=mycol[1:3],fill=mycol[1:3],print.mode = c(\"percent\", \"raw\"),ext.pos=-180,cat.cex=.7,scaled = F)
	pdf(\"venn5group.pdf\")
	grid.draw(venn.plot)
	dev.off()
	";
	#`R --no-save<cmd3.r`;
}




sub flower{
my ($input,$sort_color) = @_ ;
open Rflower,">cmd3.r"; 
print Rflower "

theta= seq(0,2*pi,pi/1000)
inf=read.table(\"$input\",header=T,check.names = FALSE)
pdf(file=\"flower.pdf\",width = $opts{w},height=$opts{h})
test=inf
aa=nrow(test) 
bb=ncol(test)
r=1;a=4;b=1;
n=4
require(graphics)
library(scales)
mycol <- c(34, 51, 142, 26, 31, 371, 36, 7, 12, 30, 84, 88, 116, 121, 77, 56, 386, 373, 423, 435, 438, 471, 512, 130, 52, 47, 6, 11, 43, 54, 367, 382, 422, 4, 8, 375, 124, 448, 419, 614, 401, 403, 613, 583, 652, 628, 633, 496, 638, 655, 132, 503, 24)
mycol <-colors()[rep(mycol,20)]
color = mycol
color <- c($sort_color)
row=nrow(color)
column=ncol(color)
space=seq(0,2*pi,length=aa+1);
rota=seq(0,360,length=aa+1);
for (i in 1:aa){
x0=cos(space[i])*(a-r);
y0=sin(space[i])*(a-r);
alpha=space[i];
#text_x1=(2*a-1+1)*cos(space[i]);
#text_y1=(2*a-1+1)*sin(space[i]);
text_x1=(2*a-1+1+.9)*cos(space[i]);
text_y1=(2*a-1+1+.9)*sin(space[i]);
#text_x2=(2*a-2)*cos(space[i]);
#text_y2=(2*a-2)*sin(space[i]);
text_x2=(2*a-1.6)*cos(space[i]);
text_y2=(2*a-1.6)*sin(space[i]);
text_xx1=(2*a-1+0.5)*cos(space[i]);
	text_yy1=(2*a-1+0.5)*sin(space[i]);
	text_xx2=(2*a-1.5)*cos(space[i]);
	text_yy2=(2*a-1.5)*sin(space[i]);
	marker1=test[i,1];
	marker2=test[i,3];
	x3=x0+a*cos(alpha)*cos(theta)-b*sin(alpha)*sin(theta);
	y3=y0+a*sin(alpha)*cos(theta)+b*cos(alpha)*sin(theta);
	par(new=T)

	plot(x3,y3,type=\"l\",xlim=c(-8,8),ylim=c(-8,8),lwd=2,axes=F,xlab=\"\",ylab=\"\",col=c($sort_color)[i])
	#polygon(x3,y3,col = rgb(255, 0, 0, 20, maxColorValue=255))
	if ($opts{type} == 2){
	polygon(x3,y3,col = alpha(c($sort_color)[i],0.3),border=NA)}
	
	name=paste(\"n=\",test[i,4])
	#		if (rota[i] <= 180){
	#	text(text_x1,text_y1,labels=marker1,srt=(rota[i]-90),xpd=NA,cex=$opts{cex},font=0.5)
	#	text(text_xx1,text_yy1,labels=name,srt=(rota[i]-90),xpd=NA,cex=$opts{cex},font=0.5)
	#	text(text_x2,text_y2,labels=marker2,srt=(rota[i]-90),xpd=NA,cex=$opts{cex},font=0.5);
	#	}else{
	#	text(text_x1,text_y1,labels=marker1,srt=(rota[i]+90),xpd=NA,cex=$opts{cex},font=0.5)
	#	text(text_xx1,text_yy1,labels=name,srt=(rota[i]+90),xpd=NA,cex=$opts{cex},font=0.5)
	#	text(text_x2,text_y2,labels=marker2,srt=(rota[i]+90),xpd=NA,cex=$opts{cex},font=0.5);
	#		}
	#	}

	if (rota[i] < 90) {
	        if(\"$opts{sample}\" == \"T\"){
		text(text_x1,text_y1,labels=paste(marker1,\"(\",name,\")\",sep=\"\"),srt=(rota[i]),xpd=NA,cex=$opts{cex},font=0.5)
		}else {text(text_x1,text_y1,labels=marker1,srt=(rota[i]),xpd=NA,cex=$opts{cex},font=0.5)}
		#text(text_xx1,text_yy1,labels=name,srt=(rota[i]),xpd=NA,cex=$opts{cex},font=0.5)
		text(text_x2,text_y2,labels=marker2,srt=(rota[i]-90),xpd=NA,cex=$opts{cex},font=0.5);
		}else if ((rota[i] >= 135) && (rota[i] <= 270)){
		if(\"$opts{sample}\" == \"T\"){
		text(text_x1,text_y1,labels=paste(marker1,\"(\",name,\")\",sep=\"\"),srt=(rota[i]+180),xpd=NA,cex=$opts{cex},font=0.5)
		}else {text(text_x1,text_y1,labels=marker1,srt=(rota[i]+180),xpd=NA,cex=$opts{cex},font=0.5)}
		#text(text_xx1,text_yy1,labels=name,srt=(rota[i]+90+90),xpd=NA,cex=$opts{cex},font=0.5)
		text(text_x2,text_y2,labels=marker2,srt=(rota[i]+90),xpd=NA,cex=$opts{cex},font=0.5);
	                    }else if ((rota[i] >= 90) && (rota[i] < 135)){
			    if(\"$opts{sample}\" == \"T\"){
		text(text_x1,text_y1,labels=paste(marker1,\"(\",name,\")\",sep=\"\"),srt=(rota[i]-180),xpd=NA,cex=$opts{cex},font=0.5)
		}else {text(text_x1,text_y1,labels=marker1,srt=(rota[i]-180),xpd=NA,cex=$opts{cex},font=0.5)}
		#text(text_xx1,text_yy1,labels=name,srt=(rota[i]-180),xpd=NA,cex=$opts{cex},font=0.5)
		text(text_x2,text_y2,labels=marker2,srt=(rota[i]+90),xpd=NA,cex=$opts{cex},font=0.5);
	                    }else{
			    if(\"$opts{sample}\" == \"T\"){
		text(text_x1,text_y1,labels=paste(marker1,\"(\",name,\")\",sep=\"\"),srt=(rota[i]),xpd=NA,cex=$opts{cex},font=0.5)
		}else {text(text_x1,text_y1,labels=marker1,srt=(rota[i]),xpd=NA,cex=$opts{cex},font=0.5)}
		#	text(text_xx1,text_yy1,labels=name,srt=(rota[i]),xpd=NA,cex=$opts{cex},font=0.5)
		text(text_x2,text_y2,labels=marker2,srt=(rota[i]+90),xpd=NA,cex=$opts{cex},font=0.5);
			}
		}







x=cos(theta);y=sin(theta);
lines(x,y,type=\"l\",col=\"black\")#NA black
polygon(x,y,col=\"white\")
text(0,0,labels=test[1,2],xpd=T,cex=1,font=1);
dev.off()

png(file = \"flower.png\",width = $opts{w}*240,height=$opts{h}*200,res=300)
for (i in 1:aa){
x0=cos(space[i])*(a-r);
y0=sin(space[i])*(a-r);
alpha=space[i];
#text_x1=(2*a-1+1)*cos(space[i]);
#text_y1=(2*a-1+1)*sin(space[i]);
text_x1=(2*a-1+1+.9)*cos(space[i]);
text_y1=(2*a-1+1+.9)*sin(space[i]);

text_x2=(2*a-1.6)*cos(space[i]);
text_y2=(2*a-1.6)*sin(space[i]);
text_xx1=(2*a-1+0.5)*cos(space[i]);
	text_yy1=(2*a-1+0.5)*sin(space[i]);
	text_xx2=(2*a-1.5)*cos(space[i]);
	text_yy2=(2*a-1.5)*sin(space[i]);
	marker1=test[i,1];
	marker2=test[i,3];
	x3=x0+a*cos(alpha)*cos(theta)-b*sin(alpha)*sin(theta);
	y3=y0+a*sin(alpha)*cos(theta)+b*cos(alpha)*sin(theta);
	par(new=T)

        plot(x3,y3,type=\"l\",xlim=c(-8,8),ylim=c(-8,8),lwd=2,axes=F,xlab=\"\",ylab=\"\",col=c($sort_color)[i])
        if ($opts{type} == 2){
             polygon(x3,y3,col = alpha(c($sort_color)[i],0.3),border=NA)}

	     #name=paste(\"n=\",test[i,4])
#	if (rota[i] <= 180){
#		text(text_x1,text_y1,labels=marker1,srt=(rota[i]-90),xpd=NA,cex=$opts{cex},font=0.5)
#		#text(text_xx1,text_yy1,labels=name,srt=(rota[i]-90),xpd=NA,cex=$opts{cex},font=0.5)
#		text(text_x2,text_y2,labels=marker2,srt=(rota[i]-90),xpd=NA,cex=$opts{cex},font=0.5);
#		}else{
#		text(text_x1,text_y1,labels=marker1,srt=(rota[i]+90),xpd=NA,cex=$opts{cex},font=0.5)
#		#text(text_xx1,text_yy1,labels=name,srt=(rota[i]+90),xpd=NA,cex=$opts{cex},font=0.5)
#		text(text_x2,text_y2,labels=marker2,srt=(rota[i]+90),xpd=NA,cex=$opts{cex},font=0.5);
#			}
#		}


	if (rota[i] < 90) {
	        if(\"$opts{sample}\" == \"T\"){
		text(text_x1,text_y1,labels=paste(marker1,\"(\",name,\")\",sep=\"\"),srt=(rota[i]),xpd=NA,cex=$opts{cex},font=0.5)
		}else {text(text_x1,text_y1,labels=marker1,srt=(rota[i]),xpd=NA,cex=$opts{cex},font=0.5)}
		#text(text_xx1,text_yy1,labels=name,srt=(rota[i]),xpd=NA,cex=$opts{cex},font=0.5)
		text(text_x2,text_y2,labels=marker2,srt=(rota[i]-90),xpd=NA,cex=$opts{cex},font=0.5);
		}else if ((rota[i] >= 135) && (rota[i] <= 270)){
		if(\"$opts{sample}\" == \"T\"){
		text(text_x1,text_y1,labels=paste(marker1,\"(\",name,\")\",sep=\"\"),srt=(rota[i]+180),xpd=NA,cex=$opts{cex},font=0.5)
		}else {text(text_x1,text_y1,labels=marker1,srt=(rota[i]+180),xpd=NA,cex=$opts{cex},font=0.5)}
		#text(text_xx1,text_yy1,labels=name,srt=(rota[i]+90+90),xpd=NA,cex=$opts{cex},font=0.5)
		text(text_x2,text_y2,labels=marker2,srt=(rota[i]+90),xpd=NA,cex=$opts{cex},font=0.5);
	                    }else if ((rota[i] >= 90) && (rota[i] < 135)){
			    if(\"$opts{sample}\" == \"T\"){
		text(text_x1,text_y1,labels=paste(marker1,\"(\",name,\")\",sep=\"\"),srt=(rota[i]-180),xpd=NA,cex=$opts{cex},font=0.5)
		}else {text(text_x1,text_y1,labels=marker1,srt=(rota[i]-180),xpd=NA,cex=$opts{cex},font=0.5)}
		#text(text_xx1,text_yy1,labels=name,srt=(rota[i]-180),xpd=NA,cex=$opts{cex},font=0.5)
		text(text_x2,text_y2,labels=marker2,srt=(rota[i]+90),xpd=NA,cex=$opts{cex},font=0.5);
	                    }else{
			    if(\"$opts{sample}\" == \"T\"){
		text(text_x1,text_y1,labels=paste(marker1,\"(\",name,\")\",sep=\"\"),srt=(rota[i]),xpd=NA,cex=$opts{cex},font=0.5)
		}else {text(text_x1,text_y1,labels=marker1,srt=(rota[i]),xpd=NA,cex=$opts{cex},font=0.5)}
		#	text(text_xx1,text_yy1,labels=name,srt=(rota[i]),xpd=NA,cex=$opts{cex},font=0.5)
		text(text_x2,text_y2,labels=marker2,srt=(rota[i]+90),xpd=NA,cex=$opts{cex},font=0.5);
			}
		}






x=cos(theta);y=sin(theta);
lines(x,y,type=\"l\",col=NA)
text(0,0,labels=test[1,2],xpd=T,cex=1,font=1);
dev.off()
svg(file = \"flower.svg\",width = $opts{w},height=$opts{h})
for (i in 1:aa){
x0=cos(space[i])*(a-r);
y0=sin(space[i])*(a-r);
alpha=space[i];
#text_x1=(2*a-1+1)*cos(space[i]);
#text_y1=(2*a-1+1)*sin(space[i]);
text_x1=(2*a-1+1+.9)*cos(space[i]);
text_y1=(2*a-1+1+.9)*sin(space[i]);

text_x2=(2*a-1.6)*cos(space[i]);
text_y2=(2*a-1.6)*sin(space[i]);
text_xx1=(2*a-1+0.5)*cos(space[i]);
	text_yy1=(2*a-1+0.5)*sin(space[i]);
	text_xx2=(2*a-1.5)*cos(space[i]);
	text_yy2=(2*a-1.5)*sin(space[i]);
	marker1=test[i,1];
	marker2=test[i,3];
	x3=x0+a*cos(alpha)*cos(theta)-b*sin(alpha)*sin(theta);
	y3=y0+a*sin(alpha)*cos(theta)+b*cos(alpha)*sin(theta);
         par(new=T)
         plot(x3,y3,type=\"l\",xlim=c(-8,8),ylim=c(-8,8),lwd=2,axes=F,xlab=\"\",ylab=\"\",col=c($sort_color)[i])
         if ($opts{type} == 2){
                 polygon(x3,y3,col = alpha(c($sort_color)[i],0.3),border=NA)}

		 #name=paste(\"n=\",test[i,4])
#	if (rota[i] <= 180){
#		text(text_x1,text_y1,labels=marker1,srt=(rota[i]-90),xpd=NA,cex=$opts{cex},font=0.5)
#		#text(text_xx1,text_yy1,labels=name,srt=(rota[i]-90),xpd=NA,cex=$opts{cex},font=0.5)
#		text(text_x2,text_y2,labels=marker2,srt=(rota[i]-90),xpd=NA,cex=$opts{cex},font=0.5);
#		}else{
#		text(text_x1,text_y1,labels=marker1,srt=(rota[i]+90),xpd=NA,cex=$opts{cex},font=0.5)
#		#text(text_xx1,text_yy1,labels=name,srt=(rota[i]+90),xpd=NA,cex=$opts{cex},font=0.5)
#		text(text_x2,text_y2,labels=marker2,srt=(rota[i]+90),xpd=NA,cex=$opts{cex},font=0.5);
#			}
#		}

	if (rota[i] < 90) {
	        if(\"$opts{sample}\" == \"T\"){
		text(text_x1,text_y1,labels=paste(marker1,\"(\",name,\")\",sep=\"\"),srt=(rota[i]),xpd=NA,cex=$opts{cex},font=0.5)
		}else {text(text_x1,text_y1,labels=marker1,srt=(rota[i]),xpd=NA,cex=$opts{cex},font=0.5)}
		#text(text_xx1,text_yy1,labels=name,srt=(rota[i]),xpd=NA,cex=$opts{cex},font=0.5)
		text(text_x2,text_y2,labels=marker2,srt=(rota[i]-90),xpd=NA,cex=$opts{cex},font=0.5);
		}else if ((rota[i] >= 135) && (rota[i] <= 270)){
		if(\"$opts{sample}\" == \"T\"){
		text(text_x1,text_y1,labels=paste(marker1,\"(\",name,\")\",sep=\"\"),srt=(rota[i]+180),xpd=NA,cex=$opts{cex},font=0.5)
		}else {text(text_x1,text_y1,labels=marker1,srt=(rota[i]+180),xpd=NA,cex=$opts{cex},font=0.5)}
		#text(text_xx1,text_yy1,labels=name,srt=(rota[i]+90+90),xpd=NA,cex=$opts{cex},font=0.5)
		text(text_x2,text_y2,labels=marker2,srt=(rota[i]+90),xpd=NA,cex=$opts{cex},font=0.5);
	                    }else if ((rota[i] >= 90) && (rota[i] < 135)){
			    if(\"$opts{sample}\" == \"T\"){
		text(text_x1,text_y1,labels=paste(marker1,\"(\",name,\")\",sep=\"\"),srt=(rota[i]-180),xpd=NA,cex=$opts{cex},font=0.5)
		}else {text(text_x1,text_y1,labels=marker1,srt=(rota[i]-180),xpd=NA,cex=$opts{cex},font=0.5)}
		#text(text_xx1,text_yy1,labels=name,srt=(rota[i]-180),xpd=NA,cex=$opts{cex},font=0.5)
		text(text_x2,text_y2,labels=marker2,srt=(rota[i]+90),xpd=NA,cex=$opts{cex},font=0.5);
	                    }else{
			    if(\"$opts{sample}\" == \"T\"){
		text(text_x1,text_y1,labels=paste(marker1,\"(\",name,\")\",sep=\"\"),srt=(rota[i]),xpd=NA,cex=$opts{cex},font=0.5)
		}else {text(text_x1,text_y1,labels=marker1,srt=(rota[i]),xpd=NA,cex=$opts{cex},font=0.5)}
		#	text(text_xx1,text_yy1,labels=name,srt=(rota[i]),xpd=NA,cex=$opts{cex},font=0.5)
		text(text_x2,text_y2,labels=marker2,srt=(rota[i]+90),xpd=NA,cex=$opts{cex},font=0.5);
			}
		}



x=cos(theta);y=sin(theta);
lines(x,y,type=\"l\",col=NA)
text(0,0,labels=test[1,2],xpd=T,cex=1,font=1);
dev.off()
";
`R --no-save < cmd3.r`;
}








#my @mycol = ("#df89ff","#0000cd","#00c4ff","#ff8805","#ff5584","#00bd94","#d3b3b0","#4b0082","#c0c0c0","#ffd700","#8b0000","#00ffff","#ff0000","#0000cd","#006400","#ffff00","#008080","#d8bfd8","#40e0d0","#00ff7f","#6a5acd","#adff2f","#00ffff","#ff00ff","#8b4513","#6495ed","#ff6347","#800080","#dc143c","#000000","#7fff00","#d2691e","#ff7f50","#6495ed","#fff8dc","#dc143c","#00ffff")x200;

my @mycol = ("#df89ff","#0000cd","#00c4ff","#ff8805","#ff5584","#00bd94","#d3b3b0","#4b0082","#c0c0c0","#ffd700","#8b0000","#00ffff","#ff0000","#006400","#ffff00","#008080","#d8bfd8","#40e0d0","#00ff7f","#6a5acd","#adff2f","#ff00ff","#8b4513","#6495ed","#ff6347","#800080","#dc143c","#000000","#7fff00","#d2691e")x200;

if ($opts{m} eq "F"){
	open INPUT,$opts{i};

	chomp(my $h=<INPUT>);
	my @hh=split /\t/,$h;
	my @group_m=();
	my $col;
	#@group_m = @hh[1..$#hh];
	my %sort_col=(); 
	my @tmp;
	if ($opts{g} eq "F"){
		@group_m = @hh[1..$#hh];
		@tmp =sort @group_m;
           }
	else
       	{
	       	open TMM,$opts{g};
		chomp(@group_m =  <TMM>);
		@tmp = @group_m;
       	} 
	my $sort_group = join("\",\"",@group_m);
	$sort_group = "\"".$sort_group."\"";
	if ($opts{color} eq "F") { 
        $col = join("\",\"",@mycol[0..($#hh-1)]);
	$col = "\"".$col."\"";
	#my %sort_col=(); 
	
	my $ok=join("\t",@group_m); 
	foreach my $j (0..$#group_m) {
		$sort_col{$group_m[$j]} = $mycol[$j];
	}
	open COLOR,">color.txt";
	foreach my $i (0..$#group_m) {
		print COLOR "$tmp[$i]\t$sort_col{$tmp[$i]}\n";
	}
	close COLOR;
       }
       else {
	open COLORr,$opts{color};
        while(<COLORr>){	
		chomp;
		my @ab=split/\t/,$_;
		$sort_col{$ab[0]} = $ab[1] ;
	}
	my @tmp_col=();
	foreach my $j (@group_m) {
                 push @tmp_col,$sort_col{$j}
    }
    $col = join("\",\"",@tmp_col);
    $col = "\"".$col."\"";
 }


       my $ok = join("\t",@group_m) ;
 #my $sort_group = join("\",\"",@group_m);
 #	$sort_group= "\"".$sort_group."\"";

	my $sam_num = scalar @hh - 1;
	#print $sam_num,"----\n";	
        if ($sam_num <= 5) {
		&calculate_less_and_eq5sam_core($opts{i},$sort_group,$col);
	}
	else {
		&calculate_more5sam_core($opts{i},$ok);
		#print $ok,">>>\n";
                &flower("venn.sets.xls",$col);
	}
	
}


else{
	&calculate_group($opts{i},$opts{m},$opts{method});
	open INPUT,"tmp.$opts{method}.txt";
		chomp(my $h=<INPUT>);
	my @hh=split /\t/,$h;
	my @group_m=();
	my $col;
	#@group_m = @hh[1..$#hh];
	my %sort_col=(); 
	my @tmp;
	if ($opts{g} eq "F"){
		@group_m = @hh[1..$#hh];
		@tmp =sort @group_m;
           }
	else
       	{
	       	open TMM,$opts{g};
		chomp(@group_m =  <TMM>);
		@tmp = @group_m;
       	} 
	my $sort_group = join("\",\"",@group_m);
	$sort_group = "\"".$sort_group."\"";
	if ($opts{color} eq "F") { 
        $col = join("\",\"",@mycol[0..($#hh-1)]);
	$col = "\"".$col."\"";
	#my %sort_col=(); 
	
	my $ok=join("\t",@group_m); 
	foreach my $j (0..$#group_m) {
		$sort_col{$group_m[$j]} = $mycol[$j];
	}
	open COLOR,">color.txt";
	foreach my $i (0..$#group_m) {
		print COLOR "$tmp[$i]\t$sort_col{$tmp[$i]}\n";
	}
	close COLOR;
       }
       else {
	open COLORr,$opts{color};
        while(<COLORr>){	
		chomp;
		my @ab=split/\t/,$_;
		$sort_col{$ab[0]} = $ab[1] ;
	}
	my @tmp_col=();
	foreach my $j (@group_m) {
                 push @tmp_col,$sort_col{$j}
    }
    $col = join("\",\"",@tmp_col);
    $col = "\"".$col."\"";
 }


       my $ok = join("\t",@group_m) ;
 #my $sort_group = join("\",\"",@group_m);
 #	$sort_group= "\"".$sort_group."\"";

	my $sam_num = scalar @hh - 1; 
        if ($sam_num <= 5) {
		&calculate_less_and_eq5sam_core("tmp.$opts{method}.txt",$sort_group,$col);
	}
	else {
		&calculate_more5sam_core("tmp.$opts{method}.txt",$ok);
		print $ok,"\n";

                &flower("venn.sets.xls",$col);
	}
}


sub row2col2 {
	my ($input) = @_ ;
open A,$input;
open B,">venn.sets.element.xls";
my $num=0;  my %h;
my %count_col=();
while(<A>){
    chomp;my @a=split/\t/,$_;$h{$.}=[@a];
    #print [@a],"\$...............\n";

    $num=$#a;
    $count_col{$num} = 0;
}
close A;

my @ok=sort {$b<=>$a} keys %count_col;
#print "@ok\n";
$num=$ok[0];
open A,$input;
while(<A>){
	chomp;my @a=split/\t/,$_;
	my $tmp=scalar @a;
	my @ok;@ok=@a;
	if ($tmp < $num) {
		push @ok,("NA") x ($num-$tmp+1);
}
# print join("\t",@ok),"\n";
$h{$.}=[@ok]
 }
                  

 #=head;
foreach my $k1(0..$num){
    my @sam;
    foreach my $k2(sort {$a<=>$b} keys %h){
	    #print "!!!!$k2...................@{$h{$k2}}...OK............\n";
        push @sam,$h{$k2}[$k1];
    }
    #my $tmp=scalar @sam;
    #  if ($tmp < $num) {
    #    push @sam,("NA") x ($num-$tmp);}
    my $sam=join "\t",@sam;
    print B "$sam\n";
}
#=cut;
close A;close B;
}




#########Geat!!!!!!!!!!GuoqiLiu2024.06.25



 `rm cmd*.r`;
&row2col2("venn2dat.tsv");

 `rm venn2dat.tsv`;

sub row2col {
my ($input,$output) = @_; 
	open A,$input;
	open B,">$output"; 
my $num=0;  my %h;
while(<A>){
    chomp;my @a=split/\t/,$_;$h{$.}=[@a];
    $num=$#a;
}
close A;

foreach my $k1(0..$num){
    my @sam;
    foreach my $k2(sort {$a<=>$b} keys %h){
        push @sam,$h{$k2}[$k1];
    }
    my $sam=join "\t",@sam;
    print B "$sam\n";
 }
 close B;
}
