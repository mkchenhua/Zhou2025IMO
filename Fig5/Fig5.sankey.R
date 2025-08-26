#https://zhuanlan.zhihu.com/p/485914410
#https://www.jianshu.com/p/d103b4b143f3
#http://cn.voidcc.com/question/p-spgsnnxm-uo.html
#https://mp.weixin.qq.com/s/P2CkrmKkFaqQi9K1SqCZ3Q
#https://zhuanlan.zhihu.com/p/138318632
#install.packages("devtools")

#devtools::install_github("fbreitwieser/sankeyD3")

library(sankeyD3)
#links <-  read.table("t1.txt",header = T,sep = "\t",check.names=F,comment.char="",quote="")
links <-  read.table("zok.txt",header = T,sep = "\t",check.names=F,comment.char="",quote="")
colnames(links)[1:2]<-c("source","target")
#nodes <- data.frame(name=c(as.character(links$source), as.character(links$target)) %>% unique())
nodes <- data.frame(name=unique(c(as.character(links$source), as.character(links$target))))
#然后基于nodes数据框构建links中节点的唯一标识符ID，而非根据节点的name
links$IDsource <- match(links$source, nodes$name)-1
links$IDtarget <- match(links$target, nodes$name)-1
links$weight <-abs(links$value)
 links$group <- "no"
links$group[links$value>=0]<-"yes"
sankeyNetwork(Links = links, Nodes = nodes, Source = "IDsource", Target = "IDtarget",
Value = "weight", NodeID = "name",nodeWidth =10,units = 'TWh',
height=300,width=300,colourScale=JS("d3.scaleOrdinal(d3.schemeCategory10);"),
numberFormat=".2f",fontSize = 8)
#NodeGroup="group",LinkGroup = "group",

#links$group <- "linkgrp"

######统一颜色###
#colourScale <- 
#  'd3.scaleOrdinal()
#     .domain(["linkgrp"])
#     .range(["blue","red"].concat(d3.schemeCategory20))'
##################################
colourScale <- 'd3.scaleOrdinal() .domain(["yes","no"]) .range(["#ff7f50", "#6495ed"])'

# Make a look up table of events and colors


#nodes$color<-sample(c("red","orange","blue","green"),nrow(nodes),replace=T) 
nodes$color <- rainbow(length(nodes$name))


#nodeStrokeWidth node边框宽度
#linkType character One of 'bezier', 'l-bezier', 'trapezoid', 'path1' and 'path2'.
#align	character Alignment of the nodes. One of 'right', 'left', 'justify', 'center', 'none'. 
#If 'none', then the labels of the nodes are always to the right of the node.
write.table(links,file = "sankey.xls",sep="\t",quote=FALSE,row.names=FALSE)
p<- sankeyNetwork(Links = links, Nodes = nodes, Source = "IDsource", Target = "IDtarget",
                  Value = "weight", NodeID = "name",nodeWidth =50,units = 'TWh',
                  height=700,width=1200,nodeStrokeWidth = 2,linkType="bezier",nodePadding = 3,
                  align="justify",showNodeValues=F,
                  colourScale=colourScale,NodeColor="color",
                  #colourScale=JS("d3.scaleOrdinal(d3.schemeCategory10);"),
                  numberFormat=".2f",fontSize = 12,
                  LinkGroup = "group")
#height=500,width=800,nodeWidth =60
#nodePadding = 5,柱子的间距

saveNetwork(p,"sankey3.html")
library(webshot)
 if(!is_phantomjs_installed()){

  install_phantomjs()

}


webshot("sankey3.html" , "sankey3.pdf")






















#nodeStrokeWidth node边框宽度
#linkType character One of 'bezier', 'l-bezier', 'trapezoid', 'path1' and 'path2'.
#align	character Alignment of the nodes. One of 'right', 'left', 'justify', 'center', 'none'. 
#If 'none', then the labels of the nodes are always to the right of the node.
write.table(links,file = "sankey.xls",sep="\t",quote=FALSE,row.names=FALSE)
s.g <- readLines("zlst.txt")
library(dplyr)
#df <- data.frame(x=1:5,y=c("c","c","b","d","a"))
nodes <- arrange(nodes,match (nodes$name ,s.g ) )
tmp <- nodes$name 
#nodes$name <- factor(nodes$name ,tmp)
#links$IDsource <- match(links$source, nodes$name)-1
#links$IDtarget <- match(links$target, nodes$name)-1

as.numeric(rownames(nodes))-1

links$IDsource <- match(links$source, nodes$name)
links$IDtarget <- match(links$target, nodes$name)
nodes$index <- as.numeric(rownames(nodes))-1
nodes$source <- nodes$name
links<-links[,c("source","target","value","weight","group")]
links <-merge(links,nodes,by="source")
colnames(links)[8] <- "IDsource"
nodes$target <- nodes$name
links <-merge(links,nodes,by="target")
colnames(links)[11] <- "IDtarget"
 colnames(links)[2] <- "source"
nodes <- nodes[,1:2]
links <- links[,c("source","target","value","IDsource","IDtarget","weight","group")]
write.table(links,file = "sankey.xls",sep="\t",quote=FALSE,row.names=FALSE)
#links <- read.table("t.txt",header = T,sep = "\t",check.names=F,comment.char="",quote="")
# 'right', 'left', 'justify', 'center', 'none'. If 'none', then the labels of the nodes are always to the right of the node.
nodes$ID <- nodes$name
write.table(nodes,file = "node.xls",sep="\t",quote=FALSE,row.names=FALSE)
p<- sankeyNetwork(Links = links, Nodes = nodes, Source = "IDsource", Target = "IDtarget",
                  Value = "weight", NodeID = "name",nodeWidth =100,units = 'TWh',
                  height=600,width=800,nodeStrokeWidth = 3,linkType="bezier",
                  showNodeValues=F,
                  colourScale=colourScale,
                   #NodeColor="color",
                  #colourScale=JS("d3.scaleOrdinal(d3.schemeCategory10);"),
                  #numberFormat=".2f",
                 fontSize = 12,
                  LinkGroup = "group",iterations =-100)

p<- sankeyNetwork(Links = links, Nodes = nodes, Source = "IDsource", Target = "IDtarget",
                  Value = "weight", NodeID = "name",
                 fontSize = 12, colourScale=colourScale,
align="left",
                  LinkGroup = "group")
library(networkD3)
library(htmlwidgets)
customJS <- 'function() { console.log(this.sankey.nodes().map(d => [d.name, d.x, d.y])); }'
onRender(p, customJS)
##################################################

htmlwidgets::onRender(p, '
      function(el) { 
        var cols_x = this.sankey.nodes().map(d => d.x).filter((v, i, a) => a.indexOf(v) === i).sort();
        var labels = ["Location", "Office", "Strategy", "Tactic",];
        cols_x.forEach((d, i) => {
          d3.select(el).select("svg")
            .append("text")
            .attr("x", d)
            .attr("y", 12)
            .text(labels[i]);

        })
      }
    ')


sankeyNetwork(Links = links, Nodes = nodes, Source = "source",
              Target = "target", Value = "weight", NodeID = "name",
              units = "TWh", fontSize = 8, nodeWidth = 30)

################################
    plot <- sankeyNetwork(Links = links, Nodes = nodes,
                  Source = "IDsource", Target = "IDtarget", Value = 'value', NodeID = 'name',
                  fontSize = 14)
    
    #Apply the manual var labels - solution from the linked stackoverflow answer
    htmlwidgets::onRender(plot, '
      function(el) { 
        var cols_x = this.sankey.nodes().map(d => d.x).filter((v, i, a) => a.indexOf(v) === i);
        var labels = ["Location", "Office", "Strategy", "Tactic", "Target", "Outcome", "Output"];
        cols_x.forEach((d, i) => {
          d3.select(el).select("svg")
            .append("text")
            .attr("x", d)
            .attr("y", 12)
            .text(labels[i]);
        })
      }
    ')
################################


 sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource",Target = "IDtarget",Value = "value",
                   NodeID = "name",

                   LinkGroup = 'group',
                   colourScale= ColourScal, 
                   # nodeWidth=40,
                   # fontSize=13
                   # nodePadding=20
)






saveNetwork(p,"sankey3.html")
library(webshot)
 if(!is_phantomjs_installed()){

  install_phantomjs()

}


webshot("sankey3.html" , "sankey3.pdf")
###################################################
library(networkD3)
library(htmlwidgets)
p <- sankeyNetwork(Links = links, Nodes =nodes, Source = 'IDsource',
              Target = 'IDtarget', Value = 'weight', NodeID = 'name',
              units = 'TWh', fontSize = 12, nodeWidth = 30)

customJS <- 'function() { console.log(this.sankey.nodes().map(d => [d.name, d.x, d.y])); }'
onRender(p, customJS)


