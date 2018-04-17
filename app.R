#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library("shiny")
library("DT")
library("clusterProfiler")
MirTarBase<-read.table("/home/lieng/Ressources/MirTarbase.MusMusculus.MirGenes.txt",
                       sep="\t",header=TRUE)
TarBase<-read.table("/home/lieng/Ressources/Tarbase.MusMusculus.MirGenes.txt",
                    sep="\t",header=TRUE)
TargetScan<-read.table("/home/lieng/Ressources/TargetScan.MusMusculus.MirGenes.txt",
                       sep="\t",header=TRUE)


OH.D<-read.table("/home/lieng/RESULTS/Forget_MiRNA_DeSEQ2_SequencageEtResequencage_DiscardSamples/TraitementEffect/One_Hour_Time.D_Localisation.TraitementEffect.txt",
                 header=TRUE,sep="\t",quote="\"")
OH.V<-read.table("/home/lieng/RESULTS/Forget_MiRNA_DeSEQ2_SequencageEtResequencage_DiscardSamples/TraitementEffect/One_Hour_Time.V_Localisation.TraitementEffect.txt",
                 header=TRUE,sep="\t",quote="\"")
TFH.D<-read.table("/home/lieng/RESULTS/Forget_MiRNA_DeSEQ2_SequencageEtResequencage_DiscardSamples/TraitementEffect/Twenty_Four_Hour_Time.D_Localisation.TraitementEffect.txt",
                  header=TRUE,sep="\t",quote="\"")
TFH.V<-read.table("/home/lieng/RESULTS/Forget_MiRNA_DeSEQ2_SequencageEtResequencage_DiscardSamples/TraitementEffect/Twenty_Four_Hour_Time.V_Localisation.TraitementEffect.txt",
                  header=TRUE,sep="\t",quote="\"")

D1.DS<-read.table("/home/lieng/RESULTS/D1_COCAINE_CURRENT/DorsalStriatum.BigResults.txt",
                  header=TRUE,sep="\t",quote="\"")
D2.DS<-read.table("/home/lieng/RESULTS/D2_COCAINE_CURRENT/DorsalStriatum.BigResults.txt",
                  header=TRUE,sep="\t",quote="\"")
D1.Nac<-read.table("/home/lieng/RESULTS/D1_COCAINE_CURRENT/NucleusAccumbens.BigResults.txt",
                   header=TRUE,sep="\t",quote="\"")
D2.Nac<-read.table("/home/lieng/RESULTS/D2_COCAINE_CURRENT/NucleusAccumbens.BigResults.txt",
                   header=TRUE,sep="\t",quote="\"")
Dir<-"/home/lieng/Ressources/EnrichRList/"
BP2Genes<-read.table(paste(Dir,"GO.BP.GeneList.txt",sep=""),header=TRUE,sep="\t")
CC2Genes<-read.table(paste(Dir,"GO.CC.GeneList.txt",sep=""),header=TRUE,sep="\t")
MF2Genes<-read.table(paste(Dir,"GO.MF.GeneList.txt",sep=""),header=TRUE,sep="\t")
Kegg2Genes<-read.table(paste(Dir,"Kegg.GeneList.txt",sep=""),header=TRUE,sep="\t")
P2Genes<-read.table(paste(Dir,"Panther.GeneList.txt",sep=""),header=TRUE,sep="\t")
R2Genes<-read.table(paste(Dir,"Reactome.DB.GeneList.txt",sep=""),header=TRUE,sep="\t")

GoTerms.BP<-read.table(paste(Dir,"GO.BP.Names.txt",sep=""),header=TRUE,sep="\t",quote="\"")
GoTerms.CC<-read.table(paste(Dir,"GO.CC.Names.txt",sep=""),header=TRUE,sep="\t",quote="\"")
GoTerms.MF<-read.table(paste(Dir,"GO.MF.Names.txt",sep=""),header=TRUE,sep="\t",quote="\"")
KeggTerms<-read.table(paste(Dir,"Kegg.Names.txt",sep=""),header=TRUE,sep="\t",quote="\"")
PTerms<-read.table(paste(Dir,"Panther.Names.txt",sep=""),header=TRUE,sep="\t",quote="\"")
ReactomeTerms<-read.table(paste(Dir,"Reactome.DB.Names.txt",sep=""),header=TRUE,sep="\t",quote="\"")





Mirs<-unique(c(as.vector(OH.D$Mir),as.vector(OH.V$Mir),
               as.vector(TFH.D$Mir),as.vector(TFH.V$Mir)))
Mirs<-sort(Mirs)
Mirs.List<-as.list(Mirs)
names(Mirs.List)<-Mirs

Mir.BaseMean.Max<-ceiling(max(log(c(OH.D$baseMean,
                                    OH.V$baseMean,
                                    TFH.D$baseMean,
                                    TFH.V$baseMean),10),na.rm = TRUE))
Mir.FC.Max<-ceiling(max(c(OH.D$log2FoldChange.Traitement_COC_vs_SAL.SAL_as_ref,
                                OH.V$log2FoldChange.Traitement_COC_vs_SAL.SAL_as_ref,
                                TFH.D$log2FoldChange.Traitement_COC_vs_SAL.SAL_as_ref,
                                TFH.V$log2FoldChange.Traitement_COC_vs_SAL.SAL_as_ref),na.rm = TRUE))
Mir.BMAxis.Value<- -1:Mir.BaseMean.Max
Mir.BMAxis.Labs<-parse(text=paste(10,"^",Mir.BMAxis.Value, sep=""))

Mir.FC.Lab<-substitute(paste(log[2]," ",frac(Num,Denom)),
                 list(Num="cocaine",
                      Denom="saline"))

RNA.BaseMean.Max<-ceiling(max(log(c(D1.DS$baseMean,
                                    D1.Nac$baseMean,
                                    D2.DS$baseMean,
                                    D2.Nac$baseMean),10),na.rm = TRUE))
RNA.FC.Max<-ceiling(max(c(D1.DS$log2FoldChange.Traitement_coc_vs_sal.sal_as_ref,
                          D2.DS$log2FoldChange.Traitement_coc_vs_sal.sal_as_ref,
                          D1.Nac$log2FoldChange.Traitement_coc_vs_sal.sal_as_ref,
                          D2.Nac$log2FoldChange.Traitement_coc_vs_sal.sal_as_ref),na.rm = TRUE))
RNA.BMAxis.Value<- -1:RNA.BaseMean.Max
RNA.BMAxis.Labs<-parse(text=paste(10,"^",RNA.BMAxis.Value, sep=""))

RNA.FC.Lab<-Mir.FC.Lab

HighlightAMir<-function(MirOfInterest="mmu-miR-1a-3p",Res=OH.D,Title="Cocaine effect"){
  Interest<-Res$Mir==MirOfInterest
  par(mar=c(5.1,7.1,4.1,2.1))
  plot(log(x=Res$baseMean,10),
       y=Res$log2FoldChange.Traitement_COC_vs_SAL.SAL_as_ref,
       col=as.vector(Res$col.Traitement_COC_vs_SAL.SAL_as_ref),
       xlab = "Mean expression (normalised counts)",
       ylab=Mir.FC.Lab,
       xlim=c(-1,Mir.BaseMean.Max),
       ylim=c(-Mir.FC.Max,Mir.FC.Max),
       xaxt='n',pch=20,main=Title)
  axis(side=1,at=Mir.BMAxis.Value,
       labels=Mir.BMAxis.Labs,cex.axis=1.25,cex.lab=1.2)
  
  XOI<-log(Res$baseMean[Interest],10)
  YOI<-Res$log2FoldChange.Traitement_COC_vs_SAL.SAL_as_ref[Interest]
  COI<-Res$col.Traitement_COC_vs_SAL.SAL_as_ref[Interest]
  points(x=XOI,y=YOI,col="blue",pch=42,cex=2)
  text(x=XOI,y=YOI,label=MirOfInterest,pos=4)
  grid(col="gray",lwd=2)
  par(mar=c(5.1,5.1,4.1,2.1))
}

HighlightSomeGenes<-function(RNAStudy=D1.DS,ListOfInterest=c("Penk","Drd1"),title="frezf"){
  OfInterest<-RNAStudy$Gene %in% ListOfInterest
  OfUninterest<-!OfInterest
  par(mar=c(5.1,7.1,4.1,2.1))
  plot(x=log(RNAStudy$baseMean[OfUninterest],10),
       y=RNAStudy$log2FoldChange.Traitement_coc_vs_sal.sal_as_ref[OfUninterest],
       col=RNAStudy$col.Traitement_coc_vs_sal.sal_as_ref[OfUninterest],
       xlim=c(-1,RNA.BaseMean.Max),ylim=c(-RNA.FC.Max,RNA.FC.Max),
       main=title,xaxt="n",
       xlab="Mean expression (normalised counts)",
       ylab=RNA.FC.Lab,pch=".")
  axis(side=1,at=RNA.BMAxis.Value,
       labels=RNA.BMAxis.Labs,cex.axis=1.25,cex.lab=1.2)
  points(x=log(RNAStudy$baseMean[OfInterest],10),
         y=RNAStudy$log2FoldChange.Traitement_coc_vs_sal.sal_as_ref[OfInterest],
         col=RNAStudy$col.Traitement_coc_vs_sal.sal_as_ref[OfInterest],pch=42,cex=2)
  grid(col="gray",lwd=2)
  
}


# Define UI for application that draws a histogram
ui <- fluidPage(
   # Application title
   titlePanel("Look at cocaine effect"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
        selectInput("Mir",
                    "Select micro rna of interest", 
                    choices = Mirs.List,
                    selected = "mmu-miR-1a-3p")
      ,checkboxGroupInput("TargetDatabase", 
                           h3("Mir target database"), 
                           choices = list("MirTarBase" = 1, 
                                          "TarBase" = 2, 
                                          "TargetScan" = 3),
                           selected = c(1,2,3)),
      selectInput("Annotation",
                  "Select annotation database", 
                  choices = list("GO Molecular Function" = "MF",
                                 "GO Biological Process" = "BP",
                                 "GO Cellular Compartment" = "CC",
                                 "KEGG" = "KEGG",
                                 "Panther" = "PNTHR",
                                 "Reactome" = "RCTM"),
                  selected = "PNTHR")
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("MirAndRNA"),
         DT::dataTableOutput("MirTargetsAnnot"),
         plotOutput("Targets")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   output$MirAndRNA <- renderPlot({
     par(mfrow=c(2,2))
     HighlightAMir(Mir=input$Mir,
                   Res=OH.D,Title="1h - Dorsal Striatum")
     HighlightAMir(Mir=input$Mir,
                   Res=OH.V,Title="1h - Ventral Striatum")
     HighlightAMir(Mir=input$Mir,
                   Res=TFH.D,Title="24h - Dorsal Striatum")
     HighlightAMir(Mir=input$Mir,
                   Res=TFH.V,Title="24h - Ventral Striatum")
     
   })
   output$MirTargetsAnnot <- DT::renderDataTable(DT::datatable({
     #cat(input$TargetDatabase,"\n")
     TargetedGenes<-c()
     if(1 %in% input$TargetDatabase){
       TargetedGenes<-c(TargetedGenes,as.vector(MirTarBase[MirTarBase$mirna==input$Mir,"GeneSymbol"]))
     }
     if(2 %in% input$TargetDatabase){
       TargetedGenes<-c(TargetedGenes,as.vector(TarBase[TarBase$mirna==input$Mir,"GeneSymbol"]))
     }
     if(3 %in% input$TargetDatabase){
       TargetedGenes<-c(TargetedGenes,as.vector(TargetScan[TargetScan$mirna==input$Mir,"GeneSymbol"]))
     }
     TargetedGenes<-unique(TargetedGenes)
     if(input$Annotation=="RCTM"){
       data<-enricher(TargetedGenes,TERM2GENE=R2Genes,TERM2NAME=ReactomeTerms,
                      minGSSize=1,pvalueCutoff=0.1,qvalueCutoff=0.1,pAdjustMethod="fdr")
     }else if(input$Annotation=="PNTHR"){
       data<-enricher(TargetedGenes,TERM2GENE=P2Genes,TERM2NAME=PTerms,
                      minGSSize=1,pvalueCutoff=0.1,qvalueCutoff=0.1,pAdjustMethod="fdr")
     }else if(input$Annotation=="KEGG"){
       data<-enricher(GeneList,TERM2GENE=Kegg2Genes,TERM2NAME=KeggTerms,
                      minGSSize=1,pvalueCutoff=0.1,qvalueCutoff=0.1,pAdjustMethod="fdr")
     }else if(input$Annotation=="CC"){
       data<-enricher(TargetedGenes,TERM2GENE=CC2Genes,TERM2NAME=GoTerms.CC,
                      minGSSize=1,pvalueCutoff=0.1,qvalueCutoff=0.1,pAdjustMethod="fdr")
     }else if(input$Annotation=="BP"){
       data<-enricher(GeneList,TERM2GENE=CC2Genes,TERM2NAME=GoTerms.BP,
                      minGSSize=1,pvalueCutoff=0.1,qvalueCutoff=0.1,pAdjustMethod="fdr")
     }else if(input$Annotation=="MF"){
       data<-enricher(TargetedGenes,TERM2GENE=CC2Genes,TERM2NAME=GoTerms.MF,
                      minGSSize=1,pvalueCutoff=0.1,qvalueCutoff=0.1,pAdjustMethod="fdr")
     }
     data<-data.frame(data)
   }))
   output$Targets <- renderPlot({
     TargetedGenes<-c()
     if(1 %in% input$TargetDatabase){
       TargetedGenes<-c(TargetedGenes,as.vector(MirTarBase[MirTarBase$mirna==input$Mir,"GeneSymbol"]))
     }
     if(2 %in% input$TargetDatabase){
       TargetedGenes<-c(TargetedGenes,as.vector(TarBase[TarBase$mirna==input$Mir,"GeneSymbol"]))
     }
     if(3 %in% input$TargetDatabase){
       TargetedGenes<-c(TargetedGenes,as.vector(TargetScan[TargetScan$mirna==input$Mir,"GeneSymbol"]))
     }
     TargetedGenes<-unique(TargetedGenes)
     par(mfrow=c(2,2))
     HighlightSomeGenes(RNAStudy=D1.DS,
                        title="Cocaine effect on D1 Dorsal Striatum")
     HighlightSomeGenes(RNAStudy=D2.DS,
                        ListOfInterest=TargetedGenes,
                        title="Cocaine effect on D2 Dorsal Striatum")
     HighlightSomeGenes(RNAStudy=D1.Nac,
                        ListOfInterest=TargetedGenes,
                        title="Cocaine effect on D1 Ventral Striatum")
     HighlightSomeGenes(RNAStudy=D2.Nac,
                        ListOfInterest=TargetedGenes,
                        title="Cocaine effect on D2 Ventral Striatum")
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

