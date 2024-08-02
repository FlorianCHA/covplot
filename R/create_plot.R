data.snp.process <- function(path.data, sheet.name, min.depth)
{
  data <- readxl::read_excel(path.data, sheet = sheet.name)
  samples <- unique(data$Sample)
  coverage <- data[order(data$Region, decreasing=TRUE),]
  coverage <- coverage[!duplicated(coverage$Region),]
  position_coverage <- coverage$Region
  coverage <- coverage[samples]
  rownames(coverage) <- position_coverage

  # Retrieve only SNP (Reference != Allele)
  data_snp <- data %>%
    filter(Reference != Allele)

  # Retrieve all sample coverage (assuming columns 58 onwards hold sample data)
  list_sample <- names(data)[59:ncol(data)-3]

  # List of headers to keep
  list_header <- c("Sample", "Region", "Allele", "Frequency",'Amino acid change in longest transcript')

  # Select desired columns for analysis
  data_snp <- data_snp %>%
    select(list_header)

  # Filter data with non-missing amino acid changes
  data_snp <- data_snp %>%
    mutate(`Amino acid change in longest transcript` = sapply(strsplit(`Amino acid change in longest transcript`, "p.",fixed=TRUE), function(x) x[2]))  # Extract amino acid change
  #    filter(!is.na(`Amino acid change in longest transcript`)) %>%

  # Identifiez les positions ayant moins de 100 de profondeur.
  val_under_100 <- which(coverage < min.depth, arr.ind=T)
  if (length(val_under_100) >0)
  {
    for (indice in 1:(length(val_under_100)/2))
    {
      indice.row <- val_under_100[indice,1]
      indice.col <- val_under_100[indice,2]
      region <- position_coverage[indice.row]
      sample <- colnames(coverage[indice.col])
      allele <- NA
      frequency <- NA
      aa <- NA
      data_snp <- rbind(data_snp,c(sample, region, allele, frequency, aa))
    }
  }

  # Keep only SNV type
  data_snp <- data_snp %>%
    mutate(Region =  sapply(str_replace(Region,fixed('^'),'..'), function(x) x)) %>%  # Split and convert Region to integer
    mutate(Region =  sapply(strsplit(Region,'..',fixed=TRUE), function(x) as.integer(x[1]))) %>%  # Split and convert Region to integer
    mutate(Sample = sapply(strsplit(Sample,'_'), function(x) x[2])) # Extract sample name
  data_snp$Region <- as.integer(data_snp$Region)
  data_snp$Frequency <- as.integer(data_snp$Frequency)
  data_snp$`Amino acid change in longest transcript`[is.na(data_snp$`Amino acid change in longest transcript`)] <- ""
  return(data_snp)
}

make.orf <- function(data, p)
{
  xmin <- c()
  xmax <- c()
  ymin <- c()
  ymax <- c()
  labs <- c()
  for (i in 1:length(unique(data$Sample)))
  {
    i = i-1
    # Add ORF1ab
    xmin <- c(xmin, 266-1)
    xmax <- c(xmax, 21555)
    ymin <- c(ymin, 0.65+i)
    ymax <- c(ymax, 1.35+i)
    labs <- c(labs,'ORF1ab')
    # Add S
    xmin <- c(xmin, 21563-1)
    xmax <- c(xmax, 25384)
    ymin <- c(ymin, 0.65+i)
    ymax <- c(ymax, 1.35+i)
    labs <- c(labs,'S')
    # Add ORF3a
    xmin <- c(xmin, 25393-1)
    xmax <- c(xmax, 26220)
    ymin <- c(ymin, 0.65+i)
    ymax <- c(ymax, 1.35+i)
    labs <- c(labs,'ORF3a')
    # Add E
    xmin <- c(xmin, 26245-1)
    xmax <- c(xmax, 26472)
    ymin <- c(ymin, 0.65+i)
    ymax <- c(ymax, 1.35+i)
    labs <- c(labs,'E')
    # Add M
    xmin <- c(xmin, 26523-1)
    xmax <- c(xmax, 27191)
    ymin <- c(ymin, 0.65+i)
    ymax <- c(ymax, 1.35+i)
    labs <- c(labs,'M')
    # Add ORF6
    xmin <- c(xmin, 27202-1)
    xmax <- c(xmax, 27387)
    ymin <- c(ymin, 0.65+i)
    ymax <- c(ymax, 1.35+i)
    labs <- c(labs,'ORF6')
    # Add ORF7a
    xmin <- c(xmin, 27394-1)
    xmax <- c(xmax, 27759)
    ymin <- c(ymin, 0.65+i)
    ymax <- c(ymax, 1.35+i)
    labs <- c(labs,'ORF7a')
    # Add ORF8
    xmin <- c(xmin, 27894-1)
    xmax <- c(xmax, 28259)
    ymin <- c(ymin, 0.65+i)
    ymax <- c(ymax, 1.35+i)
    labs <- c(labs,'ORF8')
    # Add N
    xmin <- c(xmin, 28274-1)
    xmax <- c(xmax, 29533)
    ymin <- c(ymin, 0.65+i)
    ymax <- c(ymax, 1.35+i)
    labs <- c(labs,'N')
    # Add ORF10
    xmin <- c(xmin, 29558-1)
    xmax <- c(xmax, 29674)
    ymin <- c(ymin, 0.65+i)
    ymax <- c(ymax, 1.35+i)
    labs <- c(labs,'ORF10')
  }

  rect_df <- data.frame(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                        labs = labs)


  p <- p + geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = labs),
                     alpha = 0.5, data = rect_df, inherit.aes = FALSE)

  return(p)
}


#' Render Covid plot
#'
#' Create a Covid plot for show SNP evolution in genome
#' @param path.data The path of the excel file to use
#' @param sheet.name The sheet nome of the excel to process
#' @param min.frequency The min frequency (SNPs with a frequency below the threshold will not be displayed)
#' @param min.depth The minimum depth, positions with less depth than the threshold will be displayed in white on the plot.
#' @param show This option allows you to select the items to be displayed, by default the function show all genome. You can choose between 'ORF1ab', 'S', 'ORF3a', 'E', 'M', 'ORF6', 'ORF7a', 'ORF8', 'N', 'ORF10'
#' @return A beautiful plot
#' @import ggplot2
#' @import dplyr
#' @import gtools
#' @import readxl
#' @examples
#' render.plot(path.data='/path/to/excel', sheet.name='final', min.frequency=10, min.depth=100)
#' @export
render.plot <- function(path.data, sheet.name, path.output.fig, min.frequency = 0, min.depth = 100, show='all')
{
  data <- data.snp.process(path.data, sheet.name, min.depth)
  file.name <- strsplit(basename(path.data),".",fixed=TRUE)[[1]][1]
  table.coverage <- as.data.frame(table(data$Region))
  colnames(table.coverage) <- c('Region','Nb SNP or Missing')
  write.table(table.coverage,paste(file.name,'.csv',sep=''),row.names=FALSE)
  data <- data[data$Frequency>= min.frequency | is.na(data$Frequency),]
  p <- ggplot2::ggplot(data,ggplot2::aes(x=Region, y=Sample, color=Frequency)) + ggplot2::theme_classic() + ggplot2::geom_point()+ ggplot2::scale_color_gradient(na.value = "white", low="grey82", high="black")

  #Make colored rect for orf in genome
  p <- make.orf(data,p)

  p <- p + ggplot2::geom_point() + ggplot2::theme(  axis.ticks = ggplot2::element_blank(),plot.title = ggplot2::element_text(hjust = 0.5,face='bold'),
              axis.title= ggplot2::element_blank(),
              axis.text.x = ggplot2::element_text(angle=90,size=6,face="bold",vjust = 0.5),
              axis.line.y = ggplot2::element_blank()) + ggplot2::xlim(min(data$Region)-1, 29674+1) +
    ggplot2::scale_x_continuous(breaks=data$Region,label=data$`Amino acid change in longest transcript`)+
    ggplot2::scale_y_discrete(breaks=unique(data$Sample)[gtools::mixedorder(unique(data$Sample))],
          limits=unique(data$Sample)[gtools::mixedorder(unique(data$Sample),decreasing = TRUE)])+
    ggplot2::scale_fill_manual(breaks=c('ORF1ab','S','ORF3a','E','M','ORF6','ORF7a','ORF8','N','ORF10'),
                      values=c('seagreen4','deeppink1','turquoise4','palevioletred2','lightsalmon2','cadetblue3','darkorchid2','violetred2','yellowgreen','forestgreen'),
                      limits=c('ORF1ab','S','ORF3a','E','M','ORF6','ORF7a','ORF8','N','ORF10')) +
    ggplot2::labs(fill = "")+
    ggplot2::ggtitle(paste('SNP analysis for the',file.name,'file',sep=' '))

  if (show == 'ORF1ab')
  { p <- p + ggplot2::scale_x_continuous(breaks=data$Region,label=data$`Amino acid change in longest transcript`,limits = c(266-1, 21555))+
    ggplot2::ggtitle(paste('SNP analysis for the',file.name,'file, focus on ORF1ab',sep=' '))}
  if (show == 'S')
  { p <- p + ggplot2::scale_x_continuous(breaks=data$Region,label=data$`Amino acid change in longest transcript`,limits = c(21563-1, 25384))+
    ggplot2::ggtitle(paste('SNP analysis for the',file.name,'file, focus on S',sep=' '))}
  if (show == 'ORF3a')
  { p <- p + ggplot2::scale_x_continuous(breaks=data$Region,label=data$`Amino acid change in longest transcript`,limits = c(25393-1, 26220))+
    ggplot2::ggtitle(paste('SNP analysis for the',file.name,'file, focus on ORF3',sep=' '))}
  if (show == 'E')
  { p <- p + ggplot2::scale_x_continuous(breaks=data$Region,label=data$`Amino acid change in longest transcript`,limits = c(26245-1, 26472))+
    ggplot2::ggtitle(paste('SNP analysis for the',file.name,'file, focus on E',sep=' '))}
  if (show == 'M')
  { p <- p + ggplot2::scale_x_continuous(breaks=data$Region,label=data$`Amino acid change in longest transcript`,limits = c(26523-1, 27191))+
    ggplot2::ggtitle(paste('SNP analysis for the',file.name,'file, focus on M',sep=' '))}
  if (show == 'ORF6')
  { p <- p + ggplot2::scale_x_continuous(breaks=data$Region,label=data$`Amino acid change in longest transcript`,limits = c(27202-1, 27387))+
    ggplot2::ggtitle(paste('SNP analysis for the',file.name,'file, focus on ORF6',sep=' '))}
  if (show == 'ORF7a')
  { p <- p + ggplot2::scale_x_continuous(breaks=data$Region,label=data$`Amino acid change in longest transcript`,limits = c(27394-1, 27759))+
    ggplot2::ggtitle(paste('SNP analysis for the',file.name,'file, focus on ORF7a',sep=' '))}
  if (show == 'ORF8')
  { p <- p + ggplot2::scale_x_continuous(breaks=data$Region,label=data$`Amino acid change in longest transcript`,limits = c(27894-1, 28259))+
    ggplot2::ggtitle(paste('SNP analysis for the',file.name,'file, focus on ORF8',sep=' '))}
  if (show == 'N')
  { p <- p + ggplot2::scale_x_continuous(breaks=data$Region,label=data$`Amino acid change in longest transcript`,limits = c(28274-1, 29533))+
    ggplot2::ggtitle(paste('SNP analysis for the',file.name,'file, focus on N',sep=' '))}
  if (show == 'ORF10')
  { p <- p + ggplot2::scale_x_continuous(breaks=data$Region,label=data$`Amino acid change in longest transcript`,limits = c(29558-1, 29674))+
    ggplot2::ggtitle(paste('SNP analysis for the',file.name,'file, focus on ORF10',sep=' '))}

  # Save at pdf format
  pdf(paste(file.name,'.pdf',sep=''),width=11,height=6)
  print(p)
  dev.off()
  return(p)
}


#' Render Heatmap Covid plot
#'
#' Create a Covid heatmap for show SNP evolution in genome
#' @param path.data The path of the excel file to use
#' @param sheet.name The sheet nome of the excel to process
#' @param min.depth The minimum depth, positions with less depth than the threshold will be displayed in white on the plot.
#' @param show This option allows you to select the items to be displayed (acide nucleic (ac) or amino acid (aa)), by default the function show ac information. You can choose between 'ac' and 'aa'
#' @return A beautiful plot
#' @import ggplot2
#' @import dplyr
#' @import reshape2
#' @import stringr
#' @import gtools
#' @import readxl
#' @examples
#' render.heatmap(path.data='/path/to/excel', sheet.name='final',  min.depth=100)
#' @export
render.heatmap <- function(path.data, sheet.name, min.depth = 100, show='ac')
{
  data.snp <- data.snp.process(path.data, sheet.name,min.depth=100)
  file.name <- strsplit(basename(path.data),".",fixed=TRUE)[[1]][1]

  region <- as.vector(unique(data.snp[data.snp$`Amino acid change in longest transcript` !='','Region']))$Region
  region.x.legend <- data.snp[data.snp$`Amino acid change in longest transcript` !='',c('Region','Amino acid change in longest transcript')]
  region.x.legend <- region.x.legend[order(region.x.legend$Region),]
  region.x.legend <- region.x.legend[!duplicated(region.x.legend$Region),]
  data.aa <- data.snp[data.snp$Region %in% region,]
  matrice.aa <- acast(data.aa, Sample~Region, value.var = 'Frequency', fun.aggregate=sum)
  data.aa.plot <- setNames(melt(matrice.aa), c('Sample', 'Region', 'Frequency'))

  plot.aa <- ggplot(data = data.aa.plot, aes(x=as.character(Region), y=Sample, fill=Frequency),color='black') +
    geom_tile(color='black',linewidth=0.5) +
    scale_y_discrete(breaks=unique(data.snp$Sample)[mixedorder(unique(data.snp$Sample),decreasing=TRUE)],
                     limits=unique(data.snp$Sample)[mixedorder(unique(data.snp$Sample),decreasing = TRUE)]) +
    scale_x_discrete(breaks = region.x.legend$Region, label =region.x.legend$`Amino acid change in longest transcript`) +theme_classic()+
    theme(axis.ticks = element_blank(),
          axis.title= element_blank(),
          axis.text.x = element_text(angle=90,size=6,face="bold",vjust = 0.5))+
    scale_fill_gradient(low = "white", high = "steelblue", na.value = 'grey80')

  matrice.all <- acast(data.snp, Sample~Region, value.var = 'Frequency', fun.aggregate=sum)
  data.all.plot <- setNames(melt(matrice.all), c('Sample', 'Region', 'Frequency'))
  data.all.plot$Region <- as.character(data.all.plot$Region)
  plot.ac <- ggplot(data = data.all.plot, aes(x=as.character(Region), y=Sample, fill=Frequency),color='black') +
    geom_tile(color='black',linewidth=0.5) +
    scale_y_discrete(breaks=unique(data.snp$Sample)[mixedorder(unique(data.snp$Sample),decreasing=TRUE)],
                     limits=unique(data.snp$Sample)[mixedorder(unique(data.snp$Sample),decreasing = TRUE)]) +
    scale_x_discrete(breaks = unique(data.all.plot$Region)[gtools::mixedorder(unique(data.all.plot$Region))],
                     limits = unique(data.all.plot$Region)[gtools::mixedorder(unique(data.all.plot$Region))]) +theme_classic()+
    theme(axis.ticks = element_blank(),
          axis.title= element_blank(),
          axis.text.x = element_text(angle=90,size=6,face="bold",vjust = 0.5))+
    scale_fill_gradient(low = "white", high = "steelblue", na.value = 'grey80')
  if (show == 'ac')
  {
    pdf(paste(file.name,'_ac.pdf',sep=''),width=11,height=6)
    print(plot.ac)
    dev.off()
    return(plot.ac)
  }
  if (show == 'aa')
  {
    pdf(paste(file.name,'_aa.pdf',sep=''),width=11,height=6)
    print(plot.aa)
    dev.off()
    return(plot.aa)
  }
}
