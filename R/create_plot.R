data.snp.process <- function(path.data, sheet.name, min.depth)
{
  data <- readxl::read_excel(path.data, sheet = sheet.name)
  samples <- unique(data$Sample)
  coverage <- data[samples]


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
    mutate(`Amino acid change in longest transcript` = sapply(strsplit(`Amino acid change in longest transcript`, "p."), function(x) x[2]))  # Extract amino acid change
  #    filter(!is.na(`Amino acid change in longest transcript`)) %>%

  # Identifiez les positions des valeurs supérieures à 100 dans chaque colonne
  val_under_100 <- which(coverage < min.depth, arr.ind=T)
  if (length(val_under_100) >0)
  {
    for (indice in 1:(length(val_under_100)/2))
    {
      indice.row <- val_under_100[indice,1]
      indice.col <- val_under_100[indice,2]
      region <- data$Region[indice.row]
      sample <- colnames(coverage[indice.col])
      allele <- NA
      frequency <- NA
      aa <- NA
      data_snp <- rbind(data_snp,c(sample, region, allele, frequency, aa))
    }
  }

  # Keep only SNV type
  data_snp <- data_snp %>%
    mutate(Region =  sapply(strsplit(Region,'..',fixed=TRUE), function(x) as.integer(x[1]))) %>%  # Split and convert Region to integer
    mutate(Sample = sapply(strsplit(Sample,'_'), function(x) x[2])) # Extract sample name
  data_snp$Region <- as.integer(data_snp$Region)
  data_snp$Frequency <- as.integer(data_snp$Frequency)
  data_snp$`Amino acid change in longest transcript`[is.na(data_snp$`Amino acid change in longest transcript`)] <- ""
  return(data_snp)
}



#' Render Covid plot
#'
#' Create a Covid plot for show SNP evolution in genome
#' @param path.data The path of the excel file to use
#' @param sheet.name The sheet nome of the excel to process
#' @param min.frequency The min frequency (SNPs with a frequency below the threshold will not be displayed)
#' @param min.depth The minimum depth, positions with less depth than the threshold will be displayed in white on the plot.
#' @return A beautiful plot
#' @import ggplot2
#' @import dplyr
#' @import gtools
#' @import readxl
#' @examples
#' render.plot(path.data='/path/to/excel', sheet.name='final', min.frequency=10, min.depth=100)
#' @export
render.plot <- function(path.data, sheet.name, min.frequency = 0, min.depth = 100)
{
  data <- data.snp.process(path.data, sheet.name, min.depth)
  data <- data[data$Frequency>= min.frequency | is.na(data$Frequency),]
  p <- ggplot2::ggplot(data,ggplot2::aes(x=Region, y=Sample, color=Frequency)) + ggplot2::theme_classic() + ggplot2::geom_point()+ ggplot2::scale_color_gradient(na.value = "white", low="grey82", high="black")

  # Add ORF1ab
  for (i in 1:length(unique(data$Sample)))
  {
    # Add ORF1ab
    i = i-1
    p <- p+ ggplot2::annotate("rect", xmin=266-1, xmax=21555, ymin=0.65+i, ymax=1.35+i, alpha=0.5, fill="seagreen4")
    # Add S
    p <- p+ ggplot2::annotate("rect", xmin=21563-1, xmax=25384, ymin=0.65+i, ymax=1.35+i, alpha=0.5, fill="deeppink1")
    # Add ORF3a
    p <- p+ ggplot2::annotate("rect", xmin=25393-1, xmax=26220, ymin=0.65+i, ymax=1.35+i, alpha=0.5, fill="turquoise4")
    # Add E
    p <- p+ ggplot2::annotate("rect", xmin=26245-1, xmax=26472, ymin=0.65+i, ymax=1.35+i, alpha=0.5, fill="palevioletred2")
    # Add M
    p <- p+ ggplot2::annotate("rect", xmin=26523-1, xmax=27191, ymin=0.65+i, ymax=1.35+i, alpha=0.5, fill="lightsalmon2")
    # Add ORF6
    p <- p+ ggplot2::annotate("rect", xmin=27202-1, xmax=27387, ymin=0.65+i, ymax=1.35+i, alpha=0.5, fill="cadetblue3")
    # Add ORF7a
    p <- p+ ggplot2::annotate("rect", xmin=27394-1, xmax=27759, ymin=0.65+i, ymax=1.35+i, alpha=0.5, fill="darkorchid2")
    # Add ORF8
    p <- p+ ggplot2::annotate("rect", xmin=27894-1, xmax=28259, ymin=0.65+i, ymax=1.35+i, alpha=0.5, fill="violetred2")
    # Add N
    p <- p+ ggplot2::annotate("rect", xmin=28274-1, xmax=29533, ymin=0.65+i, ymax=1.35+i, alpha=0.5, fill="yellowgreen")
    # Add ORF10
    p <- p+ ggplot2::annotate("rect", xmin=29558-1, xmax=29674, ymin=0.65+i, ymax=1.35+i, alpha=0.5, fill="forestgreen")

  }

  p <- p + ggplot2::geom_point() + ggplot2::theme(  axis.ticks = ggplot2::element_blank(),
              axis.title= ggplot2::element_blank(),
              legend.position="none",
              axis.text.x = ggplot2::element_text(angle=90,size=6,face="bold",vjust = 0.5),
              axis.line.y = ggplot2::element_blank()) + ggplot2::xlim(min(data$Region)-1, 29674+1) +
    ggplot2::scale_x_continuous(breaks=data$Region,label=data$`Amino acid change in longest transcript`)+
    ggplot2::scale_y_discrete(breaks=unique(data$Sample)[gtools::mixedorder(unique(data$Sample))],
          limits=unique(data$Sample)[gtools::mixedorder(unique(data$Sample))])

  return(p)
}

