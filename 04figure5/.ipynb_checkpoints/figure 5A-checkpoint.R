parserOption <- function() {
    suppressMessages(library(argparse))
    DESC <- paste0(
        "This script is used to plot class count rect."
    )
    # Create parser object
    parser <- ArgumentParser(description = DESC)
    # By default ArgumentParser will add a help option
    parser$add_argument(
        "-i", "--infile",
        help = "The input file [forced, eg: data_assess_pre.txt]"
    )
    parser$add_argument(
        "-o", "--outdir",
        help = "The output directory [forced, eg: current directory]"
    )
    parser$add_argument(
        "-W", "--width",
        type = "double",
        help = "The width (inch) of the output device [optional, default: %(default)s]"
    )
    parser$add_argument(
        "-H", "--height",
        type = "double",
        help = "The height (inch) of the output device [optional, default: %(default)s]"
    )
    opt <- parser$parse_args()
    return(opt)
}

opt <- parserOption()

suppressMessages(library(tidyverse))

inDF <- read_tsv(opt$infile,
    col_names = T, col_types = cols(),
    locale = locale(encoding = "UTF-8")
)

classCountDF <- inDF %>%
    group_by(`Class I`) %>%
    summarise(count = n())

# plot
data <- classCountDF
colnames(data) <- c("Class", "num")
fraction <- data$`num` / sum(data$`num`)
num_per <- paste(round(fraction * 100, digits = 2), "%", sep = "")
data$name <- paste(data$`Class`, "(", num_per, ")")
data$per <- num_per
data$fraction <- fraction
data$ymax <- cumsum(data$fraction)
data$ymin <- c(0, head(data$ymax, n = -1))

g <- ggplot(data = data, aes(fill = name, ymax = ymax, ymin = ymin, xmax = 4, xmin = 3)) +
    geom_rect(alpha = 0.8, color = "white", size = 1) +
    theme_bw() +
    xlim(0, 5) +
    coord_polar(theta = "y") +
    labs(x = "", y = "", fill = "Class") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(panel.grid = element_blank()) +
    theme(axis.text = element_blank()) +
    theme(axis.ticks = element_blank()) +
    theme(panel.border = element_blank())

png(
    file = paste(opt$outdir, "/Class_Count_Doughnut", ".png", sep = ""),
    width = opt$width, height = opt$height, res = 300, units = "in"
)
print(g)
dev.off()
pdf(
    file = paste(opt$outdir, "/Class_Count_Doughnut", ".pdf", sep = ""),
    width = opt$width, height = opt$height
)
print(g)
dev.off()

