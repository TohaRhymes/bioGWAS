# Load required libraries
library(ggplot2)
library(dplyr)
library(gridExtra)

get_only_legend <- function(plot) { 
  
  # get tabular interpretation of plot 
  plot_table <- ggplot_gtable(ggplot_build(plot))  
  
  #  Mark only legend in plot 
  legend_plot <- which(sapply(plot_table$grobs, function(x) x$name) == "guide-box")  
  
  # extract legend 
  legend <- plot_table$grobs[[legend_plot]] 
  
  # return legend 
  return(legend)  
}


setwd('./')

BARS_WIDTH=0.5
ERRORS_WIDTH=0.2



pastel_colors <- c("path_small" = "#b1decc", "path_medium" = "#a1cff0", 
                   "path_big" = "#9598f0", "path_random" = "#ebbcdb")


# Read the TPR data
data_tpr <- read.csv("data_enrich/TPR_to_draw.csv")
data_tpr$model <- recode(data_tpr$model, 
                     "linreg" = "MAGMA:\nlinreg", 
                     "mean" = "MAGMA:\nmean", 
                     "top" = "MAGMA:\ntop")
data_tpr$model <- factor(data_tpr$model, levels = c("MAGMA:\nlinreg", "MAGMA:\nmean", "MAGMA:\ntop", "PASCAL"))


# Read the FPR data
data_fpr <- read.csv("data_enrich/FPR_to_draw.csv")
data_fpr$model <- recode(data_fpr$model, 
                         "linreg" = "MAGMA:\nlinreg", 
                         "mean" = "MAGMA:\nmean", 
                         "top" = "MAGMA:\ntop")
data_fpr$model <- factor(data_fpr$model, levels = c("MAGMA:\nlinreg", "MAGMA:\nmean", "MAGMA:\ntop", "PASCAL"))

p1 <- ggplot(data_tpr, aes(x=model, y=score, fill=path)) +
  geom_bar(stat="identity", position=position_dodge(0.6), width=BARS_WIDTH, color="black", size=0.3) + 
  geom_point(position=position_dodge(0.6), size=1) +
  geom_errorbar(aes(ymin=min, ymax=max), position=position_dodge(0.6), size=0.5, width=ERRORS_WIDTH) +
  scale_fill_manual(values=pastel_colors) +
  labs(title="TPR", y="TPR", x="Model") +
  scale_y_continuous(limits=c(0,1), expand = c(0, 0)) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color="black", fill=NA, size=0.5),
    legend.position = "none"
  )

p2 <- ggplot(data_fpr, aes(x=model, y=score, fill=path)) +
  geom_bar(stat="identity", position=position_dodge(0.6), width=BARS_WIDTH/3, color="black", size=0.3) + 
  geom_point(position=position_dodge(0.6), size=1) +
  geom_errorbar(aes(ymin=min, ymax=max), position=position_dodge(0.6), size=0.5, width=ERRORS_WIDTH/3) +
  scale_fill_manual(values=pastel_colors) +
  labs(title="FPR", y="FPR", x="Model") +
  scale_y_continuous(limits=c(0, 1), expand = c(0, 0)) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color="black", fill=NA, size=0.5),
    legend.position = "none"
  )


p_merged<-grid.arrange(p1, p2, ncol=2)
p_merged_legend <- grid.arrange(get_only_legend(p1+theme(legend.position = "bottom")),
                    get_only_legend(p2+theme(legend.position = "bottom")))
p_common <- grid.arrange(p_merged, p_merged_legend, nrow = 2, heights = c(10, 1))
p_common
ggsave(filename = "images/tpr_fpr.pdf", plot = p_common, width = 10, height = 6, device = "pdf", path = "./")
