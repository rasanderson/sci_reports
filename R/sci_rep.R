# Rank abundance and other graphs for Chile Atacama data
library(tidyverse)
library(forcats)
library(cowplot)
library(vegan)
library(VennDiagram)
rm(list=ls())

rawd <-  read_csv("data/Atacamagenus.csv")

mean_LB <- apply(rawd[,28:34],1,mean)
LB_log_sort <- sort(log10(mean_LB+1), decreasing=TRUE)
mean_CAB2 <- apply(rawd[,35:38],1,mean)
CAB2_log_sort <- sort(log10(mean_CAB2+1), decreasing=TRUE)
mean_CAB3 <- apply(rawd[,39:40],1,mean)
CAB3_log_sort <- sort(log10(mean_CAB3+1), decreasing=TRUE)
mean_CHX1 <- apply(rawd[,53:56],1,mean)
CHX1_log_sort <- sort(log10(mean_CHX1+1), decreasing=TRUE)
mean_CHX2 <- apply(rawd[,57:60],1,mean)
CHX2_log_sort <- sort(log10(mean_CHX2+1), decreasing=TRUE)
mean_CHX3 <- apply(rawd[,61:65],1,mean)
CHX3_log_sort <- sort(log10(mean_CHX3+1), decreasing=TRUE)
mean_POP2 <- apply(rawd[,74:77],1,mean)
POP2_log_sort <- sort(log10(mean_POP2+1), decreasing=TRUE)
mean_VDL <- apply(rawd[,78:81],1,mean)
VDL_log_sort <- sort(log10(mean_VDL+1), decreasing=TRUE)
mean_CHX3 <- apply(rawd[,61:65],1,mean)
CHX3_log_sort <- sort(log10(mean_CHX3+1), decreasing=TRUE)
mean_Y2 <- apply(rawd[,82:85],1,mean)
Y2_log_sort <- sort(log10(mean_Y2+1), decreasing=TRUE)
mean_Y6_1 <- apply(rawd[,41:44],1,mean)
Y6_1_log_sort <- sort(log10(mean_Y6_1+1), decreasing=TRUE)
mean_Y6_2 <- apply(rawd[,45:48],1,mean)
Y6_2_log_sort <- sort(log10(mean_Y6_2+1), decreasing=TRUE)
mean_Y6_3 <- apply(rawd[,49:52],1,mean)
Y6_3_log_sort <- sort(log10(mean_Y6_3+1), decreasing=TRUE)


ranked_data <- tibble(rankno=1:391,
                      LB=LB_log_sort,
                      CAB2=CAB2_log_sort,
                      CAB3=CAB3_log_sort,
                      CHX1=CHX1_log_sort,
                      CHX2=CHX2_log_sort,
                      CHX3=CHX3_log_sort,
                      POP2=POP2_log_sort,
                      VDL=VDL_log_sort,
                      Y2=Y2_log_sort,
                      Y6_1=Y6_1_log_sort,
                      Y6_2=Y6_2_log_sort,
                      Y6_3=Y6_3_log_sort
                      )
ranked_lng <- ranked_data %>% gather(2:13, key="site", value="logabund")
site_order <- c("LB", "CAB2", "CAB3", "VDL", "Y6_1", "Y6_2", "Y6_3",
                "CHX1", "CHX2", "CHX3", "POP2", "Y2")
ranked_lng$site <- fct_relevel(ranked_lng$site, site_order)

ggplot(ranked_lng) +
   geom_line(aes(x=rankno, y=logabund, color=site), size=1.2) +
   labs(x="Abundance rank", y="Relative abundance (log)") +
   theme_classic()

mean_data <- tibble(rankno=1:391,
                      OTU=rawd$Name,
                      LB=mean_LB,
                      CAB2=mean_CAB2,
                      CAB3=mean_CAB3,
                      VDL=mean_VDL,
                      Y6_1=mean_Y6_1,
                      Y6_2=mean_Y6_2,
                      Y6_3=mean_Y6_3,
                      CHX1=mean_CHX1,
                      CHX2=mean_CHX2,
                      CHX3=mean_CHX3,
                      POP2=mean_POP2,
                      Y2=mean_Y2
)

mean_lng <- mean_data %>% gather(3:14, key="site", value="abund")

subset_OTU1 <- c("Nocardiaceae_uc",
                "Iamiaceae_uc",
                "Nocardioidaceae_uc",
                "Micromonosporaceae_uc",
                "HQ910322_f_uc",
                "Micrococcaceae_uc",
                "Microbacteriaceae_uc",
                "Geodermatophilaceae_uc",
                "Acidimicrobiaceae_uc",
                "FJ479147_f_uc")
subset_mean1 <- mean_lng %>% filter(OTU %in% subset_OTU1)
subset1 <- subset_mean1
subset1$site <- fct_relevel(subset1$site, site_order)

p <- ggplot(subset1) +
   geom_area(
      mapping=aes(x=site, y=abund, group=OTU, fill=OTU),
      position="fill") +
   labs(y = "Proportion", x="")
p <- add_sub(p, x=0.3, label="<-   Extreme hyper-arid location   ->", y=0, vjust=-2, vpadding = grid::unit(1, "lines"))
p <- add_sub(p, x=0.8, label="<-  Hyper-arid location  ->", y=1, vjust=-2, vpadding = grid::unit(0, "lines"))
ggdraw(p)   

subset_OTU2 <- c("Gordonia",
                 "Modestobacter",
                 "Microbacterium",
                 "Geodermatophilus",
                 "Verrucosispora",
                 "Arthrobacter",
                 "Blastococcus",
                 "HQ910322_g",
                 "HQ674860_g",
                 "FJ479147_g")
subset_mean2 <- mean_lng %>% filter(OTU %in% subset_OTU2)
subset2 <- subset_mean2
subset2$site  <- fct_relevel(subset2$site,  site_order)

p <- ggplot(subset2) +
   geom_area(
      mapping=aes(x=site, y=abund, group=OTU, fill=OTU),
      position="fill") +
   labs(y = "Proportion", x="")
p <- add_sub(p, x=0.3, label="<-   Extreme hyper-arid location   ->", y=0, vjust=-2, vpadding = grid::unit(1, "lines"))
p <- add_sub(p, x=0.8, label="<-  Hyper-arid location  ->", y=1, vjust=-2, vpadding = grid::unit(0, "lines"))
ggdraw(p)   


# Rarefaction
mean_data_df <- as.data.frame(mean_data)
row.names(mean_data_df) <- mean_data_df[,2]
mean_data_df <- mean_data_df[,-(1:2)]
mean_data_df <- mean_data_df[rowSums(mean_data_df>0),]
(raremax <- min(rowSums(mean_data_df)))
rareout <- rarecurve(round(t(mean_data_df)), step = 20, sample = raremax, col = "blue", cex = 0.6)
Nmax <- sapply(rareout, function(x) max(attr(x, "Subsample")))
Smax <- sapply(rareout, max)
rareout_lng <- NULL
for (i in seq_along(rareout)) {
   N <- attr(rareout[[i]], "Subsample")
   rareout_lng <- bind_rows(rareout_lng, tibble(Site=colnames(mean_data)[i+2],
                                                subsample=attr(rareout[[i]], "Subsample"),
                                                species=as.numeric(rareout[[i]])))
}
ggplot(rareout_lng) +
   geom_line(aes(x=subsample, y=species, color=Site), size=1.2) +
   labs(x="Number of sequences sampled", y="Number of OTUs observed")

venn.plot <- draw.quad.venn(
   area1 = 115,
   area2 = 127,
   area3 = 100,
   area4 = 86,
   n12 = 95,
   n13 = 68,
   n14 = 73,
   n23 = 71,
   n24 = 72,
   n34 = 65,
   n123 = 62,
   n124 = 66,
   n134 = 58,
   n234 = 59,
   n1234 = 55,
   category = c("Y6_1", "VDL", "CAB2", "LB"),
   fill = c("orange", "red", "green", "blue"),
   lwd = rep(1, 4),
   lty = "dashed",
   cex = 1,
   cat.cex = 1,
   ind=FALSE
)
grid.newpage()
grid.draw(venn.plot)

venn.plot <- draw.quad.venn(
   area1 = 123,
   area2 = 109,
   area3 = 128,
   area4 = 146,
   n12 = 83,
   n13 = 92,
   n14 = 93,
   n23 = 87,
   n24 = 93,
   n34 = 103,
   n123 = 77,
   n124 = 78,
   n134 = 83,
   n234 = 81,
   n1234 = 74,
   category = c("POP2", "CHX3", "CHX2", "CHX1"),
   fill = c("orange", "red", "green", "blue"),
   lwd = rep(1, 4),
   lty = "dashed",
   cex = 1,
   cat.cex = 1,
   ind=FALSE
)
grid.newpage()
grid.draw(venn.plot)

