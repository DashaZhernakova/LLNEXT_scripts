args <- commandArgs(trailingOnly = TRUE)

library(ggplot2)
library(tidyr)
library(dplyr)
library(ggrepel)

infile = args[1]

d <- read.delim(infile, sep = "\t", as.is = T, check.names = F, header = F)

colnames(d) <- c("CHROM","POS","ID","REF","ALT","A1","TEST","OBS_CT","BETA","SE","T_STAT","P")
d$POS <- as.numeric(d$POS)
d$P <- as.numeric(d$P)


don <- d %>% 
  
  # Compute chromosome size
  group_by(CHROM) %>% 
  summarise(chr_len=max(POS)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(d, ., by=c("CHROM"="CHROM")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHROM, POS) %>%
  mutate( BPcum=POS+tot)


axisdf <- don %>% group_by(CHROM) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png(paste0(infile,".manhattan_noannot.png"), width = 15, height = 5, units = 'in', res = 400)
ggplot(don, aes(x=BPcum, y=-log10(P))) +
  # Show all points
  geom_point( aes(color=as.factor(CHROM)),  size=1) +
  scale_color_manual(values = rep(c("#495DA0", "#74AFDF"), 22 )) +
  geom_hline(yintercept=-1*log(1.79e-09,10), color = "#EF3B2C") +
  # custom X axis:
  xlab("Chromosome") +
  scale_x_continuous( label = axisdf$CHROM, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) , limits = c(0,45)) +
  
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line.y = element_line(color="lightgrey", size = 0.5),
    axis.ticks.x = element_blank(),
    axis.text = element_text(size=12)
    
  )
dev.off()