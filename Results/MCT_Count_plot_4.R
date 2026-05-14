setwd("../2_results/Results_Simulation_Excel/")

pacman::p_load(ggplot2, writexl, readxl, dplyr, tidyr, purrr, janitor, 
               forcats, gridExtra, patchwork, cowplot)

## Type-I Error
TI4 <- read_excel("MCT_Count_typeI_4groups.xlsx")
## Power
P4 <- read_excel("MCT_Count_power_4groups.xlsx")

## change directory to the folder for plots
setwd("../Results_Plots/")

############################# Type-I Error ####################################
TI4_long <- TI4 %>%  rename(Type_I_Error = global) %>%
  mutate(Method = fct_relevel(Method,  "boot", "hom", "het", "het_log", "het_sqrt", 
                              "poi", "nb", "q_poi"))

#################################################################################
labels <- c("(1) POI balanced homoscedastic", 
            "(2) POI unblanced homoscedastic",
            expression(paste("(3) CMP balanced homoscedastic ", italic(nu)[1],"=",italic(nu)[2],"=",italic(nu)[3],"=",italic(nu)[4],"=0.5")),
            expression(paste("(4) CMP balanced homoscedastic ", italic(nu)[1],"=",italic(nu)[2],"=",italic(nu)[3],"=",italic(nu)[4],"=2")), 
            expression(paste("(5) CMP balanced heterosc. ", italic(nu)[1],"=0.2, ", italic(nu)[2],"=0.2, ", italic(nu)[3],"=0.5, ", italic(nu)[4],"=0.5")), 
            expression(paste("(6) CMP balanced heteroscedastic ", italic(nu)[1],"=0.2, ", italic(nu)[2],"=0.5, ", italic(nu)[3],"=2, ", italic(nu)[4],"=2")),
            expression(paste("(7) CMP positive pairing ", italic(nu)[1],"=0.5, ", italic(nu)[2],"=0.5, ", italic(nu)[3],"=0.2, ", italic(nu)[4],"=0.2")),
            expression(paste("(8) CMP positive pairing ", italic(nu)[1],"=2, ", italic(nu)[2],"=2, ", italic(nu)[3],"=0.5, ", italic(nu)[4],"=0.2")),
            expression(paste("(9) CMP negative pairing ", italic(nu)[1],"=0.2, ", italic(nu)[2],"=0.2, ", italic(nu)[3],"=0.5, ", italic(nu)[4],"=0.5")),
            expression(paste("(10) CMP negative pairing ", italic(nu)[1],"=0.2, ", italic(nu)[2],"=0.5, ", italic(nu)[3],"=2, ", italic(nu)[4],"=2")),
            expression(paste("(11) NB balanced homoscedastic ", size[1],"=",size[2],"=",size[3],"=",size[4],"=3")), 
            expression(paste("(12) NB positive pairing ", size[1],"=5, ", size[2],"=3, ", size[3],"=3, ", size[4],"=2")),
            expression(paste("(13) NB negative pairing ", size[1],"=0.75, ", size[2],"=2, ", size[3],"=2, ", size[4],"=3")),
            expression(paste("(14) NB negative pairing ", size[1],"=2, ", size[2],"=3, ", size[3],"=3, ", size[4],"=5")))

plots <- map(as.list(unique(TI4_long$no)),
                 function(b) {
                   TI4_long %>% filter(no == b) %>% 
                     ggplot(aes(y = Type_I_Error,  x = m, fill = Method, color = Method, linetype=Method)) +
                     geom_line()+
                     scale_y_continuous(limits = c(0,0.25), breaks = c(seq(0,0.25,0.05)))+
                     labs(y = "Type I Error",
                          subtitle = labels[b]) +
                     geom_hline(yintercept = 0.05, lwd = 0.5, color = "black") +
                     theme_minimal(base_size = 12) +
                     theme(legend.position = "bottom",
                           plot.subtitle = element_text(hjust = 0.5)) +
                     guides(color = guide_legend(nrow = 1)) +
                     facet_grid(lambda ~ Contrast, labeller = label_both)
                 })  

### for plotting individual plots
#for (i in 1:14) {
#  png(filename = paste("Plot_",i,".png"), 
#       width = 2000, height = 2000, res = 300)
#  plot(plots_[[i]])
#  dev.off()
#}

################################################################################
########################## combine Graphics ####################################
plot0 <- plots[[1]]
g <- ggplotGrob(plot0)
leg <- g$grobs[which(sapply(g$grobs, function(x) x$name) == "guide-box")][[1]]

### POI --- Figure A1 ---
p1 <- plots[[1]]+theme(legend.position = "none")
p2 <- plots[[2]]+theme(legend.position = "none")   

png(filename = "Plot_POI_4.png",
    width = 3000, height = 1600, res = 300)
(p1 | p2 ) / 
  leg + plot_layout(heights = c(1,0.1))
dev.off()

### CMP
# 1 --- Figure 1 ---
p1 <- plots[[3]]+theme(legend.position = "none")
p2 <- plots[[4]]+theme(legend.position = "none")   
p3 <- plots[[5]]+theme(legend.position = "none")
p4 <- plots[[6]]+theme(legend.position = "none")   

png(filename = "Plot_CMP1_4.png",
    width = 3000, height = 3200, res = 300)
(p1 | p3 ) / 
(p2 | p4) /
  leg + plot_layout(heights = c(1,1,0.1))
dev.off()

# 2 --- Figure 2 ---
p1 <- plots[[7]]+theme(legend.position = "none")
p2 <- plots[[8]]+theme(legend.position = "none")   
p3 <- plots[[9]]+theme(legend.position = "none")
p4 <- plots[[10]]+theme(legend.position = "none")   

png(filename = "Plot_CMP2_4.png",
    width = 3000, height = 3200, res = 300)
(p1 | p3 ) /
(p2 | p4) /
  leg + plot_layout(heights = c(1,1,0.1))
dev.off()

### NB --- Figure A2 ---
p1 <- plots[[11]]+theme(legend.position = "none")
p2 <- plots[[12]]+theme(legend.position = "none")   
p3 <- plots[[13]]+theme(legend.position = "none")
p4 <- plots[[14]]+theme(legend.position = "none")   

png(filename = "Plot_NB_4.png",
    width = 3000, height = 3200, res = 300)
(p1 | p3 ) / 
(p2 | p4) /
  leg + plot_layout(heights = c(1,1,0.1))
dev.off()


################################################################################
############################ Power #############################################
################################################################################
P4_long <- P4 %>%  rename(Power = global) %>%
  filter(Method!="poi" & Method!="het_log" & Method!="het_sqrt") %>%
  mutate(Method = fct_relevel(Method,  "boot", "hom", "het", "nb", "q_poi"))


################################################################################
#default <- hue_pal()(8)
#show_col(default)
custom_colors <- c("#F8766D", "#CD9600", "#7CAE00", "#C77CFF", "#FF61CC")
custom_lines <- c(1,3,2,5,6)


plots <- map(as.list(unique(P4_long$no)),
                 function(b) {
                   P4_long %>% filter(no == b) %>% 
                     ggplot(aes(y = Power,  x = delta1, fill = Method, color = Method, linetype=Method)) +
                     geom_line()+
                     scale_y_continuous(limits = c(0,1), breaks = c(seq(0,1,0.2)))+
                     scale_x_continuous(breaks = c(seq(0,3,0.5)))+
                     scale_color_manual(values = custom_colors)+
                     scale_linetype_manual(values = custom_lines)+
                     labs(y = "Power",
                          x = expression(italic(delta)),
                          subtitle = labels[[b]]) +
                     geom_hline(yintercept = 0.8, lwd = 0.5, color = "black") +
                     theme_minimal(base_size = 12) +
                     theme(legend.position = "bottom",
                           plot.subtitle = element_text(hjust = 0.5)) +
                     guides(color = guide_legend(nrow = 1)) +
                     facet_grid(lambda ~ Contrast, labeller = label_both) 
                 })

################################################################################
########################## combine Graphics ####################################
plot0 <- plots[[1]]
g <- ggplotGrob(plot0)
leg <- g$grobs[which(sapply(g$grobs, function(x) x$name) == "guide-box")][[1]]

### POI --- Figure A3 ---
p1 <- plots[[1]]+theme(legend.position = "none")
p2 <- plots[[2]]+theme(legend.position = "none")   

png(filename = "PlotPower_POI_4.png",
    width = 3000, height = 1600, res = 300)
(p1 | p2 ) / 
  leg + plot_layout(heights = c(1,0.1))
dev.off()

### CMP
# 1 --- Figure A4 ---
p1 <- plots[[3]]+theme(legend.position = "none")
p2 <- plots[[4]]+theme(legend.position = "none")   
p3 <- plots[[5]]+theme(legend.position = "none")
p4 <- plots[[6]]+theme(legend.position = "none")   

png(filename = "PlotPower_CMP1_4.png",
    width = 3000, height = 3200, res = 300)
(p1 | p3 ) / 
(p2 | p4) /
  leg + plot_layout(heights = c(1,1,0.1))
dev.off()

# 2 --- Figure 3 ---
p1 <- plots[[7]]+theme(legend.position = "none")
p2 <- plots[[8]]+theme(legend.position = "none")   
p3 <- plots[[9]]+theme(legend.position = "none")
p4 <- plots[[10]]+theme(legend.position = "none")   

png(filename = "PlotPower_CMP2_4.png",
    width = 3000, height = 3200, res = 300)
(p1 | p3 ) /
(p2 | p4) /
  leg + plot_layout(heights = c(1,1,0.1))
dev.off()

### NB --- Figure A5 ---
p1 <- plots[[11]]+theme(legend.position = "none")
p2 <- plots[[12]]+theme(legend.position = "none")   
p3 <- plots[[13]]+theme(legend.position = "none")
p4 <- plots[[14]]+theme(legend.position = "none")   

png(filename = "PlotPower_NB_4.png",
    width = 3000, height = 3200, res = 300)
(p1 | p3 ) / 
(p2 | p4) /
  leg + plot_layout(heights = c(1,1,0.1))
dev.off()
