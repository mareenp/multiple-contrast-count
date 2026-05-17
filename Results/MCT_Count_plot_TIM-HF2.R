pacman::p_load(ggplot2, writexl, readxl, dplyr, tidyr, purrr, janitor, 
               forcats, gridExtra, patchwork, cowplot)

## Type-I Error
TI4 <- read_excel("MCT_Count_typeI_TIM-HF2.xlsx")

############################# Type-I Error ####################################
TI4_long <- TI4 %>%  rename(Type_I_Error = global) %>%
  mutate(Method = fct_relevel(Method,  "boot", "hom", "het", "het_log", "het_sqrt", 
                              "poi", "nb", "q_poi"))

#################################################################################
plots <- map(as.list(unique(TI4_long$no)),
                 function(b) {
                   TI4_long %>% filter(no == b) %>% 
                     ggplot(aes(y = Type_I_Error,  x = n1, fill = Method, color = Method, linetype=Method)) +
                     geom_line()+
                     scale_y_continuous(limits = c(0,0.25), breaks = c(seq(0,0.25,0.05)))+
                     scale_x_continuous(limits = c(60, 400), breaks = c(60, seq(100,400, 50)))+
                     labs(y = "Type I Error",
                          x = expression(paste("Sample size ", n[i], " per group"))) + 
                     geom_hline(yintercept = 0.05, lwd = 0.5, color = "black") +
                     theme_minimal(base_size = 12) +
                     theme(legend.position = "bottom",
                           plot.subtitle = element_text(hjust = 0.5)) +
                     guides(color = guide_legend(nrow = 1))
                 })  


png(filename = "Plot_sim_TIM-HF2.png",
    width = 3000, height = 3000, res = 660)
plots[[1]]
dev.off()

