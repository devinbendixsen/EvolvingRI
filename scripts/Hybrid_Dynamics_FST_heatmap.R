# ==============================================================================
# 
# ==============================================================================
H1_fst <- read.table(file='results/SNV/FST/H1_Fst.txt',row.names = NULL)
H1_fst <- separate_wider_delim(H1_fst, cols = row.names, delim = ";", names = c("P1", "P2"))
H1_fst <- H1_fst %>%
  mutate( across(
    .cols = everything(),
    ~str_replace( ., "N_Founder", "N_Founder_R1" )
  ) )

H1_fst$P1 <- str_replace_all(H1_fst$P1, "N_Founder", "N_Founder_R1")
H1_fst$P2 <- str_replace_all(H1_fst$P2, "N_Founder", "N_Founder_R1")
H1_fst$P2 <- str_replace_all(H1_fst$P2, "LE_Founder", "LE_Founder_R1")
H1_fst$P1 <- str_replace_all(H1_fst$P1, "LE_Founder", "LE_Founder_R1")

H1_fst <- separate_wider_delim(H1_fst, cols = P1, delim = "_", names = c("p1env", "p1gen",'p1rep'),cols_remove = FALSE)
H1_fst <- separate_wider_delim(H1_fst, cols = P2, delim = "_", names = c("p2env", "p2gen",'p2rep'),cols_remove = FALSE)
H1_fst_r <- H1_fst
H1_fst_r$P1_r <- H1_fst_r$P2
H1_fst_r$P2 <- H1_fst_r$P1
H1_fst_r$P1 <- H1_fst_r$P1_r
H1_fst_r = subset(H1_fst_r, select = -c(P1_r) )
H1_fst <- rbind(H1_fst,H1_fst_r)


H1_plot <- ggplot(H1_fst, aes(P2, P1, fill= Fst.Estimate)) + 
  geom_tile() +
  coord_fixed()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text=element_text(family="EB Garamond"),
        legend.position = "bottom")+
  scale_fill_distiller(palette = "YlGnBu")+
  labs(title='NaCl x LiAc 0.01', x='', y='')

# ==============================================================================
# 
# ==============================================================================
H2_fst <- read.table(file='results/SNV/FST/H2_Fst.txt',row.names = NULL)
H2_fst <- separate_wider_delim(H2_fst, cols = row.names, delim = ";", names = c("P1", "P2"))
H2_fst <- H2_fst %>%
  mutate( across(
    .cols = everything(),
    ~str_replace( ., "N_Founder", "N_Founder_R1" )
  ) )

H2_fst$P1 <- str_replace_all(H2_fst$P1, "N_Founder", "N_Founder_R1")
H2_fst$P2 <- str_replace_all(H2_fst$P2, "N_Founder", "N_Founder_R1")
H2_fst$P2 <- str_replace_all(H2_fst$P2, "LE_Founder", "LE_Founder_R1")
H2_fst$P1 <- str_replace_all(H2_fst$P1, "LE_Founder", "LE_Founder_R1")

H2_fst <- separate_wider_delim(H2_fst, cols = P1, delim = "_", names = c("p1env", "p1gen",'p1rep'),cols_remove = FALSE)
H2_fst <- separate_wider_delim(H2_fst, cols = P2, delim = "_", names = c("p2env", "p2gen",'p2rep'),cols_remove = FALSE)
H2_fst_r <- H2_fst
H2_fst_r$P1_r <- H2_fst_r$P2
H2_fst_r$P2 <- H2_fst_r$P1
H2_fst_r$P1 <- H2_fst_r$P1_r
H2_fst_r = subset(H2_fst_r, select = -c(P1_r) )
H2_fst <- rbind(H2_fst,H2_fst_r)


H2_plot <- ggplot(H2_fst, aes(P2, P1, fill= Fst.Estimate)) + 
  geom_tile() +
  coord_fixed()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text=element_text(family="EB Garamond"),
        legend.position = "bottom")+
  scale_fill_distiller(palette = "YlGnBu")+
  labs(title='NaCl x LiAc 0.02', x='', y='')
# ==============================================================================
# 
# ==============================================================================
H3_fst <- read.table(file='results/SNV/FST/H3_Fst.txt',row.names = NULL)
H3_fst <- separate_wider_delim(H3_fst, cols = row.names, delim = ";", names = c("P1", "P2"))
H3_fst <- H3_fst %>%
  mutate( across(
    .cols = everything(),
    ~str_replace( ., "N_Founder", "N_Founder_R1" )
  ) )

H3_fst$P1 <- str_replace_all(H3_fst$P1, "N_Founder", "N_Founder_R1")
H3_fst$P2 <- str_replace_all(H3_fst$P2, "N_Founder", "N_Founder_R1")
H3_fst$P2 <- str_replace_all(H3_fst$P2, "LE_Founder", "LE_Founder_R1")
H3_fst$P1 <- str_replace_all(H3_fst$P1, "LE_Founder", "LE_Founder_R1")

H3_fst <- separate_wider_delim(H3_fst, cols = P1, delim = "_", names = c("p1env", "p1gen",'p1rep'),cols_remove = FALSE)
H3_fst <- separate_wider_delim(H3_fst, cols = P2, delim = "_", names = c("p2env", "p2gen",'p2rep'),cols_remove = FALSE)
H3_fst_r <- H3_fst
H3_fst_r$P1_r <- H3_fst_r$P2
H3_fst_r$P2 <- H3_fst_r$P1
H3_fst_r$P1 <- H3_fst_r$P1_r
H3_fst_r = subset(H3_fst_r, select = -c(P1_r) )
H3_fst <- rbind(H3_fst,H3_fst_r)


H3_plot <- ggplot(H3_fst, aes(P2, P1, fill= Fst.Estimate)) + 
  geom_tile() +
  coord_fixed()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text=element_text(family="EB Garamond"),
        legend.position = "bottom")+
  scale_fill_distiller(palette = "YlGnBu")+
  labs(title='NaCl x Ethanol', x='', y='')

H1_plot + H2_plot + H3_plot
ggsave('figures/Hybrid_Dynamics_Fst.pdf',width=21,height=8,dpi = 900)
