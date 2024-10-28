# Load libraries
library(tidyverse)
library(phytools)
library(ape)
library(ggtree)
library(ggnewscale)
library(colorRamps)
library(cowplot)
library(scatterpie)
library(vegan)
library(fitdistrplus)
library(stats)
library(scales)
library(ComplexHeatmap)

setwd('~/Desktop/my_pc/Vibrio/Github/Vibrio_2024/')

########
#Fig. 1#
########
# Read input matrices
meta_df1 <- read_csv('Input/Fig1/combined_meta.csv', col_names = TRUE)
meta_df2 <- read_csv('Input/Fig1/newset_meta.csv', col_names = TRUE)


# Read trees
tree1 <-
  read.tree('Input/Fig1/MLtree_whole_set.tre') %>%
  makeNodeLabel(method = 'number', prefix = '')

tree2 <- 
  read.tree('Input/Fig1/MLtree_new_set.tre') %>%
  drop.tip(c('1331650','D17181568','D17181396','D17181426','D17181545'))


# Set colour scales
lineage_col_scale <-
  scale_discrete_manual('colour',
                        values = c('grey', 'dodgerblue3',
                                   'maroon1', 'firebrick2', 'firebrick4',
                                   'purple'),
                        guide = guide_legend(order=1),
                        na.translate = FALSE)

geo_col_palette <- c(RColorBrewer::brewer.pal(5,'Set1'),'#F781BF','#A65628')
geo_col_scale <- scale_discrete_manual(
                'fill',
                values = geo_col_palette,
                name = 'Region',
                guide = guide_legend(order=3),
                na.translate = FALSE)

lineage_col_scale2 <-
  scale_discrete_manual('colour',
                        values = c('dodgerblue3',
                                   'firebrick2', 'firebrick4',
                                   'purple'),
                        guide = 'none')

loc_col_scale <- scale_discrete_manual("fill",values = geo_col_palette,
                                       breaks = c('1','2','4','5','10','15','17'))


# Plot
## P1_1 ##
hilight_df <- data.frame(id=c('137', '39'), `Sub-lineage`=c('BD-1.2', 'BD-1.1'))
tree_plot1 <- 
  ggtree(tree1, layout = 'rectangular', alpha = 0.7) %<+%
  meta_df1 +
  geom_tippoint(aes(colour=Lineage), show.legend = TRUE) +
  geom_cladelab(node='244', label='BD-1 this study', offset=0.008) +
  geom_cladelab(node='644', label='BD-2 this study', offset=0.007) +
  geom_cladelab(node='894', label='BD-3 this study', offset=0.031) +
  lineage_col_scale +
  geom_hilight(hilight_df, aes(node=id, fill=`Sub.lineage`), alpha=0.3) +
  scale_fill_manual(values = c('red', 'blue'), name = 'Sub-lineage', 
                    labels = c('BD-1.1', 'BD-1.2'),
                    guide = guide_legend(order = 2)) +
  new_scale_fill()

hmap_df1_1 <- data.frame('Region'=meta_df1$Region, 
                       row.names=meta_df1$Isolate)

hmap1_1 <- 
  gheatmap(tree_plot1,
           hmap_df1_1, width = 0.055, font.size = 3, 
           offset = 0.004,
           colnames = TRUE,
           colnames_position = "top",
           colnames_offset_y = 10,
           legend_title = "Region") + 
  ylim(0, 1100) +
  geo_col_scale
addi_hmap_col1 <- hmap1_1 + new_scale_fill()

hmap_df1_2 <- data.frame('Year'=as.numeric(meta_df1$Year), 
                         row.names=meta_df1$Isolate)

p1_1 <- 
  gheatmap(addi_hmap_col1,
           hmap_df1_2, width = 0.055, font.size = 3, 
           offset = 0.00825,
           colnames = TRUE,
           colnames_position = "top",
           colnames_offset_y = 10,
           legend_title = "Year") + 
  ylim(0, 1100) +
  scale_fill_continuous(low='yellow', high='blue', name='Year',
                        guide = guide_legend(order=4))

## P1_2 ##
hmap_df2_1 <- data.frame(`ICE_typing` = meta_df2$ICE_typing,
                         row.names = meta_df2$`Isolate`) %>%
  set_names(c('ICE'))

hmap_df2_2 <- data.frame(`Location` = factor(meta_df2$Location),
                         row.names = meta_df2$Isolate) %>%
  set_names(c('Location'))

tree_plot2 <- ggtree(tree2,
                     layout="rectangular",
                     alpha=0.7) %<+% meta_df2 +
  geom_tippoint(aes(colour=`Lineage`)) +
  lineage_col_scale2 +
  labs(colour="Lineage") 

hmap2_1 <- gheatmap(tree_plot2, 
                    hmap_df2_1, width=0.05, font.size=3, 
                    offset = 0,
                    colnames = TRUE,
                    colnames_position = "top",
                    colnames_offset_y = 1.1,
                    legend_title = "ICE") +
  scale_fill_viridis_d(option = "D", 'ICE', 
                       guide=guide_legend(order=2),
                       na.translate = FALSE)
addi_hmap_col2_1 <- hmap2_1 + new_scale_fill()

p1_2 <- gheatmap(addi_hmap_col2_1,
                    hmap_df2_2, width = 0.05, font.size=3,
                    offset = 0.009,
                    colnames = TRUE,
                    colnames_position = "top",
                    colnames_offset_y = 1.1) +
  loc_col_scale + labs(fill='Bangladesh site') +
  ylim(0,290) 

## Combine plots for fig. 1 ##
fig1 <- 
  cowplot::plot_grid(p1_1, p1_2,
                     ncol=1, nrow=2, labels='AUTO',
                     align = 'v')


########
#Fig. 2#
########
# Read input matrices
map_df <- read_csv('Input/Fig2/map_data.csv', col_names = TRUE)
shannon_loc_df <- read_csv('Input/Fig2/shannon_loc.csv', col_names = TRUE) %>%
  mutate(Location=factor(Location, levels = c('4','2','17','1','15','5','10')))
shannon_time_df <- read_csv('Input/Fig2/shannon_time.csv', col_names = TRUE)


# Set colour scales
map_col_scale <- scale_discrete_manual(
  "fill",
   values = c(RColorBrewer::brewer.pal(4,'Set1'),"gold"))


## P2_1 ##
world_data <- map_data('world')
Bangle_map <- ggplot(world_data, aes(x = long,
                                     y = lat))+
  geom_map(map=world_data,aes(map_id=region),fill='white',colour='black')+
  coord_quickmap()+
  coord_sf(xlim = c(88,93), ylim = c(20.5,27), expand = FALSE)

scatter_cols <- names(map_df)[4:8]

p2_1 <-
  Bangle_map+
  geom_scatterpie(aes(x=longitude,
                      y=latitude,
                      group=distinct_location,
                      r=0.2),
                  data=map_df,
                  cols=scatter_cols,
                  colour=NA,
                  alpha=1,
                  show.legend=TRUE,
                  legend_name = 'Lineage')+
  map_col_scale+
  geom_text(aes(x=longitude+0.3,
                y=latitude,
                label=distinct_location),
                data=map_df)+
  theme(plot.margin = margin(l=0, r=0)) +
  ylab('Latitude') +
  xlab('Longitude') +
  theme(text = element_text(size=16))

## P2_2 ##
p2_2 <- 
  ggplot(shannon_loc_df) +
  geom_col(aes(x=Location, y=Shannon)) +
  ylab('Shannon diversity index') +
  theme(aspect.ratio = 1, 
        plot.margin = margin(l=0, r=0)) +
  ylim(0,2.25) +
  theme(text = element_text(size=16),
        axis.text = element_text(size=11))

## P2_3 ##
p2_3 <- 
  ggplot(shannon_time_df) +
  geom_col(aes(x=Time, y=Shannon)) +
  ylab('Shannon diversity index') +
  theme(aspect.ratio = 1, 
        plot.margin = margin(l=0, r=0)) +
  ylim(0,2.25) +
  theme(text = element_text(size=16),
        axis.text = element_text(size=11))

## Combine plots for fig. 2 ##
fig2 <-
  cowplot::plot_grid(p2_1, NA, p2_2, p2_3,
                     nrow = 2, ncol = 2, 
                     labels = c('A', NA,'B','C'),
                     align = 'hv', axis = 'tblr',
                     label_size = 16)


########
#Fig. 3#
########
# Function
plot_pan_vs_genes <- function(input_table){
  df <- read.table(input_table,header = TRUE, sep = ',', check.names = FALSE) %>%
    dplyr::mutate(`Number of Isolates` = as.numeric(`Number of Isolates`)) %>%
    dplyr::mutate(`Pangenome Size` = as.numeric(`Pangenome Size`))
  
  group_by_size <- group_by(df, `Number of Isolates`) %>%
    summarise(mid_size = median(`Pangenome Size`))
  
  power_model <- stats::lm(log(group_by_size$`mid_size`)~log(group_by_size$`Number of Isolates`))
  
  g_value <- power_model[['coefficients']][['log(group_by_size$`Number of Isolates`)']] %>%
    scales::scientific(digits=3)
  
  my_label <- options('scipen'=999)
  my_label <- paste('gamma==',g_value)
  
  ggplot() + geom_point(df, mapping = aes(x=`Number of Isolates`, y=`Pangenome Size`),
                        size = 1, alpha = 0.7) +
    geom_smooth(group_by_size, mapping=aes(x=`Number of Isolates`, y=mid_size),
                colour = 'red') + 
    geom_text(aes(x=2/3*max(df$`Number of Isolates`), 
                  y=min(df$`Pangenome Size`)+1/4*(max(df$`Pangenome Size`)-min(df$`Pangenome Size`))),
              label=my_label, parse = TRUE) +
    scale_y_continuous(limits = c(3500, 4100)) +
    scale_x_continuous(limits = c(0, 60), labels = label_comma()) +
    theme(text=element_text(size = 16))
}


# Read input matrices
pan_bd1_pth <- 'Input/Fig3/pangenome_BD1.csv'
pan_bd2_pth <- 'Input/Fig3/pangenome_BD2.csv'
gene_dist_df <- read_csv('Input/Fig3/gene_dist.csv', col_names = TRUE)


## P3_1 ##
p3_1 <- cowplot::plot_grid(plot_pan_vs_genes(pan_bd1_pth),
                           plot_pan_vs_genes(pan_bd2_pth),
                           nrow = 1,
                           ncol = 2,
                           labels = c('BD-1', 'BD-2'),
                           label_size = 16)

## P3_2 ##
p3_2 <- 
  ggboxplot(gene_dist_df,
            x = 'core_dist_bin',
            y = 'accessory_dist',
            group = 'Lineages',
            fill = 'Lineages'
  ) + 
  stat_compare_means(aes(group=Lineages),
                     method = 't.test',
                     label = 'p.signif') +
  ylab('Accessory genome Jaccard distances') +
  xlab(expression(paste('Bin by core distances (ceil (', pi%*%1e+04, '))'))) +
  theme(text = element_text(size = 16))

## Combine plots for fig. 3 ##
fig3 <-
  cowplot::plot_grid(p3_1, p3_2,
                     nrow = 2, ncol = 1, labels = c('A', 'B'),
                     align = 'v', axis = 'tblr')


########
#Fig. 4#
########
# Read input matrices
amr_gene_df <- read_csv('Input/Fig4/card_rgi.csv')

amr_hmap_df <- 
  amr_gene_df %>%
  dplyr::select(-c('Isolate', 'Lineage')) %>%
  rename(c('vanY'='vanY gene in vanB cluster',
           'vanT'='vanT gene in vanG cluster',
           'parE'='Escherichia coli parE conferring resistance to fluoroquinolones',
           'varG'='Vibrio cholerae varG')) %>%
  mutate_all(function(x)
  {
    x=case_when(
      x == '1' ~ 'Yes',
      x == '0' ~ 'No')
  })

rownames(amr_hmap_df) <- amr_gene_df$Isolate


dfs_gene_df <- read_csv('Input/Fig4/dfs_system.csv')

dfs_hmap_df <-
  dfs_gene_df %>%
  dplyr::select(-c('Isolate', 'Lineage')) %>%
  mutate_all(function(x)
  {
    x=case_when(
      x == '1' ~ 'Yes',
      x == '0' ~ 'No')
  }) %>%
  rename(c('dCTP'='dCTPdeaminase', 'Lamassu'='Lamassu-Fam'))

rownames(hmap_df) <- dfs_df$Isolate


row_split <- c(rep('BD-1', 185), rep('BD-2a', 45),
               rep('BD-2b', 34), rep('BD-3', 9))


## p4_1 ##
p4_1 <- 
  Heatmap(amr_hmap_df, show_row_names = FALSE, 
          show_column_dend = FALSE, show_row_dend = FALSE,
          column_names_rot = 0, column_names_centered = TRUE,
          name = 'AMR gene presence          ',
          heatmap_legend_param = 
            list(labels_gp=gpar(fontsize=12),
                 title_gp=gpar(fontsize=12,fontface='bold')),
          col = c('black', 'grey'),
          rect_gp = gpar(col='grey', lwd=0.1),
          row_split = row_split)

## p4_2 ##
p4_2 <-
  Heatmap(dfs_hmap_df, show_row_names = FALSE, 
          show_column_dend = FALSE, show_row_dend = FALSE,
          column_names_rot = 0, column_names_centered = TRUE,
          name = 'Defence system presence',
          heatmap_legend_param = 
            list(labels_gp=gpar(fontsize=12),
                 title_gp=gpar(fontsize=12,fontface='bold')),
          col = c('black', 'grey'),
          rect_gp = gpar(col='grey', lwd=0.1),
          row_split = row_split)

## Combine plots for fig.4 ##
fig4 <-
  cowplot::plot_grid(grid.grabExpr(draw(p4_1)), 
                     grid.grabExpr(draw(p4_2)),
                     nrow = 2, ncol = 1,
                     labels = 'AUTO',
                     align = 'v')


## Save figures ##
ggsave(plot=fig1, filename='Plots/fig1.png', width=28, height=32, units='cm')
ggsave(plot=fig2, filename='Plots/fig2.png', width=36, height=36, units='cm')
ggsave(plot=fig3, filename='Plots/fig3.png', width=28, height=32, units='cm')
ggsave(plot=fig4, filename='Plots/fig4.png', width=30, height=30, units='cm')



