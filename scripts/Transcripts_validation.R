#Validate detected transcripts with CAGE-seq, 3'-seq, and splice junctions data


library(ggplot2)
library(ggtranscript)
library(here)
library(cowplot)
library(dplyr)


valid_gtf <- import("E:/Encode_annotation/all/locus_select_CLN3_my_id_filtered.gtf") %>% 
  as_tibble()
#valid_gtf <- my_id_gtf_filtered

trans <- makeGRangesFromDataFrame(valid_gtf,keep.extra.columns = T)



get_locus <- function(seqnames, start, end, strand, gene_id) {
  
  locus_subset <- 
    GRanges(seqnames = seqnames,
            ranges = IRanges(start = start, 
                             end = end),
            strand = strand)
  
}



#zoomed in for top 15
cln3_locus <- get_locus("chr16", 28474111,28495575,"-", "ENSG00000188603" )
CAGE_seq <- 
  rtracklayer::import.bed("E:/FANTOM5/refTSS_v3.3_human_coordinate.hg38.bed.gz") %>% 
  subsetByOverlaps(., cln3_locus)


CAGE_plot <-
  CAGE_seq %>% 
  data.frame() %>% 
  ggplot(aes(xmin = start, xmax = end, ymin = 0, ymax = 1)) +
  geom_rect(colour = "black", fill = "black", show.legend = F, alpha=0.8) +
  xlim(start(cln3_locus), end(cln3_locus)) +
  #labs(y = "",#"CAGE peaks",
  #     x = "") +
  theme_classic() +
  theme(#axis.title = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_line(colour = "white"),
        axis.title.y = element_text(angle = 1),
        #strip.text.y = element_text(face = "bold",
        #                            size = 10),
        strip.background = element_rect(fill ="lightgrey"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.spacing = unit(0, "cm"))





#check distance of TSS to CAGE peaks
full_trans <- import.gff3("E:/Encode_annotation/all/locus_select_cln3_my_id_filtered.gff3")

trans_tibble <- as_tibble(full_trans)


dist_to_CAGE <- function(data, CAGE) {
  
  data_transcript <- 
    data %>% plyranges::filter(type == "transcript")
  
  table <- data.frame(Transcript = data_transcript$ID,
                      Dist_to_CAGE = as.numeric(0))
  
  for (i in table$Transcript) {
    
    gr <-
      data_transcript %>% 
      plyranges::filter(ID == i)
    
    if (as.character(strand(gr)) == "-") {
      
      start(gr) <- end(gr)
      GenomicRanges::distance(gr, CAGE)
      table$Dist_to_CAGE[table$Transcript == i] <- min(GenomicRanges::distance(gr, CAGE))
      
    }else {
      
      end(gr) <- start(gr)
      GenomicRanges::distance(gr, CAGE) 
      table$Dist_to_CAGE[table$Transcript == i] <- min(GenomicRanges::distance(gr, CAGE))
    } 
  }
  
  return(table)
}



test <- dist_to_CAGE(full_trans, CAGE_seq) 

test_50 <- test %>% filter(test$Dist_to_CAGE <= 50)



# 3'PAS seq


pas_seq <- read_tsv(file = "E:/PAS_seq/polyAsite_chr16.tsv", col_names = T)


colnames(pas_seq)[1] <- "chr"
colnames(pas_seq)[2] <- "start"
colnames(pas_seq)[3] <- "end"



pas_seq$chr <- "chr16"

pas_seq <- makeGRangesFromDataFrame(pas_seq)

pas_seq <- pas_seq %>% 
  subsetByOverlaps(., cln3_locus)



dist_to_PAS <- function(data, PAS) {
  
  data_transcript <- 
    data %>% plyranges::filter(type == "transcript")
  
  table <- data.frame(Transcript = data_transcript$ID,
                      Dist_to_PAS = as.numeric(0))
  
  for (i in table$Transcript) {
    
    gr <-
      data_transcript %>% 
      plyranges::filter(ID == i)
    
    if (as.character(strand(gr)) == "-") {
      
       end(gr) <- start(gr)
      GenomicRanges::distance(gr, PAS)
      table$Dist_to_PAS[table$Transcript == i] <- min(GenomicRanges::distance(gr, PAS))
      
    }else {
      
       start(gr) <- end(gr)
      GenomicRanges::distance(gr, PAS) 
      table$Dist_to_PAS[table$Transcript == i] <- min(GenomicRanges::distance(gr, PAS))
    } 
  }
  
  return(table)
}


test_pas <- dist_to_PAS(full_trans, pas_seq) 



pas_plot <-
  pas_seq %>% 
  data.frame() %>% 
  ggplot(aes(xmin = start, xmax = end, ymin = 0, ymax = 1)) +
  geom_rect(colour = "black", fill = "black", show.legend = F, alpha=0.8) +
  xlim(start(cln3_locus), end(cln3_locus)) +
  #labs(y = "",   #"PolyA sites",
  #     x = "") +
  theme_classic() +
  theme(#axis.title = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_line(colour = "white"),
        axis.title.y = element_text(angle = 1),
        #strip.text.y = element_text(face = "bold",
         #                           size = 10),
        strip.background = element_rect(fill ="lightgrey"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.spacing = unit(0, "cm"))

#pas_plot





###splicing junctions

junc <- read.csv("E:/splicing junctions/junction.csv")


rjun <- junc %>% separate(Junction.location, c(NA,'coord', NA), sep = ":")

junc$Junction.location <- sub(":", "/", junc$Junction.location)
junc <- junc %>% separate(Junction.location, c("seqnames", "position"), sep = "/")
junc <- junc %>% separate(position, c("coord", "strand"), sep = ":") 


junc <- junc %>% separate(coord, c("start", "end"), sep = "\\|")


#check all junctions

export(my_id_gtf_filtered, "E:/splicing junctions/my_id_filtered.gtf", format = "gtf")
full_gtf <- gread::read_format("E:/splicing junctions/my_id_filtered.gtf")
full_introns <- gread::construct_introns(full_gtf) %>% as_tibble()

full_junc <- full_introns %>% filter(feature == "intron") %>% dplyr::select(start,end) %>% unique()

#get junction coordinates from intron coordinates
full_junc$start <- full_junc$start - 1
full_junc$end <- full_junc$end + 1


full_junc$junc <- paste(full_junc$start, full_junc$end, sep = "|")

full_in_rjun <- full_junc %>% filter(full_junc$junc %in% rjun$coord)

rjun_in_encode <- rjun %>% filter(coord %in% full_in_rjun$junc)


#use SpliceVault data
splicevault <- read.table(file="E:/splicing junctions/splicevault/splicevault_subset.txt", sep = "\t")
colnames(splicevault) <- c("splice_site_pos", "ss_type", "exon_no", "splicing_event_class", "event_rank",
                           "in_gtex", "in_sra", "missplicing_inframe", "gtex_sample_count", "max_uniq_reads", 
                           "sra_sample_count", "sample_count", "skipped_exons_count",
                           "skipped_exons_id", "cryptic_distance", "chr", "donor_pos", "acceptor_pos", "gene_tx_id")
cln3_known <- read.table(file = "E:/splicing junctions/splicevault/cln3_rows.txt", header = F, sep = "\t")

#Keep only CLN3 junctions using tx ids
splicevault_cln3 <- splicevault %>% filter(splicevault$gene_tx_id %in% cln3_known$V1)

splicevault_cln3$coord <- paste(splicevault_cln3$acceptor_pos - 1, "|", splicevault_cln3$donor_pos + 1, sep = "")

junc_not_in_rjun <- unique(full_not_in_rjun$junc)


#For 17 absent junctions
detected <- splicevault_cln3 %>% filter(coord %in% junc_not_in_rjun)



full_junc_vector <- unique(full_junc$junc)
detected_full <- splicevault_cln3 %>% filter(coord %in% full_junc_vector)

full_cryptic <- detected_full %>% filter(splicing_event_class %in% c("cryptic donor", "cryptic acceptor")) %>%
  #dplyr::select(splicing_event_class, coord) %>%
  unique()


full_cryptic$transcript_id <- cln3_known$V3[match(full_cryptic$gene_tx_id, cln3_known$V1)]

full_cryptic <- full_cryptic %>% dplyr::select(-coord) %>% unique()
