######## This R script works best on Mac OS ######## 

pacman::p_load(tidyr, dplyr, stringr, ggplot2, forcats, ggtext)

#INPUT = snakemake@input[[1]]
INPUT = "../all.combined.tsv"
KEEP = c("MeV", "E2", "HCoV_OC43", "FLUAV", "DENV-1","HAdV-5", "HuAHV3; VZV", "MPV")

# first clean up - filter to whats needed
dat = read.csv(INPUT, sep = "\t") %>%
    # keep only hits with a consensus
    filter(!is.na(consensus_completeness)) %>%
    # keep only viruses spiked in
    filter(name %in% KEEP)

# second clean up - renaming samplename
dat2 = dat %>%
  separate(col = id, into = c("sample_name", "replicate"), sep = "-") %>%
  filter(!grepl("[cC]ontrol", dat$id))

# split sample_name into needed fields
dat3 = dat2 %>%
  separate(col = sample_name, into = c("molecule", "matrix", "dilution", "subsample_limit_short"), sep = "_", remove = F, extra = "merge") %>%
  mutate(matrix = case_when(matrix == "Comp" ~ "Composite", 
                            matrix == "NFW" ~ "Nuclease-Free Water", 
                            matrix == "PS" ~ "Passive Sampler", 
                            TRUE ~ matrix)) %>%
  mutate(name = case_when(name == "FLUAV" ~ "IAV", 
                            name == "HCoV_OC43" ~ "OC43", 
                            name == "DENV-1" ~ "DENV", 
                            name == "E2" ~ "EV",
                            name == "HAdV-5" ~ "AdV",
                            name == "HuAHV3; VZV" ~ "VZV",
                            name == "MPV" ~ "MPXV",
                            TRUE ~ name)
         ) %>%
  mutate(segment = case_when(segment == "RNA1" ~ "PB2", 
                             segment == "RNA2" ~ "PB1",
                             segment == "RNA3" ~ "PA",
                             segment == "RNA4" ~ "HA",
                             segment == "RNA5" ~ "NP",
                             segment == "RNA6" ~ "NA",
                             segment == "RNA7" ~ "MP",
                             segment == "RNA8" ~ "NS",
                             TRUE ~ segment)
        ) %>%
  mutate(replicate = gsub("replicate", "", replicate)) %>%
  mutate(name_segments = case_when(
    name == "IAV" ~ paste(name, segment, sep ="-"), 
    TRUE ~ name)
    )

# flu HA NA MP segments stats analysis  
flu_stats = dat3 %>% 
  filter(name == "IAV") %>%
  filter(segment %in% c("HA", "NA", "MP")) %>%
  group_by(molecule, matrix, dilution, subsample_limit_short, name_segments) %>%
  mutate(
    std_cc = sd(consensus_completeness),
    mean_cc = mean(consensus_completeness),
    lower_cc = mean_cc - std_cc,
    upper_cc = mean_cc + std_cc,
    std_rc = sd(remap_coverage),
    mean_rc = mean(remap_coverage),
    lower_rc = mean_rc - std_rc,
    upper_rc = mean_rc + std_rc,
  ) %>%
  ungroup %>%
  # refactor levels for plotting purposes
  mutate(subsample_limit_short = fct_relevel(subsample_limit_short, c("50k","100k","200k","300k","500k","800k","1M","1.5M","2M","3M"))) %>%
  # now that we calculate stats, we can drop replicates
  distinct(molecule, matrix, dilution, subsample_limit_short, name_segments, segment, .keep_all = T)


# flu all internal segment stats analysis
flu_segments_stats = dat3 %>% 
  filter(name == "IAV") %>%
  group_by(molecule, matrix, dilution, subsample_limit_short, name, segment) %>%
  mutate(
    std_cc = sd(consensus_completeness),
    mean_cc = mean(consensus_completeness),
    lower_cc = mean_cc - std_cc,
    upper_cc = mean_cc + std_cc,
    std_rc = sd(remap_coverage),
    mean_rc = mean(remap_coverage),
    lower_rc = mean_rc - std_rc,
    upper_rc = mean_rc + std_rc,
  ) %>%
  ungroup %>%
  # refactor levels for plotting purposes
  mutate(subsample_limit_short = fct_relevel(subsample_limit_short, c("50k","100k","200k","300k","500k","800k","1M","1.5M","2M","3M"))) %>%
  # now that we calculate stats, we can drop replicates
  distinct(molecule, matrix, dilution, subsample_limit_short, name, segment, .keep_all = T)

# all other stats analysis
allother_stats = dat3 %>%
  filter(name != "IAV") %>%
  group_by(molecule, matrix, dilution, subsample_limit_short, name_segments) %>%
  mutate(
    std_cc = sd(consensus_completeness),
    mean_cc = mean(consensus_completeness),
    lower_cc = mean_cc - std_cc,
    upper_cc = mean_cc + std_cc,
    std_rc = sd(remap_coverage),
    mean_rc = mean(remap_coverage),
    lower_rc = mean_rc - std_rc,
    upper_rc = mean_rc + std_rc,
  ) %>%
  ungroup %>%
  # refactor levels for plotting purposes
  mutate(subsample_limit_short = fct_relevel(subsample_limit_short, c("50k","100k","200k","300k","500k","800k","1M","1.5M","2M","3M"))) %>%
  # now that we calculate stats, we can drop replicates
  distinct(molecule, matrix, dilution, subsample_limit_short, name, .keep_all = T)

# combine datasets
stats = bind_rows(allother_stats, flu_stats) %>%
  mutate(matrix = fct_relevel(matrix, "Nuclease-Free Water", "Composite", "Grab", "Passive Sampler")) %>%
  mutate(name_segments = fct_relevel(name_segments, "OC43", "DENV", "EV", "IAV-HA", "IAV-NA", "IAV-MP", "MeV", "AdV", "VZV", "MPXV"))

#set up mytheme for ggplot
my_colors_new = c('#e31a1c' = "1e106",
              '#1f78b4' = "2e105",
              '#33a02c' = "4e104",
              '#ff7f00' = "8e103")

mytheme = list(
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.text = element_markdown()
    ),
  scale_color_manual(
      values = names(my_colors_new), 
      breaks = my_colors_new,
      labels = c("1e106" = "1x10<sup>6</sup>", "2e105" = "2x10<sup>5</sup>", "4e104" = "4x10<sup>4</sup>", "8e103" = "8x10<sup>3</sup>" )
    ),
  scale_shape_manual(
      guide = guide_legend(override.aes = list(size = 5, alpha = 1, linetype = 0)),
      values = c(15,16,17,18),
      labels = c("1e106" = "1x10<sup>6</sup>", "2e105" = "2x10<sup>5</sup>", "4e104" = "4x10<sup>4</sup>", "8e103" = "8x10<sup>3</sup>")
    )
)

#ggplot 
plot_cc = stats %>%
  ggplot(aes(
    x = subsample_limit_short,
    y = mean_cc,
    group = dilution,
    color = dilution
  )) +
  geom_line(alpha=0.5) +
  geom_point(
    aes(shape = dilution), 
    size = 2, 
    alpha = 0.5) +
  geom_errorbar(aes(ymin = lower_cc, ymax = upper_cc), width = 0.1) +
  facet_grid(name_segments ~ matrix) +
  ylab("Consensus Completeness (%)") +
  xlab("Read Subsampling") +
  mytheme
plot_cc

ggsave(filename = "subsamp_cc.png", plot = plot_cc, device = "png", width = 297, height = 210, units = "mm")


plot_rc = stats %>%
  ggplot(aes(
    x = subsample_limit_short,
    y = mean_rc,
    group = dilution,
    color = dilution
  )) +
  geom_line(alpha=0.5) +
  geom_point(
    aes(shape = dilution), 
    size = 2, 
    alpha = 0.5) +
  geom_errorbar(aes(ymin = lower_rc, ymax = upper_rc), width = 0.1) +
  facet_grid(name_segments ~ matrix) +
  ylab("Genome Coverage (%)") +
  xlab("Read Subsampling") +
  mytheme
plot_rc

ggsave(filename = "subsamp_rc.png", plot = plot_cc, device = "png", width = 297, height = 210, units = "mm")

plot_flu_segments = flu_segments_stats %>%
  ggplot(aes(
    x = subsample_limit_short, 
    y = mean_cc,
    group = dilution,
    color = dilution
  )) +
  geom_line(alpha=0.5) + 
  geom_point(
    aes(shape = dilution),
    size = 2, 
    alpha = 0.5) +
  geom_errorbar(aes(ymin = lower_cc, ymax = upper_cc), width = 0.1) +
  facet_grid(name_segments ~ matrix) +
  ylab("Consensus Completeness (%)") +
  xlab("Read Subsampling") +
  mytheme
plot_flu_segments
