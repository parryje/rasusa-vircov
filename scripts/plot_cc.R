pacman::p_load(tidyr, dplyr, stringr, ggplot2, forcats)

# do_something <- function(data_path, out_path, threads, myparam) {
#     # R code
# }
# 
# do_something(snakemake@input[[1]], snakemake@output[[1]], snakemake@threads, snakemake@config[["myparam"]])

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
    mutate(sample_name = str_replace(sample_name, pattern = "Control_", replacement = "Control")) %>%
    mutate(sample_name = str_replace(sample_name, pattern = "_cDNA", replacement = ""))

# split sample_name into needed fields
dat3 = dat2 %>%
    separate(col = sample_name, into = c("molecule", "matrix", "dilution", "subsample_limit_short"), sep = "_", remove = F, extra = "merge") %>%
    mutate(replicate = gsub("replicate", "", replicate))

flu = dat3 %>% filter(name == "FLUAV") %>%
    group_by(sample_name, molecule, matrix, dilution, subsample_limit_short, name, replicate) %>%
    reframe(consensus_completeness = mean(consensus_completeness)) %>%
    ungroup()

# combine back flu into main df
dat_final = dat3 %>% filter(name != "FLUAV") %>%
    bind_rows(flu)

#compute std bars for error bars from replicates
stats = dat_final %>%
    group_by(molecule, matrix, dilution, subsample_limit_short, name) %>%
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

#plot
plot_cc = stats %>%
    ggplot(aes(x = subsample_limit_short, y = mean_cc)) +
    geom_line(aes(group = dilution, color = dilution)) +
    geom_point(aes(group = dilution, color = dilution)) + 
    geom_errorbar(aes(ymin = lower_cc, ymax = upper_cc), width = 0.1) +
    facet_grid(name ~ matrix) + 
    ylab("Consensus Completeness (%)") + 
    xlab("Read Subsampling") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(filename = "subsamp_cc.png", plot = plot_cc, device = "png", width = 297, height = 210, units = "mm")

plot_rc = stats %>%
    ggplot(aes(x = subsample_limit_short, y = mean_rc)) +
    geom_line(aes(group = dilution, color = dilution)) +
    geom_point(aes(group = dilution, color = dilution)) + 
    geom_errorbar(aes(ymin = lower_rc, ymax = upper_rc), width = 0.1) +
    geom_hline(yintercept=80, linetype='dotted', col = 'red') +
    facet_grid(name ~ matrix) + 
    ylab("Consensus Completeness (%)") + 
    xlab("Read Subsampling") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(filename = "subsamp_rc.png", plot = plot_rc, device = "png", width = 297, height = 210, units = "mm")
