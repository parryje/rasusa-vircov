pacman::p_load(tidyr, dplyr, stringr, ggplot2)

# do_something <- function(data_path, out_path, threads, myparam) {
#     # R code
# }
# 
# do_something(snakemake@input[[1]], snakemake@output[[1]], snakemake@threads, snakemake@config[["myparam"]])

#INPUT = snakemake@input[[1]]
INPUT = "../all.combined.tsv"

# clean up
dat = read.csv(INPUT, sep = "\t") %>% 
    filter(!is.na(consensus_completeness)) %>%
    separate(col = id, into = c("sample_name", "replicate"), sep = "-") %>%
    separate(col = sample_name, into = c("molecule", "matrix", "dilution", "subsample_limit_short"), sep = "_", remove = F) %>%
    mutate(replicate = gsub("replicate", "", replicate))

flu = dat %>% filter(name == "FLUAV") %>%
    group_by(sample_name, molecule, matrix, dilution, subsample_limit_short, name, replicate) %>%
    reframe(consensus_completeness = mean(consensus_completeness)) %>%
    ungroup()
# combine back flu into main df

dat = dat %>% filter(name != "FLUAV") %>%
    bind_rows(flu)

#compute std bars for error bars from replicates
stats = dat %>%
    group_by(molecule, matrix, dilution, subsample_limit_short, name) %>%
    mutate(
        std = sd(consensus_completeness),
        mean_cc = mean(consensus_completeness),
        lower = mean_cc - std,
        upper = mean_cc + std
        )

# now that we calculate stats, we can drop replicates
stats = stats %>%
    distinct(molecule, matrix, dilution, subsample_limit_short, name, .keep_all = T)
    
#plot
stats %>%
    ggplot(aes(x = subsample_limit_short, y = mean_cc)) +
    geom_line(aes(group = dilution, color = dilution)) +
    geom_point(aes(group = dilution, color = dilution)) + 
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
    facet_wrap(matrix~name)
    