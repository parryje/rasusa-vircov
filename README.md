Current 4/07/2025 version:
    - Subsampling is working
    - vircov is working
    - vircov collate is working
    - plotting is working (`scripts/plot.R`) but:
        - tested on limited input samples and replicates.
        - Needs sprucing up (ugly atm)

Work to do:
    - Calculate the best subsampling read levels for input library. Some libraries have less reads, which affects their  upper limit. eg one sample might have 5M, while another might have 2M. Subsampling down to 5M on a 2m sample produces 2M reads.
    - Easiest way to do this is to use seqkit stats on the inputs. Find the lower limit,