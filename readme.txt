Codes for project of "GPHR: Geometry-Preserving Harmonization and Regression of Brain Connectivity via Reduced Log-Cholesky Representation"

## Folder Descriptions

├── src/               # Source code files
    ├── combat   # combat harmonization (Yu et al. 2018, Johnson et al. 2007)
    ├── comp     # competitive approaches (Desai et al. 2023, Zhang et al. 2023)
    ├── common   # code for proposed mapping
        ├── gh_map.m     # isometric map from SPSD with k0 rank to its log-Cholesky factor
        ├── gh_inv_map.m # inverse of the isometric map
        ├── vech.m       # vectorization of lower triangular matrices
    ├── spsd_test_q.m    # GPHR based hypothesis testing for Hub-centered Edges Detection with model selection for different number of latent factor
    ├── spsd_test.m      # GPHR based hypothesis testing for Hub-centered Edges Detection without model selection
├── RDA/                 # Code for Real Data Analysis
    ├── harmonization.m     # RLC-ComBat Harmonization applied in real data
    ├── RDA_ABIDE505_rcc_q.m  # Hub-centered Edges based on harmonized data