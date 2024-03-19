A tutorial on Lipid Structure Enrichment Analysis (LSEA)
================

## load package and example data and perform clr transformation

``` r
library(lsea)

data(comp_data)
data(labels)

#transform composition data (features x samples) using centered log-ratio transformation
# note: this is not for normalized abundance data!!
clr_data <- lsea::clr.transform(comp_data)
```

## perform differential composition analysis on clr-transformed data

``` r
# to perform a paired t-test (samples must be ordered by pair)
de_df <- lsea::two.group.row.test(clr_data, labels, test = "t", paired = TRUE)
# to perform unpaired wilcox test
# note: wilcox test will throw warnings for every row with ties, so suppress those warnings
w_df <- suppressWarnings(lsea::two.group.row.test(clr_data, labels, test = "w", paired = FALSE))
```

### inspect the results

t-test results

|         |       stat |      mean1 |      mean2 |         dm |    pvalue |      padj |
|:--------|-----------:|-----------:|-----------:|-----------:|----------:|----------:|
| CE 12:0 | -0.2614509 | -0.1935301 |  0.1935301 | -0.3870601 | 0.7958048 | 0.9892482 |
| CE 14:0 |  1.3159110 |  0.9237276 | -0.9237276 |  1.8474552 | 0.1996890 | 0.9551987 |
| CE 14:1 |  0.4184487 |  0.2633849 | -0.2633849 |  0.5267698 | 0.6790554 | 0.9892482 |
| CE 15:0 | -0.9589883 | -0.5506678 |  0.5506678 | -1.1013357 | 0.3463982 | 0.9785673 |
| CE 16:0 | -0.9461001 | -0.7103156 |  0.7103156 | -1.4206311 | 0.3528080 | 0.9785673 |
| CE 16:1 |  0.7659404 |  0.4240928 | -0.4240928 |  0.8481855 | 0.4506099 | 0.9785673 |

Wilcoxon test results

|         |  stat |      mean1 |      mean2 |         dm |    pvalue |      padj |
|:--------|------:|-----------:|-----------:|-----------:|----------:|----------:|
| CE 12:0 | 385.5 | -0.1935301 |  0.1935301 | -0.3870601 | 0.7104152 | 0.9929020 |
| CE 14:0 | 416.0 |  0.9237276 | -0.9237276 |  1.8474552 | 0.3686215 | 0.9857487 |
| CE 14:1 | 345.0 |  0.2633849 | -0.2633849 |  0.5267698 | 0.7406416 | 0.9929020 |
| CE 15:0 | 345.0 | -0.5506678 |  0.5506678 | -1.1013357 | 0.7368086 | 0.9929020 |
| CE 16:0 | 246.0 | -0.7103156 |  0.7103156 | -1.4206311 | 0.0335855 | 0.7923619 |
| CE 16:1 | 369.0 |  0.4240928 | -0.4240928 |  0.8481855 | 0.9439534 | 0.9929020 |

## lipid structure enrichment analysis

``` r
# we can use the t-statistic to rank lipid species and perform enrichment analysis using the GSEA algorithm
# note: the returned Wilcoxon statistic is not directional
res <- lsea::lsea(de_df, rnk_name = "stat")
```

Top 5 positive results

|     | pathway      |      pval |      padj |        ES |      NES |
|:----|:-------------|----------:|----------:|----------:|---------:|
| 104 | TG_UFA_12-16 | 0.0061168 | 0.4342936 | 0.4911824 | 1.709831 |
| 191 | UFA_12-16    | 0.0061168 | 0.4342936 | 0.4911824 | 1.709831 |
| 281 | TG_UFA_16    | 0.0381238 | 0.8893568 | 0.5313508 | 1.559604 |
| 40  | CE_SFA_22-26 | 0.0099840 | 0.5670927 | 0.9674523 | 1.467142 |
| 230 | PE.O_MUFA_18 | 0.0571199 | 0.9594808 | 0.8365140 | 1.428722 |

Top 5 negative results

|     | pathway           |      pval |      padj |         ES |       NES |
|:----|:------------------|----------:|----------:|-----------:|----------:|
| 153 | FA_12-16          | 0.0354213 | 0.8893568 | -0.7230966 | -1.550059 |
| 61  | HexCER_MUFA_22-26 | 0.0246148 | 0.8893568 | -0.7448246 | -1.596637 |
| 157 | HexCER_22-26      | 0.0246148 | 0.8893568 | -0.7448246 | -1.596637 |
| 277 | TG_SFA_18         | 0.0259688 | 0.8893568 | -0.6808962 | -1.605585 |
| 193 | CE_PUFA_18        | 0.0051659 | 0.4342936 | -0.9356382 | -1.606246 |
| 144 | CE_17-20          | 0.0038641 | 0.4342936 | -0.6904010 | -1.849383 |

### lipid species structure annotation

``` r
# LSEA works by annotating the provided lipid species based on their LIPIDMAPS-style name
# and sorts them into multiple lipid sets based on class, tail length, and double bonds
lipid_anno <- lsea::annotate.lipid.species(rownames(de_df))
```

| Species | Class | Category | Total.Carbons | Longest.Tail | Total.DBs | Saturation | Chain |
|:--------|:------|:---------|--------------:|-------------:|----------:|:-----------|:------|
| CE 12:0 | CE    | Sterol   |            12 |           12 |         0 | SFA        | MCFA  |
| CE 14:0 | CE    | Sterol   |            14 |           14 |         0 | SFA        | LCFA  |
| CE 14:1 | CE    | Sterol   |            14 |           14 |         1 | MUFA       | LCFA  |
| CE 15:0 | CE    | Sterol   |            15 |           15 |         0 | SFA        | LCFA  |
| CE 16:0 | CE    | Sterol   |            16 |           16 |         0 | SFA        | LCFA  |
| CE 16:1 | CE    | Sterol   |            16 |           16 |         1 | MUFA       | LCFA  |
