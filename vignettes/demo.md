A tutorial on Lipid Structure Enrichment Analysis (LSEA)
================

### load package and example data

``` r
library(lsea)

data(comp_data)
data(labels)
```

## CLR transformation

``` r
#transform composition data (features x samples) using centered log-ratio transformation
# note: this is not for normalized abundance data!!
clr_data <- lsea::clr.transform(comp_data)
```

## Differential composition analysis

``` r
# to perform a paired t-test (samples must be ordered by pair)
t_df <- lsea::two.group.row.test(clr_data, labels, test = "t", paired = TRUE)
# to perform unpaired wilcox test
# note: wilcox test will throw warnings for every row with ties, so suppress those warnings
w_df <- suppressWarnings(lsea::two.group.row.test(clr_data, labels, test = "w", paired = FALSE))
```

### inspect the results

``` r
# t-test results ordered by p-value
knitr::kable(head(t_df[order(t_df$pvalue),]))
```

|                |      stat |     mean1 |     mean2 |        dm |    pvalue |      padj |
|:---------------|----------:|----------:|----------:|----------:|----------:|----------:|
| TG 56:4-FA20:3 |  4.334250 |  2.168452 | -2.168452 |  4.336905 | 0.0001947 | 0.1738581 |
| TG 49:3-FA18:3 |  3.825307 |  2.280831 | -2.280831 |  4.561661 | 0.0007361 | 0.2324484 |
| TG 54:5-FA16:0 | -3.639901 | -1.953318 |  1.953318 | -3.906636 | 0.0011869 | 0.2324484 |
| FA 20:3        |  3.615850 |  1.678991 | -1.678991 |  3.357982 | 0.0012623 | 0.2324484 |
| TG 52:6-FA20:5 |  3.603903 |  1.856168 | -1.856168 |  3.712336 | 0.0013015 | 0.2324484 |
| PC 18:2_20:2   | -3.225378 | -1.399152 |  1.399152 | -2.798303 | 0.0033828 | 0.5034804 |

``` r
# Wilcoxon test results ordered by p-value
knitr::kable(head(w_df[order(w_df$pvalue),]))
```

|                |  stat |      mean1 |      mean2 |        dm |    pvalue |      padj |
|:---------------|------:|-----------:|-----------:|----------:|----------:|----------:|
| PC 16:1_18:1   | 624.0 |  1.6414554 | -1.6414554 |  3.282911 | 0.0000062 | 0.0055207 |
| FA 16:0        | 139.5 | -1.5086985 |  1.5086985 | -3.017397 | 0.0000831 | 0.0266629 |
| CE 18:1        | 140.0 | -1.6537323 |  1.6537323 | -3.307465 | 0.0000896 | 0.0266629 |
| PE P-16:0/16:1 | 581.0 |  1.3571488 | -1.3571488 |  2.714298 | 0.0001750 | 0.0390769 |
| TG 46:1-FA16:1 | 558.5 |  1.4270360 | -1.4270360 |  2.854072 | 0.0007628 | 0.1362407 |
| PE 18:0_18:1   | 554.5 |  0.6344112 | -0.6344112 |  1.268822 | 0.0009795 | 0.1457757 |

## lipid structure enrichment analysis

``` r
# we can use the t-statistic to rank lipid species and perform enrichment analysis using the GSEA algorithm
# note: the returned Wilcoxon statistic is not directional
res <- lsea::lsea(de_df, rnk_name = "stat")
```

Top 5 positive results

|     | pathway      |      pval |      padj |        ES |      NES |
|:----|:-------------|----------:|----------:|----------:|---------:|
| 104 | TG_UFA_12-16 | 0.0084202 | 0.5481303 | 0.4911824 | 1.707807 |
| 191 | UFA_12-16    | 0.0084202 | 0.5481303 | 0.4911824 | 1.707807 |
| 281 | TG_UFA_16    | 0.0367009 | 0.8951155 | 0.5313508 | 1.557362 |
| 40  | CE_SFA_22-26 | 0.0096502 | 0.5481303 | 0.9674523 | 1.479860 |
| 230 | PE.O_MUFA_18 | 0.0540984 | 0.9336334 | 0.8365140 | 1.442940 |

Top 5 negative results

|     | pathway           |      pval |      padj |         ES |       NES |
|:----|:------------------|----------:|----------:|-----------:|----------:|
| 153 | FA_12-16          | 0.0374257 | 0.8951155 | -0.7230966 | -1.541114 |
| 61  | HexCER_MUFA_22-26 | 0.0257426 | 0.8951155 | -0.7448246 | -1.587423 |
| 157 | HexCER_22-26      | 0.0257426 | 0.8951155 | -0.7448246 | -1.587423 |
| 277 | TG_SFA_18         | 0.0287253 | 0.8951155 | -0.6808962 | -1.599495 |
| 193 | CE_PUFA_18        | 0.0042918 | 0.5481303 | -0.9356382 | -1.607803 |
| 144 | CE_17-20          | 0.0033838 | 0.5481303 | -0.6904010 | -1.842860 |

### lipid species structure annotation

``` r
# LSEA works by annotating the provided lipid species based on their LIPIDMAPS-style name
# and sorts them into multiple lipid sets based on class, tail length, and double bonds
lipid_anno <- lsea::annotate.lipid.species(rownames(t_df))
```

| Species | Class | Category | Total.Carbons | Longest.Tail | Total.DBs | Saturation | Chain |
|:--------|:------|:---------|--------------:|-------------:|----------:|:-----------|:------|
| CE 12:0 | CE    | Sterol   |            12 |           12 |         0 | SFA        | MCFA  |
| CE 14:0 | CE    | Sterol   |            14 |           14 |         0 | SFA        | LCFA  |
| CE 14:1 | CE    | Sterol   |            14 |           14 |         1 | MUFA       | LCFA  |
| CE 15:0 | CE    | Sterol   |            15 |           15 |         0 | SFA        | LCFA  |
| CE 16:0 | CE    | Sterol   |            16 |           16 |         0 | SFA        | LCFA  |
| CE 16:1 | CE    | Sterol   |            16 |           16 |         1 | MUFA       | LCFA  |
