# T2D positive-control MR

This directory contains the runnable entry point (`run_mr.R`) and receives all
analysis outputs. The exposure file remains the sibling file
`../Mahajan.NatGenet2018b.T2D.European.txt` and is never overwritten.

Run from the project root:

```bash
Rscript T2D/Mahajan.NatGenet2018b.T2D.European_results/run_mr.R \
  --plink /absolute/path/to/plink
```

The analysis uses the Mahajan et al. 2018b European DIAGRAM T2D GWAS,
unadjusted for BMI, local EUR-reference LD clumping, the raw BCAC 2021 survival
results, and Slope-Hunter-adjusted survival results.
