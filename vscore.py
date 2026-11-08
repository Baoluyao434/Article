"""V-score computation (Treat vs. Control) for drug-perturbation
signature construction, compatible with drugreflector.

The v-score is (mean_treat - mean_control) /
               (sqrt(var_treat + var_control) + 1_{vars are zero}),
using population variance (ddof=0).
"""
import numpy as np
import pandas as pd
import anndata as ad


def compute_vscores(control_data, treat_data):
    """control_data / treat_data: shape (n_samples, n_genes)"""
    mean_control = np.mean(control_data, axis=0)
    mean_treat   = np.mean(treat_data,   axis=0)
    var_control  = np.var(control_data,  axis=0, ddof=0)
    var_treat    = np.var(treat_data,    axis=0, ddof=0)
    denom = np.sqrt(var_control + var_treat) + (var_control + var_treat == 0)
    return (mean_treat - mean_control) / denom


def run_vscore(df):
    """
    df : DataFrame (rows = genes, columns = samples), sample names
         containing "_Control" or "_PMOP".

    Returns
    -------
    adata       : AnnData with one observation "Treat_vs_Control"
    vscores_df  : DataFrame (gene, v_score, mean_control, mean_treat)
    """
    control_samples = [c for c in df.columns if "_Control" in c]
    treat_samples   = [c for c in df.columns if "_PMOP"    in c]
    if not control_samples or not treat_samples:
        raise ValueError("Column names must contain '_Control' or '_PMOP'")

    control_data = df[control_samples].values.T
    treat_data   = df[treat_samples].values.T

    vscores = compute_vscores(control_data, treat_data)
    genes   = df.index.tolist()

    adata = ad.AnnData(
        X=vscores.reshape(1, -1),
        obs=pd.DataFrame(
            index=["Treat_vs_Control"],
            data={
                "comparison":       "Treat_vs_Control",
                "reference":        "Control",
                "comparison_group": "Treat",
            },
        ),
        var=pd.DataFrame(index=genes, data={"gene_symbol": genes}),
    )

    vscores_df = pd.DataFrame({
        "gene":                     genes,
        "v_score_Treat_vs_Control": vscores,
        "mean_control":             np.mean(control_data, axis=0),
        "mean_treat":               np.mean(treat_data,   axis=0),
    }).sort_values("v_score_Treat_vs_Control", ascending=False)

    return adata, vscores_df