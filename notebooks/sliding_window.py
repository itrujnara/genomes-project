#!/usr/bin/env python
# coding: utf-8

# In[27]:


from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from safetensors.torch import load_file
import seaborn as sns
from scipy.stats import entropy, iqr, ks_2samp, shapiro, wilcoxon
import polars as pl
import plotly.express as px

import torch


# ## Exploration

# In[28]:


sns.set_theme()


# In[29]:


def plot_file(path, **kwargs):
    tensors = load_file(path)
    probs = tensors["probs"]
    seq_entropy = entropy(probs.to(torch.float32).numpy(), axis=1, base=2)
    seq_ic = 2 - seq_entropy
    df = pl.DataFrame({"ic": seq_ic}).with_row_index()
    return sns.lineplot(df, x="index", y="ic", **kwargs)


# In[30]:


basepath = Path("/users/cn/itrujnara/genomes-project/infer_chromosome/one_chunk/tensors")


# In[31]:


#plot_file(basepath / "probs_context8192_stride16.safetensors")


# In[32]:


#plot_file(basepath / "probs_context8192_stride512.safetensors")


# In[33]:


plot_file(basepath / "evo2_7b_probs_context8192_stride1.safetensors")


# In[34]:


embeddings = load_file(basepath / "evo_embeddings_context8192_stride1024.safetensors")


# In[35]:


embeddings["embeddings"].shape


# ## Stride comparisons

# In[36]:


ics = {}
probs = {}

for stride in [1,16,32,64,128,256,512,1024,2048,4096]:
    tensors = load_file(basepath/f"evo2_7b_probs_context8192_stride{stride}.safetensors")
    seq_probs = tensors["probs"]
    probs[f"prob_{stride}"] = seq_probs
    seq_entropy = entropy(seq_probs.to(torch.float32).numpy(), axis=1, base=2)
    seq_ic = 2 - seq_entropy
    ics[f"ic_{stride}"] = seq_ic

df_ic = pl.DataFrame(ics).with_row_index()


# In[37]:


df_ic


# In[38]:


diff = (df_ic["ic_256"] - df_ic["ic_1"])
diff


# In[39]:


print(diff.mean())
print(diff.std())
print()
print(diff.max())
print(np.median(diff))


# In[40]:


df_diff = {"diff": diff}
#px.histogram(df_diff, x="diff")


# In[41]:


shapiro(diff)


# In[42]:


wilcoxon(diff)


# In[43]:


wilcoxon([1., 1., 1., 1., 1., 1., 1., 1., 1., 1.], [0.5, 0.5, 0.5, 0.5, 0.5, 1.5, 1.5, 1.5, 1.5, 1.5])


# In[44]:


df_ic_pivot = df_ic.unpivot(index=["index"])


# In[45]:


ax = sns.boxplot(df_ic_pivot, x="variable", y="value", color="#43C7CF")

ax.set_title("Information content distribution by stride", fontsize=15)
ax.set_xlabel("Stride", fontsize=14)
ax.set_ylabel("Information content", fontsize=14)
ax.set_xticks(list(ics.keys()), [1,16,32,64,128,256,512,1024,2048,4096])


# In[46]:


cols = ["ic_16", "ic_32", "ic_64", "ic_128", "ic_256", "ic_512", "ic_1024", "ic_2048", "ic_4096"]

tests = []

for col1 in cols:
    for col2 in cols:
        if col1 == col2:
            continue
        p = wilcoxon(df_ic[col1], df_ic[col2]).pvalue
        tests.append({"col1": col1, "col2": col2, "pvalue": p})

df_diff_tests = pl.DataFrame(tests).sort("col1", "col2")


# In[47]:


df_tests_pivot = (df_diff_tests
                  .pivot(index="col1", on="col2")
                  .sort(pl.col("col1").cast(pl.Enum(cols)))
                  .select(["col1"] + cols))

df_tests_matrix = df_tests_pivot.select(df_tests_pivot.columns[1:]).to_numpy()
df_tests_pivot


# In[48]:


df_tests_pivot.columns[1:]


# In[49]:


fig = px.imshow(df_tests_matrix, 
                x=df_tests_pivot.columns[1:], 
                y=df_tests_pivot["col1"], 
                text_auto=".3f",
                title="Wilcoxon p-value on IC vector",
                width=800)
fig.update_xaxes(side="top")
fig.show()


# In[50]:


def format_prob(p):
    if p < 1e-4:
        return "****"
    elif p < 1e-3:
        return "***"
    elif p < 1e-2:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "ns"


# In[51]:


def monte_carlo_ic(seed=42):
    alpha = [0.3,0.2,0.2,0.3]
    rng = np.random.default_rng(seed)
    p_arr = []
    q_arr = []
    np.random.seed(seed)
    for _ in range(131_072):
        p = rng.dirichlet(alpha)
        p_arr.append(p)
        q = rng.dirichlet(alpha)
        q_arr.append(q)
    p_arr = np.array(p_arr)
    q_arr = np.array(q_arr)
    p_ic = 2 - entropy(p_arr, axis=1, base=2)
    q_ic = 2 - entropy(q_arr, axis=1, base=2)
    return np.abs(p_ic - q_ic)


# In[52]:


monte_carlo_ic_vec = monte_carlo_ic()


# In[53]:


ic_diffs = {}

for col in cols:
    diff = np.abs(df_ic[col] - df_ic["ic_16"])
    ic_diffs[col] = diff

df_ic_diff = pl.DataFrame(ic_diffs)


# In[54]:


df_ic_diff


# In[55]:


px.box(df_ic_diff.unpivot(), x="variable", y="value")


# In[56]:


print("Mean difference from ic_16")
for col in cols:
    print(f"{col}: {df_ic_diff[col].mean():.4f}")


# In[57]:


print("Monte Carlo information content mean quantile")
for col in cols:
    col_mean = df_ic_diff[col].to_numpy().mean()
    p = (monte_carlo_ic_vec < col_mean).mean()
    print(f"{col}: {p:.4f} ({format_prob(p)})")


# In[58]:


print("Monte Carlo information content Kolmogorov p-value")
for col in cols:
    p = ks_2samp(df_ic_diff[col], monte_carlo_ic_vec, alternative="greater").pvalue
    print(f"{col}: {p:.4f} ({format_prob(p)})")


# ### Percentage difference

# In[ ]:





# ### Distance analysis

# In[59]:


dist_p16_p32 = ((probs["prob_32"] - probs["prob_16"])**2).sum(axis=1)
dist_p16_p32


# In[60]:


def euclidean_dist(arr1, arr2):
    return np.sqrt(((arr1 - arr2) ** 2).sum(axis=1))


# In[61]:


euclidean_dist(probs["prob_16"].numpy(), probs["prob_64"].numpy())


# In[62]:


df_probs = pl.DataFrame(probs)
df_probs


# In[63]:


cols_prob = ["prob_1", "prob_16", "prob_32", "prob_64", "prob_128", "prob_256", "prob_512", "prob_1024", "prob_2048", "prob_4096"]

dist_tests = []

for col1 in cols_prob:
    for col2 in cols_prob:
        if col1 == col2:
            continue
        arr1 = df_probs[col1].to_numpy()
        arr2 = df_probs[col2].to_numpy()
        dist = euclidean_dist(arr1, arr2)
        mean_dist = dist.mean()
        dist_tests.append({"col1": col1, "col2": col2, "dist": mean_dist})

df_dist_tests = pl.DataFrame(dist_tests).sort("col1", "col2")


# In[64]:


df_dist_tests


# In[65]:


df_dist_pivot = (df_dist_tests
                  .pivot(index="col1", on="col2")
                  .sort(pl.col("col1").cast(pl.Enum(cols_prob)))
                  .select(["col1"] + cols_prob))

df_dist_matrix = df_dist_pivot.select(df_dist_pivot.columns[1:])


# In[66]:


fig = px.imshow(df_dist_matrix, 
                x=df_dist_pivot.columns[1:], 
                y=df_dist_pivot["col1"], 
                text_auto=".3f", 
                title="Average Euclidean distance on probability vectors",
                width=800)
fig.update_xaxes(side="top")
fig.show()


# In[67]:


dist_vectors = {}

for col in cols_prob:
    dist_vectors[col] = euclidean_dist(probs["prob_1"], probs[col])

df_dist_long = pl.DataFrame(dist_vectors)

df_dist_long


# In[68]:


px.box(df_dist_long.unpivot(), x="variable", y="value", title="Euclidean distance from prob_11", labels={
    "variable": "Stride",
    "value": "Distance to stride 1"
}, width=800)


# In[69]:


upper_fence = np.quantile(df_dist_long["prob_32"], 0.75) + 1.5 * iqr(df_dist_long["prob_32"])
sum(df_dist_long["prob_32"] > upper_fence)


# In[70]:


def percent_outliers(arr, factor=1.5):
    fence = np.quantile(arr, 0.75) + factor * iqr(arr)
    return sum(arr > fence) / arr.shape[0]


# In[71]:


print("Factor 1.5 (standard)")
for col in cols_prob:
    print(f"{col}: {percent_outliers(df_dist_long[col]):.4f}")


# In[72]:


print("Factor 2.5")
for col in cols_prob:
    print(f"{col}: {percent_outliers(df_dist_long[col], 2.5):.4f}")


# In[73]:


print("Mean Euclidean distance from stride 1")
for col in cols_prob:
    col_mean = df_dist_long[col].mean()
    print(f"{col}: {col_mean:.4f}")


# In[74]:


def monte_carlo_euclidean(seed=42):
    alpha = [0.3,0.2,0.2,0.3]
    rng = np.random.default_rng(seed)
    p_arr = []
    q_arr = []
    np.random.seed(seed)
    for _ in range(131_072):
        p = rng.dirichlet(alpha)
        p_arr.append(p)
        q = rng.dirichlet(alpha)
        q_arr.append(q)
    p_arr = np.array(p_arr)
    q_arr = np.array(q_arr)
    return euclidean_dist(p_arr, q_arr)


# In[75]:


monte_carlo_vec_eucl = monte_carlo_euclidean()


# In[76]:


mean_256 = df_dist_long["prob_256"].mean()
mean_256


# In[77]:


(monte_carlo_vec_eucl < mean_256).mean()


# In[78]:


df_dist_long["prob_128"].mean()


# In[79]:


monte_carlo_vec_eucl.mean()


# In[80]:


cols_prob = ["prob_1", "prob_16", "prob_32", "prob_64", "prob_128", "prob_256", "prob_512", "prob_1024", "prob_2048", "prob_4096"]

print("Monte Carlo Euclidean mean quantiles")
for col in cols_prob:
    col_mean = df_dist_long[col].mean()
    p = (monte_carlo_vec_eucl < col_mean).mean()
    print(f"{col}: {p:.4f} ({format_prob(p)})")


# In[81]:


print("Monte Carlo Euclidean Kolmogorov p-values")
for col in cols_prob:
    p = ks_2samp(df_dist_long[col].to_numpy(), monte_carlo_vec_eucl, alternative="greater").pvalue
    print(f"{col}: {p} ({format_prob(p)})")


# In[82]:


def plot_cdfs(sample1, sample2, label1="Sample 1", label2="Sample 2"):
    from statsmodels.distributions.empirical_distribution import ECDF
    ecdf1 = ECDF(sample1)
    ecdf2 = ECDF(sample2)

    x = np.linspace(0, max(max(sample1), max(sample2)), 1000)
    plt.plot(x, ecdf1(x), label=label1)
    plt.plot(x, ecdf2(x), label=label2)
    plt.title("Empirical CDFs")
    plt.xlabel("Distance")
    plt.ylabel("CDF")
    plt.legend()
    plt.grid(True)
    plt.show()


# In[83]:


plot_cdfs(df_dist_long["prob_512"], monte_carlo_vec_eucl, "Model", "Monte Carlo")


# In[84]:


px.scatter(x=df_ic["ic_16"], y=df_dist_long["prob_256"],
           log_y=True, trendline="rolling", trendline_color_override="firebrick", trendline_options=dict(window=100),
           labels={'x': "Information content at stride 16", 'y': "Euclidean distance to stride 16 (log10)"},
           title="Euclidean distance between stride 16 and 256 vs. information content",
           width=800
          )


# ### KL divergence

# In[85]:


div_1_512 = entropy(probs["prob_1"], probs["prob_512"], axis=1, base=2)


# In[86]:


div_1_512.shape


# In[87]:


div_1_512.mean()


# In[88]:


div_1_4096 = entropy(probs["prob_1"], probs["prob_4096"], axis=1, base=2)


# In[89]:


div_1_4096.mean()


# In[90]:


def monte_carlo_divergence(seed=42):
    alpha = [0.3,0.2,0.2,0.3]
    rng = np.random.default_rng(seed)
    p_arr = []
    q_arr = []
    np.random.seed(seed)
    for _ in range(131_072):
        p = rng.dirichlet(alpha)
        p_arr.append(p)
        q = rng.dirichlet(alpha)
        q_arr.append(q)
    p_arr = np.array(p_arr)
    q_arr = np.array(q_arr)
    return entropy(p_arr, q_arr, axis=1, base=2)


# In[91]:


monte_carlo_vec_div = monte_carlo_divergence()


# In[92]:


monte_carlo_vec_div


# In[93]:


monte_carlo_vec_div.shape


# In[94]:


monte_carlo_vec_div.mean()


# In[95]:


torch.random.manual_seed(42)

idx = torch.randperm(122880)

perm_probs = probs["prob_1"][idx, :]

perm_vec_div = entropy(probs["prob_1"], perm_probs, axis=1, base=2)

perm_vec_div


# In[96]:


perm_vec_div.mean()


# In[97]:


div_vectors = {}

for col in cols_prob:
    div_vectors[col] = entropy(probs["prob_1"], probs[col], axis=1, base=2)

df_div_long = pl.DataFrame(div_vectors)

df_div_long


# In[98]:


div_vectors = {}

for col in cols_prob:
    div_vectors[col] = entropy(probs[col], probs["prob_1"], axis=1, base=2)

df_div_rev = pl.DataFrame(div_vectors)

df_div_rev


# In[99]:


print("Mean KL divergence from stride 1")
for col in cols_prob:
    col_mean = df_div_long[col].mean()
    print(f"{col}: {col_mean:.4f}")


# In[100]:


print("Mean KL divergence from stride 1 (reverse)")
for col in cols_prob:
    col_mean = df_div_rev[col].mean()
    print(f"{col}: {col_mean:.4f}")


# In[101]:


print("Monte Carlo KL mean quantiles")
for col in cols_prob:
    col_mean = df_div_long[col].mean()
    p = (monte_carlo_vec_div < col_mean).mean()
    print(f"{col}: {p:.6f} ({format_prob(p)})")


# In[102]:


print("Monte Carlo KL mean quantiles (reverse)")
for col in cols_prob:
    col_mean = df_div_rev[col].mean()
    p = (monte_carlo_vec_div < col_mean).mean()
    print(f"{col}: {p:.6f} ({format_prob(p)})")


# In[103]:


print("Permutation KL mean quantiles")
for col in cols_prob:
    col_mean = df_div_long[col].mean()
    p = (perm_vec_div < col_mean).mean()
    print(f"{col}: {p:.6f} ({format_prob(p)})")


# In[104]:


print("Permutation KL mean quantiles (reverse)")
for col in cols_prob:
    col_mean = df_div_rev[col].mean()
    p = (perm_vec_div < col_mean).mean()
    print(f"{col}: {p:.6f} ({format_prob(p)})")


# In[105]:


print("Monte Carlo KL Kolmogorov p-values")
for col in cols_prob:
    p = ks_2samp(df_div_long[col].to_numpy(), monte_carlo_vec_div, alternative="greater").pvalue
    print(f"{col}: {p} ({format_prob(p)})")


# In[106]:


print("Monte Carlo KL Kolmogorov p-values (reverse)")
for col in cols_prob:
    p = ks_2samp(df_div_rev[col].to_numpy(), monte_carlo_vec_div, alternative="greater").pvalue
    print(f"{col}: {p} ({format_prob(p)})")


# In[107]:


print("Permutation KL Kolmogorov p-values")
for col in cols_prob:
    p = ks_2samp(df_div_long[col].to_numpy(), perm_vec_div, alternative="greater").pvalue
    print(f"{col}: {p} ({format_prob(p)})")


# In[108]:


print("Permutation KL Kolmogorov p-values (reverse)")
for col in cols_prob:
    p = ks_2samp(df_div_rev[col].to_numpy(), perm_vec_div, alternative="greater").pvalue
    print(f"{col}: {p} ({format_prob(p)})")


# In[113]:


ax = sns.boxplot(df_div_long.unpivot().filter(pl.col("variable") != "prob_1"), x="variable", y="value")

ax.set_title("Per-position KL divergence to stride 1", fontsize=15)
ax.set_xlabel("Stride", fontsize=14)
ax.set_ylabel("KL divergence to stride 1", fontsize=14)
ax.set_xticks(list(probs.keys())[1:], [16,32,64,128,256,512,1024,2048,4096])


# ### Cosine distance

# In[ ]:


def cos_similarity(a, b):
    normA = np.linalg.norm(a, axis=1)
    normB = np.linalg.norm(b, axis=1)
    prod = a * b
    cos_sim = prod.sum(axis=1)
    cos_sim = cos_sim / normA
    cos_sim = cos_sim / normB
    return cos_sim


# In[ ]:


cos_16_32 = cos_similarity(probs["prob_16"], probs["prob_32"])


# In[ ]:


#px.histogram(cos_16_32)


# In[ ]:


cols_prob = ["prob_16", "prob_32", "prob_64", "prob_128", "prob_256", "prob_512", "prob_1024", "prob_2048", "prob_4096"]

cos_tests = []

for col1 in cols_prob:
    for col2 in cols_prob:
        if col1 == col2:
            continue
        arr1 = df_probs[col1].to_numpy()
        arr2 = df_probs[col2].to_numpy()
        cos_sim = cos_similarity(arr1, arr2).mean()
        cos_tests.append({"col1": col1, "col2": col2, "cos_sim": cos_sim})

df_cos_tests = pl.DataFrame(cos_tests).sort("col1", "col2")


# In[ ]:


df_cos_pivot = (df_cos_tests
                  .pivot(index="col1", on="col2")
                  .sort(pl.col("col1").cast(pl.Enum(cols_prob)))
                  .select(["col1"] + cols_prob))

df_cos_matrix = df_cos_pivot.select(df_cos_pivot.columns[1:])


# In[ ]:


fig = px.imshow(df_cos_matrix, 
                x=df_cos_pivot.columns[1:], 
                y=df_cos_pivot["col1"], 
                text_auto=".3f", 
                title="Average cosine similarity on probability vectors",
                width=800)
fig.update_xaxes(side="top")
fig.show()


# In[ ]:


cos_vectors = {}

for col in cols_prob:
    cos_vectors[col] = cos_similarity(probs["prob_16"].numpy(), probs[col].numpy())

df_cos_long = pl.DataFrame(cos_vectors)


# In[ ]:


print("Mean cosine similarity")
for col in cols_prob:
    col_mean = df_cos_long[col].mean()
    print(f"{col}: {col_mean:.4f}")


# In[ ]:


def monte_carlo_cosine(seed=42):
    alpha = [0.3,0.2,0.2,0.3]
    rng = np.random.default_rng(seed)
    p_arr = []
    q_arr = []
    np.random.seed(seed)
    for _ in range(1_000_000):
        p = rng.dirichlet(alpha)
        p_arr.append(p)
        q = rng.dirichlet(alpha)
        q_arr.append(q)
    p_arr = np.array(p_arr)
    q_arr = np.array(q_arr)
    return cos_similarity(p_arr, q_arr)


# In[ ]:


monte_carlo_vec_cosine = monte_carlo_cosine()


# In[ ]:


cols_prob = ["prob_32", "prob_64", "prob_128", "prob_256", "prob_512", "prob_1024", "prob_2048", "prob_4096"]

print("Monte Carlo cosine mean quantiles")
for col in cols_prob:
    col_mean = df_cos_long[col].mean()
    p = (monte_carlo_vec_cosine > col_mean).mean()
    print(f"{col}: {p:.4f} ({format_prob(p)})")


# In[ ]:


cols_prob = ["prob_32", "prob_64", "prob_128", "prob_256", "prob_512", "prob_1024", "prob_2048", "prob_4096"]

print("Monte Carlo cosine Kolmogorov p-values")
for col in cols_prob:
    p = ks_2samp(df_cos_long[col].to_numpy(), monte_carlo_vec_cosine, alternative="less").pvalue
    print(f"{col}: {p} ({format_prob(p)})")


# In[ ]:


plot_cdfs(df_cos_long["prob_256"], monte_carlo_vec_cosine, "Model", "Monte Carlo")


# ## Context size

# In[ ]:


context_ics = {}
context_probs = {}

for context in [512,1024,2048,4096,8192,16384]:
    tensors = load_file(basepath/f"evo2_7b_probs_context{context}_stride128.safetensors")
    seq_probs = tensors["probs"][-114688:,:]
    context_probs[f"prob_{context}"] = seq_probs
    seq_entropy = entropy(seq_probs.to(torch.float32).numpy(), axis=1, base=2)
    seq_ic = 2 - seq_entropy
    context_ics[f"ic_{context}"] = seq_ic

df_context = pl.DataFrame(context_ics).with_row_index()


# In[ ]:


df_context


# In[ ]:


df_context_pivot = df_context.unpivot(index=["index"])


# In[ ]:


px.box(df_context_pivot, x="variable", y="value", title="Information content distribution by context size", width=800)


# In[ ]:


cols_ic = ["ic_512", "ic_1024", "ic_2048", "ic_4096", "ic_8192", "ic_16384"]
context_diffs = {}

for col in cols_ic:
    diff = np.abs(df_context[col] - df_context["ic_16384"])
    context_diffs[col] = diff

df_context_diff = pl.DataFrame(context_diffs)


# In[ ]:


df_context_diff


# In[ ]:


ax = sns.boxplot(df_context_diff.unpivot(), x="variable", y="value")
ax.set_title("IC difference from context 16384")


# In[ ]:


print("Mean IC difference")
for col in cols_ic:
    col_mean = df_context_diff[col].mean()
    print(f"{col}: {col_mean:.4f}")


# In[ ]:


print("Monte Carlo mean IC difference quantile")
for col in cols_ic:
    col_mean = df_context_diff[col].mean()
    p = (monte_carlo_ic_vec < col_mean).mean()
    print(f"{col}: {p:.4f} ({format_prob(p)})")


# In[ ]:


print("Monte Carlo IC difference Kolmogorov p-value")
for col in cols_ic:
    p = ks_2samp(df_context_diff[col], monte_carlo_ic_vec, alternative="greater").pvalue
    print(f"{col}: {p:.4f} ({format_prob(p)})")


# In[ ]:


cols_prob = ["prob_512", "prob_1024", "prob_2048", "prob_4096", "prob_8192", "prob_16384"]

dist_vectors_context = {}

for col in cols_prob:
    dist_vectors_context[col] = euclidean_dist(context_probs["prob_16384"], context_probs[col])

df_context_long = pl.DataFrame(dist_vectors_context)

df_context_long


# In[ ]:


px.box(df_context_long.unpivot(), x="variable", y="value", title="Euclidean distance from prob_8192", labels={
    "variable": "Stride",
    "value": "Distance to context 8192"
}, width=800)


# In[ ]:


"""
px.scatter(x=df_context["ic_8192"], y=df_context_long["prob_4096"],
           log_y=True, trendline="rolling", trendline_color_override="firebrick", trendline_options=dict(window=100),
           labels={'x': "Information content at context 8192", 'y': "Euclidean distance to context 8192 (log10)"},
           title="Euclidean distance between context 4096 and 8192 vs. information content"
          )
"""


# In[ ]:


mean_4096 = df_context_long["prob_4096"].mean()
mean_4096


# In[ ]:


(monte_carlo_vec_eucl < mean_4096).mean()


# In[ ]:


print("Mean Euclidean distance")
for col in cols_prob:
    col_mean = df_context_long[col].mean()
    print(f"{col}: {col_mean:.4f}")


# In[ ]:


print("Monte Carlo Euclidean distance mean quantile")
for col in cols_prob:
    col_mean = df_context_long[col].mean()
    p = (monte_carlo_vec_eucl < col_mean).mean()
    print(f"{col}: {p:.4f} ({format_prob(p)})")


# In[ ]:


print("Monte Carlo Euclidean distance Kolmogorov p-value")
for col in cols_prob:
    p = ks_2samp(df_context_long[col], monte_carlo_vec_eucl, alternative="greater").pvalue
    print(f"{col}: {p:.4f} ({format_prob(p)})")

