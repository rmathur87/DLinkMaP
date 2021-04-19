# libraries & dataset
import seaborn as sns
import matplotlib.pyplot as plt
# set a grey background (use sns.set_theme() if seaborn version 0.11.0 or above)
# sns.set(style="darkgrid")
df = sns.load_dataset("iris")

# creating a figure composed of two matplotlib.Axes objects (ax_box and ax_hist)
f, (ax_box, ax_hist) = plt.subplots(2, sharex=True, gridspec_kw={"height_ratios": (.1, .9)})

# assigning a graph to each ax
sns.boxplot(df["sepal_length"], ax=ax_box)
sns.histplot(data=df, x="sepal_length", ax=ax_hist)

print( df["sepal_length"] )

# Remove x axis name for the boxplot
ax_box.set(xlabel='')
plt.show()
