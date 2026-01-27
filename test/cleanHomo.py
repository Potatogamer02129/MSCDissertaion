import pandas as pd

df = pd.read_csv("insect_ligand_targets_raw.csv")

# Keep rows that look insect-related
df_insect = df[
    df["Species"].str.contains(
        "Drosophila|Spodoptera|Aedes|Anopheles|Bombyx|Helicoverpa|Plutella|Nilaparvata|insect",
        case=False,
        na=False
    )
]

print("Before:", len(df))
print("After insect filter:", len(df_insect))

df_insect.to_csv(
    "insect_ligand_targets_insect_only.csv",
    index=False
)
