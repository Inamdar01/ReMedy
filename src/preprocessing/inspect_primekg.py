import pandas as pd

path = "data/raw/primekg/primekg.csv"

print("Loading PrimeKG...")
df = pd.read_csv(path, low_memory=False)

print("\n========== BASIC INFORMATION ==========")
print("Rows:", len(df))
print("Columns:", len(df.columns))

print("\n========== COLUMN NAMES ==========")
for i, col in enumerate(df.columns):
    print(i, ":", col)

print("\n========== FIRST 5 ROWS ==========")
print(df.head())

print("\n========== DATA TYPES ==========")
print(df.dtypes)

print("\n========== MISSING VALUES ==========")
print(df.isna().sum().sort_values(ascending=False).head(20))

print("\n========== UNIQUE VALUES ==========")

for col in df.columns:
    if df[col].nunique() < 30:
        print(f"\n{col}:")
        print(df[col].unique())

print("\nInspection complete.")