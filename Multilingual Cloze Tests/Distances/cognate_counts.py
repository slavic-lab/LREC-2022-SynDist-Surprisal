import pandas as pd
from collections import defaultdict

cognate_data = pd.read_csv('northwindsun_matched_words_annotated.csv')

cognate_token_counts = defaultdict(lambda:0)
unique_cognates = defaultdict(lambda:[])

for i, row in cognate_data.iterrows():
    l1, l2 = row.L1, row.L2
    if type(row.Cognates) == str:
        if 'x' in row.Cognates:
            cognate_token_counts[(l1, l2)] += 1
            w1, w2 = row.Word1, row.Word2
            unique_cognates[(l1, l2)].append((w1.lower().strip(), w2.lower().strip()))
      
unique_cognate_counts = {}
for pair in unique_cognates:
    unique_cognate_counts[pair] = len(set(unique_cognates[pair]))