
import dask.dataframe as dd
import networkx as nx
from community import community_louvain
import matplotlib.pyplot as plt

# Function to add edges to the graph incrementally
def add_edges_to_graph(df_chunk, graph):
    for _, row in df_chunk.iterrows():
        graph.add_edge(row['read1'], row['read2'], weight=row['distance'])
    return graph

# Load the data using Dask
df = dd.read_csv('closest_knn.txt', delim_whitespace=True, header=None, names=['read1', 'read2', 'distance'])

# Randomly sample 1/4 of the rows
sampled_df = df.sample(frac=0.5, random_state=42)

# Initialize an empty graph
G = nx.Graph()

# Process data in chunks to build the graph incrementally
num_chunks = len(sampled_df.to_delayed())
for idx, chunk in enumerate(sampled_df.to_delayed()):
    if (idx % 10) == 0:
        print("Adding chunk {0} of {1}...".format(idx, num_chunks))
    chunk = chunk.compute()
    G = add_edges_to_graph(chunk, G)

# Apply the Louvain method for community detection
partition = community_louvain.best_partition(G, weight='distance')

# Save clusters to a file
with open('clusters.txt', 'w') as f:
    for node, cluster in partition.items():
        f.write(f"{node}\t{cluster}\n")

# Optionally visualize a small sample of the graph
subgraph = G.subgraph(list(G.nodes)[:1000])  # Visualize only a subset if graph is too large
pos = nx.spring_layout(subgraph)
nx.draw(subgraph, pos, node_size=10, with_labels=False)
plt.savefig('graph_plot.png')
