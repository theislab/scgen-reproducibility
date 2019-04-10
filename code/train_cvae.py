import scgen
import scanpy as sc
import numpy as np

train = sc.read("../data/train_pbmc.h5ad")
valid = sc.read("../data/valid_pbmc.h5ad")
train = train[~((train.obs["cell_type"] == "CD4T") & (train.obs["condition"] == "stimulated"))]
z_dim = 20
network = scgen.CVAE(x_dimension=train.X.shape[1], z_dimension=z_dim, alpha=0.1, model_path="../models/CVAE/pbmc/all/models/scgen")
network.train(train, use_validation=True, valid_data=valid, n_epochs=100)
labels, _ = scgen.label_encoder(train)
train = sc.read("../data/train_pbmc.h5ad")
CD4T = train[train.obs["cell_type"] == "CD4T"]
unperturbed_data = train[((train.obs["cell_type"] == "CD4T") & (train.obs["condition"] == "control"))]
fake_labels = np.ones((len(unperturbed_data), 1))
predicted_cells = network.predict(unperturbed_data, fake_labels)
adata = sc.AnnData(predicted_cells, obs={"condition": ["pred"]*len(fake_labels)})
adata.var_names = CD4T.var_names
all_adata = CD4T.concatenate(adata)
all_adata.write("../data/reconstructed/CVAE_CD4T.h5ad")