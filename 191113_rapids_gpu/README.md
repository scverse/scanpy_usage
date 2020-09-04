# Speeding up Scanpy using GPU-enabled RAPIDS

We will use a machine on Google Cloud Platform, since it makes it easy
to get a GPU-enabled machine with all the software pre-configured.
(See https://cloud.google.com/deep-learning-vm/docs/images.)

You may wish to change the zone specified here:

```bash
export IMAGE_FAMILY="rapids-latest-gpu-experimental"
export ZONE="us-central1-a"
export INSTANCE_NAME="$USER-scanpy"

gcloud compute instances create $INSTANCE_NAME \
  --zone=$ZONE \
  --image-family=$IMAGE_FAMILY \
  --image-project=deeplearning-platform-release \
  --machine-type=n1-standard-16 \
  --maintenance-policy=TERMINATE \
  --accelerator="type=nvidia-tesla-t4,count=1" \
  --metadata="install-nvidia-driver=True"
```

Connect to the instance with the following command. Note that it may
take a few minutes before it is ready to accept connections.

```bash
gcloud compute ssh --zone $ZONE $INSTANCE_NAME
```

Once connected, type `nvidia-smi` to check that there is indeed a GPU.

```
$ nvidia-smi
+-----------------------------------------------------------------------------+
| NVIDIA-SMI 410.104      Driver Version: 410.104      CUDA Version: 10.0     |
|-------------------------------+----------------------+----------------------+
| GPU  Name        Persistence-M| Bus-Id        Disp.A | Volatile Uncorr. ECC |
| Fan  Temp  Perf  Pwr:Usage/Cap|         Memory-Usage | GPU-Util  Compute M. |
|===============================+======================+======================|
|   0  Tesla T4            Off  | 00000000:00:04.0 Off |                    0 |
| N/A   65C    P0    29W /  70W |      0MiB / 15079MiB |      0%      Default |
+-------------------------------+----------------------+----------------------+

+-----------------------------------------------------------------------------+
| Processes:                                                       GPU Memory |
|  GPU       PID   Type   Process name                             Usage      |
|=============================================================================|
|  No running processes found                                                 |
+-----------------------------------------------------------------------------+
```

Install scanpy (with RAPIDS GPU support):

```bash
sudo /opt/anaconda3/bin/pip uninstall -y enum34 # https://github.com/iterative/dvc/issues/1995
pip install --user git+https://github.com/theislab/scanpy
```

Install louvain (to run on a CPU, for comparison):

```bash
sudo apt-get install -y libz-dev libxml2-dev && pip install --user louvain # takes a while to install from source
```

Checkout this repo in the instance to give access to the scripts.

```bash
git clone https://github.com/theislab/scanpy_usage
cd scanpy_usage/191113_rapids_gpu
```

Copy test data to the instance. Note that you will have to change the user
(`-u`) to be your GCP project since the data is in a requester pays bucket.

```bash
gsutil -u hca-scale cp gs://ll-sc-data/10x/1M_neurons_filtered_gene_bc_matrices_h5.h5 1M_neurons_filtered_gene_bc_matrices_h5.h5
```

Run the analyses:

```bash
# Regular (CPU only)
python cluster_130K.py 1M_neurons_filtered_gene_bc_matrices_h5.h5
# GPU: neighbors only
python cluster_130K_gpu_neighbors.py 1M_neurons_filtered_gene_bc_matrices_h5.h5
# GPU: louvain only
python cluster_130K_gpu_louvain.py 1M_neurons_filtered_gene_bc_matrices_h5.h5
# GPU: umap only
python cluster_130K_gpu_umap.py 1M_neurons_filtered_gene_bc_matrices_h5.h5
# GPU: neighbors, louvain, and umap
python cluster_130K_gpu.py 1M_neurons_filtered_gene_bc_matrices_h5.h5
```

## Timings

| Step      | CPU time (s) | GPU time (s) | Speedup |
| --------- | ------------ | ------------ | ------- |
| Neighbors | 47           | 15           | 3x      |
| Louvain   | 70           | 1            | 70x     |
| UMAP      | 186          | 15           | 12x     |

Note that Scanpy only reports times to the nearest second, so the speedup may be quite inaccurate.

## Figures

The following is the UMAP clustering produced by running regular Scanpy:

![Regular (CPU only)](figures/umap_130K.png)

The following is equivalent, but using RAPIDS for nearest neighbor calculations.
Notice that there is one more cluster (community), and the layout is different - probably
because RAPIDS uses exact nearest neighbors, but UMAP uses approximate nearest neighbors.

![GPU: neighbors only](figures/umap_130K_gpu_neighbors.png)

This uses RAPIDS for Louvain (and regular Scanpy defaults for the rest).
Notice that the number of clusters (communities) increases from 27 to 44.

![GPU: louvain only](figures/umap_130K_gpu_louvain.png)

This uses RAPIDS for UMAP (and regular Scanpy defaults for the rest).
Notice that the figure is bunched in the corner, possibly because of outlying points?

![GPU: umap only](figures/umap_130K_gpu_umap.png)

To view the figures, run the following from your local machine to copy
them locally:

```bash
gcloud compute scp --zone $ZONE --recurse $INSTANCE_NAME:scanpy_usage/191113_rapids_gpu/figures figures_local
```

## Finishing up

```bash
gcloud compute instances delete --zone $ZONE --quiet $INSTANCE_NAME
```
