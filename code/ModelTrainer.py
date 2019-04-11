import subprocess
import sys


def main():
    if len(sys.argv) == 1:
        model_to_train = "all"
    else:
        model_to_train = sys.argv[1]
    if model_to_train == "all":
        command = "python ./vec_arth_pca.py"
        subprocess.call([command], shell=True)

        command = "python ./vec_arth.py"
        subprocess.call([command], shell=True)

        command = "python ./st_gan.py train"
        subprocess.call([command], shell=True)

        command = "python ./train_cvae.py"
        subprocess.call([command], shell=True)

        command = "python ./train_scGen.py"
        subprocess.call([command], shell=True)

    elif model_to_train == "PCA":
        command = "python ./vec_arth_pca.py"
        subprocess.call([command], shell=True)
    elif model_to_train == "VecArithm":
        command = "python ./vec_arth.py"
        subprocess.call([command], shell=True)
    elif model_to_train == "STGAN":
        command = "python ./st_gan.py train"
        subprocess.call([command], shell=True)
    elif model_to_train == "CVAE":
        command = "python ./train_cvae.py"
        subprocess.call([command], shell=True)
    elif model_to_train == "scGen":
        command = "python ./train_scGen.py"
        subprocess.call([command], shell=True)
