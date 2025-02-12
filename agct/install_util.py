# import context
import requests

import yaml
from random import random
import tarfile 
import importlib.resources as pkg_resources
import io
import os
import config
from .file_util import create_folder
from .container import VEBenchmarkContainer


def sample_user_scores(ve_bm_container):
    user_scores_df = ve_bm_container._score_repo.get(
        "CANCER", "REVEL")
    random_idxs = random.sample(list(range(len(user_scores_df))), 2000)
    user_variants = user_scores_df.iloc[random_idxs]
    user_variants['RANK_SCORE'] = user_variants['RANK_SCORE'].apply(
        lambda scor: scor + (random.uniform(0.05, 0.15) *
                             random.sample([1, -1], 1)[0])
        if scor < 0.84 and scor > 0.16 else scor)
    return user_variants


def get_sample_config() -> dict:
    with pkg_resources.open_text(config, 'agct.yaml.sample') as sample:
    # with open("./config/agct.yaml.sample") as sample:
        return yaml.safe_load(sample)


def init_config_file(config_dir: str = ".", data_dir: str = ".",
                     output_dir: str = "./analysis_output",
                     log_dir: str = "./log"):
    # create_folder(config_dir)
    config_file = os.path.join(config_dir, "agct.yaml")
    config = get_sample_config()
    config["repository"]["root_dir"] = data_dir
    config["output_dir"] = output_dir
    config["log"]["dir"] = log_dir
    with (open(config_file, "w") as
            conf_file):
        yaml.dump(config, conf_file)


def init_db(dir: str = "."):
    config = get_sample_config()
    url = config["repository"]["source_url"]
    response = requests.get(url, stream=True)
    tar_file = os.path.join(dir, "repo1.tar.gz")
    with open(tar_file, "wb") as file:
        for chunk in response.iter_content(chunk_size=1024):
            if chunk:
                file.write(chunk)
    # with tarfile.open(tar_file) as archive:
    #     archive.extractall(dir)


def check_install():
    init_config_file()
    create_folder("./demo/output")
    container = VEBenchmarkContainer()
    user_test_vep_scores = sample_user_scores(container)
    analyzer = container.analyzer
    metrics = analyzer.compute_metrics(
                "CANCER", user_test_vep_scores, vep_min_overlap_percent=50,
                variant_vep_retention_percent=1, list_variants=True)

    container.reporter.write_summary(metrics)
    container.plotter.plot_results(metrics)
    container.exporter.export_results(metrics, "./demo/output")








