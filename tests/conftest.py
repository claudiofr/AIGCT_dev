import pytest
import random
import context  # noqa: F401
from agct.container import VEBenchmarkContainer
from agct.repository import (
    VARIANT_PK_COLUMNS
)


@pytest.fixture
def ve_bm_container():
    return VEBenchmarkContainer()


@pytest.fixture
def ve_analyzer(ve_bm_container):
    return ve_bm_container.analyzer


@pytest.fixture
def ve_reporter(ve_bm_container):
    return ve_bm_container.reporter


@pytest.fixture
def ve_plotter(ve_bm_container):
    return ve_bm_container.plotter


def generate_random_floats(n, start, end):
    return [random.uniform(start, end) for _ in range(n)]


@pytest.fixture
def sample_user_scores1(ve_bm_container):
    user_scores_df = ve_bm_container._score_repo.get_all_by_task(
        "CANCER")
    user_variants = user_scores_df[VARIANT_PK_COLUMNS].drop_duplicates()
    random_idxs = random.sample(list(range(len(user_variants))), 2000)
    user_variants = user_variants.iloc[random_idxs]
    random_scores = generate_random_floats(len(user_variants), 0.01, 0.99)
    user_variants['RANK_SCORE'] = random_scores
    return user_variants


@pytest.fixture
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


@pytest.fixture
def sample_user_scores_col_name_map(ve_bm_container, sample_user_scores):
    user_scores_df = ve_bm_container._score_repo.get("CANCER",
                                                          "REVEL")
    user_scores_df["RANK_SCORE"] = user_scores_df["RANK_SCORE"].apply(
        lambda sc: sc + 0.1 if sc < 0.9 else 0.95)
    col_name_map = {"CHROMOSOME": "chr", "POSITION": "pos",
                    "REFERENCE_NUCLEOTIDE": "ref",
                    "ALTERNATE_NUCLEOTIDE": "alt",
                    "GENOME_ASSEMBLY": "assembly",
                    "RANK_SCORE": "score"}
    user_scores_df.rename(columns=col_name_map, inplace=True)
    col_name_map = {val: key for key, val in col_name_map.items()}
    return user_scores_df, col_name_map


@pytest.fixture
def ve_bm_query_mgr(ve_bm_container):
    return ve_bm_container.query_mgr


@pytest.fixture
def ve_analysis_exporter(ve_bm_container):
    return ve_bm_container.exporter


@pytest.fixture
def ve_data_validator(ve_bm_container):
    return ve_bm_container.data_validator
 
