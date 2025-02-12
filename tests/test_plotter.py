import context  # noqa: F401
from agct.analyzer import VEAnalyzer
from agct.model import VEQueryCriteria
from agct.query import VEBenchmarkQueryMgr
from agct.plotter import VEAnalysisPlotter
from agct.repository import (
    VARIANT_PK_COLUMNS
)

# from agct.model import VariantId  # noqa: F401


def test_plot_results(
        ve_analyzer: VEAnalyzer,
        sample_user_scores,
        ve_plotter: VEAnalysisPlotter):
    metrics = ve_analyzer.compute_metrics(
        "CANCER", sample_user_scores, "UserVep",
        vep_min_overlap_percent=50,
        variant_vep_retention_percent=1, list_variants=True)
    ve_plotter.plot_results(metrics)
    pass


def test_plot_results_user_vep_only(
        ve_analyzer: VEAnalyzer,
        sample_user_scores,
        ve_plotter: VEAnalysisPlotter):
    metrics = ve_analyzer.compute_metrics(
        "CANCER", sample_user_scores, "UserVep",
        include_variant_effect_sources=False,
        list_variants=True)
    ve_plotter.plot_results(metrics)
    pass


def test_plot_results_system_veps_only(
        ve_analyzer: VEAnalyzer,
        sample_user_scores,
        ve_plotter: VEAnalysisPlotter):
    metrics = ve_analyzer.compute_metrics(
        "CANCER",
        vep_min_overlap_percent=50,
        variant_vep_retention_percent=1, list_variants=True)
    ve_plotter.plot_results(metrics)
    pass


def test_plot_results_file(
        ve_analyzer: VEAnalyzer,
        sample_user_scores,
        ve_plotter: VEAnalysisPlotter):
    metrics = ve_analyzer.compute_metrics(
        "CANCER", sample_user_scores, "UserVep",
        vep_min_overlap_percent=50,
        variant_vep_retention_percent=1, list_variants=True)
    ve_plotter.plot_results(metrics, dir="./demo/output")
    pass


