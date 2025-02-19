import context  # noqa: F401
from aigct.analyzer import VEAnalyzer
from aigct.model import VEQueryCriteria
from aigct.query import VEBenchmarkQueryMgr
from aigct.reporter import VEAnalysisReporter
from aigct.repository import (
    VARIANT_PK_COLUMNS
)

# from aigct.model import VariantId  # noqa: F401


def test_write_summary_stdout(
        ve_analyzer: VEAnalyzer,
        sample_user_scores,
        ve_reporter: VEAnalysisReporter):
    metrics = ve_analyzer.compute_metrics(
        "CANCER", sample_user_scores, "UserVep",
        vep_min_overlap_percent=50,
        variant_vep_retention_percent=1, list_variants=True)
    ve_reporter.write_summary(metrics)
    pass


def test_write_summary_file(
        ve_analyzer,
        sample_user_scores,
        ve_reporter: VEAnalysisReporter):
    metrics = ve_analyzer.compute_metrics(
        "CANCER", sample_user_scores, "UserVep",
        list_variants=True)
    ve_reporter.write_summary(metrics, "./demo/output")
    pass

