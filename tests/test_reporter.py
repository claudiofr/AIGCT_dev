import context  # noqa: F401
from agct.analyzer import VEAnalyzer
from agct.model import VEQueryCriteria
from agct.query import VEBenchmarkQueryMgr
from agct.reporter import VEAnalysisReporter
from agct.repository import (
    VARIANT_PK_COLUMNS
)

# from agct.model import VariantId  # noqa: F401


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

