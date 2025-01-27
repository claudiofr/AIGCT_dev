import context  # noqa: F401
from agct.analyzer import VEAnalyzer
from agct.exporter import VEAnalysisExporter

# from agct.model import VariantId  # noqa: F401


def test_export(
        variant_bm_analyzer: VEAnalyzer,
        sample_user_scores,
        ve_analysis_exporter: VEAnalysisExporter):
    metrics = variant_bm_analyzer.compute_metrics(
        "cancer", sample_user_scores, vep_min_overlap_percent=50,
        variant_vep_retention_percent=1, list_variants=True)
    ve_analysis_exporter.export_results(metrics, "./demo/output")
    pass

