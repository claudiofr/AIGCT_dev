import context  # noqa: F401
from agct.analyzer import VEAnalyzer
from agct.model import VEQueryCriteria
from agct.query import VEBenchmarkQueryMgr
from agct.reporter import VEAnalysisReporter
from agct.pd_util import filter_dataframe_by_list

# from agct.model import VariantId  # noqa: F401


def test_compute_metrics_basic(variant_bm_analyzer,
                               sample_user_scores):
    metrics = variant_bm_analyzer.compute_metrics(
        "cancer", sample_user_scores, list_variants=True)


def test_compute_metrics_col_name_map_include_ve_sources(
        variant_bm_analyzer: VEAnalyzer,
        sample_user_scores_col_name_map):
    user_scores, col_name_map = sample_user_scores_col_name_map
    metrics = variant_bm_analyzer.compute_metrics(
        "cancer", user_scores, col_name_map, 
        ['REVEL', 'EVE'], list_variants=True)
    metrics


def test_compute_metrics_exclude_ve_sources(
        variant_bm_analyzer: VEAnalyzer,
        sample_user_scores_col_name_map):
    user_scores, col_name_map = sample_user_scores_col_name_map
    metrics = variant_bm_analyzer.compute_metrics(
        "cancer", user_scores, col_name_map, 
        ['REVEL', 'EVE'], False, list_variants=True)
    metrics


def test_compute_metrics_query_params_genes(
        variant_bm_analyzer: VEAnalyzer,
        variant_query_mgr: VEBenchmarkQueryMgr,
        sample_user_scores):
    qry = VEQueryCriteria(['MTOR', 'PTEN'])
    metrics = variant_bm_analyzer.compute_metrics(
        "cancer", sample_user_scores, None,
        ['REVEL', 'EVE'], False, qry, list_variants=True)
    qry = VEQueryCriteria(variant_ids=metrics.variants_included)
    variants = variant_query_mgr.get_variants_by_task("cancer", qry)
    genes = variants['GENE_SYMBOL'].unique()
    metrics


def test_compute_metrics_query_params_exclude_genes(
        variant_bm_analyzer: VEAnalyzer,
        variant_query_mgr: VEBenchmarkQueryMgr,
        sample_user_scores):
    qry = VEQueryCriteria(['MTOR', 'PTEN'], include_genes=False)
    metrics = variant_bm_analyzer.compute_metrics(
        "cancer", sample_user_scores, None, None,
        False, qry, list_variants=True)
    assert metrics.num_variants_included > 0
    qry = VEQueryCriteria(variant_ids=metrics.variants_included)
    variants = variant_query_mgr.get_variants_by_task("cancer", qry)
    genes = variants['GENE_SYMBOL'].unique()
    assert all([gene not in ['MTOR', 'PTEN'] for gene in genes])


def test_compute_metrics_query_params_filter(
        variant_bm_analyzer: VEAnalyzer,
        variant_query_mgr: VEBenchmarkQueryMgr,
        sample_user_scores):
    qry = VEQueryCriteria(filter_name="Oncogene")
    metrics = variant_bm_analyzer.compute_metrics(
        "cancer", sample_user_scores, None, None,
        False, qry, list_variants=True)
    qry = VEQueryCriteria(variant_ids=metrics.variants_included)
    variants = variant_query_mgr.get_variants_by_task("cancer", qry)
    genes = variants['GENE_SYMBOL'].unique()
    filter = variant_query_mgr.get_variant_filter("cancer",
                                                  "Oncogene")
    filter_genes = filter.filter_genes['GENE_SYMBOL']
    intersect = set(genes) & set(filter_genes)
    assert len(intersect) == len(genes)
    metrics


def test_compute_metrics_query_params_variant_ids(
        variant_bm_analyzer: VEAnalyzer,
        variant_query_mgr: VEBenchmarkQueryMgr,
        sample_user_scores):
    qry = VEQueryCriteria(variant_ids=sample_user_scores[:1500])
    metrics = variant_bm_analyzer.compute_metrics(
        "cancer", sample_user_scores, None, None,
        False, qry, list_variants=True)
    assert metrics.num_variants_included == 1500


def test_compute_metrics_query_params_exclude_variant_ids(
        variant_bm_analyzer: VEAnalyzer,
        variant_query_mgr: VEBenchmarkQueryMgr,
        sample_user_scores):
    qry = VEQueryCriteria(variant_ids=sample_user_scores[:1500],
                          include_variant_ids=False)
    metrics = variant_bm_analyzer.compute_metrics(
        "cancer", sample_user_scores, None, None,
        False, qry, list_variants=True)
    assert metrics.num_variants_included == 500


def test_compute_metrics_query_params_allele_freq(
        variant_bm_analyzer: VEAnalyzer,
        variant_query_mgr: VEBenchmarkQueryMgr,
        sample_user_scores,
        variant_bm_reporter: VEAnalysisReporter):
    qry = VEQueryCriteria(allele_frequency_operator=">",
                          allele_frequency=1.0e-7)
    metrics = variant_bm_analyzer.compute_metrics(
        "cancer", sample_user_scores, None, None,
        False, qry, list_variants=True)
    qry1 = VEQueryCriteria(variant_ids=metrics.variants_included)
    variants = variant_query_mgr.get_variants_by_task("cancer", qry1)
    assert 0 == len(
        variants.query('ALLELE_FREQUENCY<=1.0e-7'))
    variant_bm_reporter.write_summary(metrics)


def test_compute_metrics_query_params_multiple(
        variant_bm_analyzer: VEAnalyzer,
        variant_query_mgr: VEBenchmarkQueryMgr,
        sample_user_scores,
        variant_bm_reporter: VEAnalysisReporter):
    qry = VEQueryCriteria(allele_frequency_operator=">",
                          allele_frequency=1.0e-8,
                          variant_ids=sample_user_scores[:1500],
                          filter_name="Oncogene",
                          gene_symbols=['MTOR', 'PTEN'],
                          include_genes=False)
    metrics = variant_bm_analyzer.compute_metrics(
        "cancer", sample_user_scores,
        variant_effect_sources=['REVEL', 'EVE'],
        variant_query_criteria=qry, vep_min_overlap_percent=50,
        variant_vep_retention_percent=10, list_variants=True)
    variant_bm_reporter.write_summary(metrics)


def test_compute_metrics_query_params_multiple2(
        variant_bm_analyzer: VEAnalyzer,
        variant_query_mgr: VEBenchmarkQueryMgr,
        sample_user_scores,
        variant_bm_reporter: VEAnalysisReporter):
    qry = VEQueryCriteria(
                          variant_ids=sample_user_scores[:1500],
                          filter_name="Oncogene",
                          gene_symbols=['MTOR', 'PTEN'],
                          include_genes=False)
    metrics = variant_bm_analyzer.compute_metrics(
        "cancer", sample_user_scores,
        variant_effect_sources=['REVEL', 'EVE'],
        variant_query_criteria=qry, vep_min_overlap_percent=50,
        variant_vep_retention_percent=10, list_variants=True)
    variant_bm_reporter.write_summary(metrics)


