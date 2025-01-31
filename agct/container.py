"""
Class to simulate a proper Dependency Injection container.
It could be reimplemented in the future if we decide to use
a proper one. The interface, however, would remain the same.
"""

from .repository import (
    VariantEffectScoreRepository,
    VariantEffectLabelRepository,
    RepoSessionContext,
    VariantFilterRepository,
    VariantRepository,
    VariantEffectSourceRepository,
    VariantTaskRepository,
    TABLE_DEFS
)
from .analyzer import VEAnalyzer
from .query import VEBenchmarkQueryMgr
from .reporter import VEAnalysisReporter
from .plotter import VEAnalysisPlotter
from .exporter import VEAnalysisExporter
from .util import Config

import yaml
import os


class VEBenchmarkContainer:

    def __init__(self, app_root: str = "."):
        with (open(os.path.join(app_root, "config", "config.yaml"), "r") as
              config_file):
            self.config = Config(yaml.safe_load(config_file))
        self._repo_session_context = RepoSessionContext(
            self.config.repository.root_dir, TABLE_DEFS)
        self._variant_task_repo = VariantTaskRepository(
            self._repo_session_context)
        self._variant_repo = VariantRepository(self._repo_session_context)
        self._variant_filter_repo = VariantFilterRepository(
            self._repo_session_context)
        self._label_repo = VariantEffectLabelRepository(
            self._repo_session_context,
            self._variant_repo,
            self._variant_filter_repo)
        self._score_repo = VariantEffectScoreRepository(
            self._repo_session_context,
            self._variant_repo,
            self._variant_filter_repo)
        self._variant_effect_source_repo = VariantEffectSourceRepository(
            self._repo_session_context,
            self._score_repo)
        self._variant_filter_repo = VariantFilterRepository(
            self._repo_session_context
        )
        self._analyzer = VEAnalyzer(
            self._score_repo,
            self._label_repo,
            self._variant_effect_source_repo)
        self._query_mgr = VEBenchmarkQueryMgr(self._label_repo,
                                              self._variant_repo,
                                              self._variant_task_repo,
                                              self._variant_effect_source_repo,
                                              self._score_repo,
                                              self._variant_filter_repo)
        self._reporter = VEAnalysisReporter()
        self._plotter = VEAnalysisPlotter(self.config.plot)
        self._exporter = VEAnalysisExporter()

    @property
    def analyzer(self):
        return self._analyzer

    @property
    def query_mgr(self):
        return self._query_mgr

    @property
    def reporter(self):
        return self._reporter

    @property
    def plotter(self):
        return self._plotter

    @property
    def exporter(self):
        return self._exporter
