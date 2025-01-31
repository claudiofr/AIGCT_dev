"""
Data access layer methods for accessing variant repository.
Classes here provide an encapsulation layer to hide the internal
details of the repository structure.
"""

import os
import pandas as pd
from typing import List, Dict
from .util import ParameterizedSingleton
import threading
from dataclasses import dataclass, field
from .pd_util import (
    filter_dataframe_by_list,
    build_dataframe_where_clause
)
from .model import VariantFilter, VEQueryCriteria


TASK_SUBFOLDER = {
    "cancer": "cancer"
}

DATA_FOLDER = "data"
TASK_FOLDERS = [os.path.join(DATA_FOLDER, task) for task in ["cancer", "adrd",
                                                             "chd", "ddd"]]


@dataclass
class TableDef:
    folder: str
    file_name: str
    pk_columns: list[str]
    non_pk_columns: list[str]
    columns: list[str] = field(init=False),
    full_file_name: str = field(init=False)

    def __post_init__(self):
        self.columns = self.pk_columns + self.non_pk_columns
        self.full_file_name = os.path.join(self.folder, self.file_name)
    


VARIANT_PK_COLUMNS = [
    "GENOME_ASSEMBLY",
    "CHROMOSOME",
    "POSITION",
    "REFERENCE_NUCLEOTIDE",
    "ALTERNATE_NUCLEOTIDE"
]
VARIANT_NON_PK_COLUMNS = [
    "PRIOR_GENOME_ASSEMBLY",
    "PRIOR_CHROMOSOME",
    "PRIOR_POSITION",
    "PRIOR_PRIOR_GENOME_ASSEMBLY",
    "PRIOR_PRIOR_CHROMOSOME",
    "PRIOR_PRIOR_POSITION",
    "REFERENCE_AMINO_ACID",
    "ALTERNATE_AMINO_ACID",
    "AMINO_ACID_POSITION",
    "RS_DBSNP",
    "GENE_SYMBOL",
    "ENSEMBL_GENE_ID",
    "ENSEMBL_TRANSCRIPT_ID",
    "ENSEMBL_PROTEIN_ID",
    "ALLELE_FREQUENCY_SOURCE",
    "ALLELE_FREQUENCY"
]
VARIANT_TABLE_DEF = TableDef(DATA_FOLDER, "variant.csv",
                             VARIANT_PK_COLUMNS,
                             VARIANT_NON_PK_COLUMNS)

VARIANT_LABEL_NON_PK_COLUMNS = [
    "LABEL_SOURCE",
    "RAW_LABEL",
    "BINARY_LABEL"
]
VARIANT_EFFECT_LABEL_TABLE_DEF = TableDef(DATA_FOLDER,
                                          "variant_effect_label.csv",
                                          VARIANT_PK_COLUMNS,
                                          VARIANT_LABEL_NON_PK_COLUMNS)

VARIANT_EFFECT_SCORE_PK_COLUMNS = VARIANT_PK_COLUMNS + ["SCORE_SOURCE"]
VARIANT_EFFECT_SCORE_NON_PK_COLUMNS = [
    "RAW_SCORE",
    "RANK_SCORE"
]
VARIANT_EFFECT_SCORE_TABLE_DEF = TableDef(DATA_FOLDER,
                                          "variant_effect_score.csv",
                                          VARIANT_EFFECT_SCORE_PK_COLUMNS,
                                          VARIANT_EFFECT_SCORE_NON_PK_COLUMNS)

VARIANT_TASK_TABLE_DEF = TableDef(DATA_FOLDER,
                                  "variant_task.csv",
                                  ["CODE"], ["NAME", "DESCRIPTION"])

VARIANT_EFFECT_SOURCE_TABLE_DEF =\
    TableDef(DATA_FOLDER,
             "variant_effect_source.csv", ["CODE"],
             ["NAME", "SOURCE_TYPE", "DESCRIPTION"])

VARIANT_DATA_SOURCE_TABLE_DEF =\
    TableDef(DATA_FOLDER,
             "variant_data_source.csv", ["CODE"],
             ["NAME", "DESCRIPTION"])

VARIANT_FILTER_TABLE_DEF =\
    TableDef(DATA_FOLDER,
             "variant_filter.csv", ["CODE"],
             ["NAME", "DESCRIPTION", "INCLUDE_GENES",
              "INCLUDE_VARIANTS"])

VARIANT_FILTER_GENE_TABLE_DEF =\
    TableDef(DATA_FOLDER,
             "variant_filter_gene.csv", ["FILTER_CODE, GENE_SYMBOL"],
             [])

VARIANT_FILTER_VARIANT_TABLE_DEF =\
    TableDef(DATA_FOLDER,
             "variant_filter_variant.csv",
             ["FILTER_CODE"] + VARIANT_PK_COLUMNS, [])

TABLE_DEFS = {
    "VARIANT_TASK": VARIANT_TASK_TABLE_DEF,
    "VARIANT_EFFECT_SOURCE": VARIANT_EFFECT_SOURCE_TABLE_DEF,
    "VARIANT": VARIANT_TABLE_DEF,
    "VARIANT_EFFECT_LABEL": VARIANT_EFFECT_LABEL_TABLE_DEF,
    "VARIANT_DATA_SOURCE": VARIANT_DATA_SOURCE_TABLE_DEF,
    "VARIANT_EFFECT_SCORE": VARIANT_EFFECT_SCORE_TABLE_DEF,
    "VARIANT_FILTER": VARIANT_FILTER_TABLE_DEF,
    "VARIANT_FILTER_GENE": VARIANT_FILTER_GENE_TABLE_DEF,
    "VARIANT_FILTER_VARIANT": VARIANT_FILTER_VARIANT_TABLE_DEF,
}


class RepoSessionContext:

    def __init__(self, data_folder_root: str,
                 table_defs: dict[str, TableDef]):
        self._data_folder_root = data_folder_root
        self._table_defs = table_defs

    @property
    def data_folder_root(self):
        return self._data_folder_root

    def table_def(self, table_name: str):
        return self._table_defs[table_name]

    def table_file(self, table_name: str, task: str = None):
        if task:
            return os.path.join(self._data_folder_root,
                                DATA_FOLDER, task,
                                self._table_defs[table_name].file_name)
        else:
            return os.path.join(self._data_folder_root,
                                DATA_FOLDER,
                                self._table_defs[table_name].file_name)


class VariantEffectLabelCache(ParameterizedSingleton):
    """
    Caches the variant csv file in a dataframe. Implements the singleton
    pattern to ensure there is only one instance of the cached dataframe.
    We use an _init_once method rather than the normal __init__ method
    as required by the ParameterizedSingleton class.
    """

    def _init_once(self, data_folder_root: str):
        self._data_folder_root = data_folder_root
        self._lock = threading.Lock()
        self._cache = dict()

    def get_data_frame(self, task_name: str):
        if task_name not in self._cache:
            with self._lock:
                if task_name not in self._cache:
                    self._cache[task_name] = pd.read_csv(
                        os.path.join(self._data_folder_root, DATA_FOLDER,
                                     task_name,
                                     VARIANT_EFFECT_LABEL_TABLE_DEF.file_name))
        return self._cache[task_name]


class DataCache(ParameterizedSingleton):
    """
    Caches the variant csv file in a dataframe. Implements the singleton
    pattern to ensure there is only one instance of the cached dataframe.
    We use an _init_once method rather than the normal __init__ method
    as required by the ParameterizedSingleton class.
    """

    def _init_once(self, data_folder_root: str, table_def: TableDef):
        self._data_folder_root = data_folder_root
        self._table_def = table_def
        self._lock = threading.Lock()
        self._data_frame = None

    @property
    def data_frame(self):
        if self._data_frame is None:
            with self._lock:
                if self._data_frame is None:
                    self._data_frame = pd.read_csv(
                        os.path.join(self._data_folder_root,
                                     self._table_def.full_file_name))
        return self._data_frame


class TaskBasedDataCache(ParameterizedSingleton):
    """
    Caches the variant csv file in a dataframe. Implements the singleton
    pattern to ensure there is only one instance of the cached dataframe.
    We use an _init_once method rather than the normal __init__ method
    as required by the ParameterizedSingleton class.
    """

    def _init_once(self, data_folder_root: str, table_def: TableDef):
        self._data_folder_root = data_folder_root
        self._table_def = table_def
        self._cache = dict()
        self._lock = threading.Lock()

    def get_data_frame(self, task_name: str):
        if task_name not in self._cache:
            with self._lock:
                if task_name not in self._cache:
                    self._cache[task_name] = pd.read_csv(
                        os.path.join(self._data_folder_root, DATA_FOLDER,
                                     task_name,
                                     self._table_def.file_name))
        return self._cache[task_name]


class TaskDataCache(DataCache):
    """
    Caches the variant csv file in a dataframe. Implements the singleton
    pattern to ensure there is only one instance of the cached dataframe.
    We use an _init_once method rather than the normal __init__ method
    as required by the ParameterizedSingleton class.
    """

    def _init_once(self, data_folder_root: str):
        super()._init_once(data_folder_root, VARIANT_TASK_TABLE_DEF)


class TaskDataCache1(ParameterizedSingleton):
    """
    Caches the variant csv file in a dataframe. Implements the singleton
    pattern to ensure there is only one instance of the cached dataframe.
    We use an _init_once method rather than the normal __init__ method
    as required by the ParameterizedSingleton class.
    """

    def _init_once(self, data_folder_root: str, table_def: TableDef):
        self._data_folder_root = data_folder_root
        self._table_def = table_def
        self._lock = threading.Lock()

    def get_data_frame(self, task_name: str) -> pd.DataFrame:
        if self._data_frame is None:
            with self._lock:
                if self._data_frame is None:
                    self._data_frame = pd.read_csv(
                        os.path.join(self._data_folder_root,
                                     self._table_def.full_file_name))
        return self._data_frame


class VariantEffectScoreCache(TaskBasedDataCache):
    """
    Caches the variant csv file in a dataframe. Implements the singleton
    pattern to ensure there is only one instance of the cached dataframe.
    We use an _init_once method rather than the normal __init__ method
    as required by the ParameterizedSingleton class.
    """

    def _init_once(self, data_folder_root: str):
        super()._init_once(data_folder_root, VARIANT_EFFECT_SCORE_TABLE_DEF)


class VariantCache(DataCache):
    """
    Caches the variant csv file in a dataframe. Implements the singleton
    pattern to ensure there is only one instance of the cached dataframe.
    We use an _init_once method rather than the normal __init__ method
    as required by the ParameterizedSingleton class.
    """

    def _init_once(self, data_folder_root: str):
        super()._init_once(data_folder_root, VARIANT_TABLE_DEF)


class VariantCache1(ParameterizedSingleton):
    """
    Caches the variant csv file in a dataframe. Implements the singleton
    pattern to ensure there is only one instance of the cached dataframe.
    We use an _init_once method rather than the normal __init__ method
    as required by the ParameterizedSingleton class.
    """

    def _init_once(self, data_folder_root: str):
        self._data_folder_root = data_folder_root
        self._lock = threading.Lock()

    @property
    def data_frame(self):
        if self._data_frame is None:
            with self._lock:
                if self._data_frame is None:
                    self._data_frame = pd.read_csv(
                        os.path.join(self._data_folder_root,
                                     VARIANT_TABLE_DEF.full_file_name))
        return self._data_frame


class VariantTaskCache(DataCache):
    """
    Caches the variant csv file in a dataframe. Implements the singleton
    pattern to ensure there is only one instance of the cached dataframe.
    We use an _init_once method rather than the normal __init__ method
    as required by the ParameterizedSingleton class.
    """

    def _init_once(self, data_folder_root: str):
        super()._init_once(data_folder_root, VARIANT_TASK_TABLE_DEF)


class VariantTaskCache1(ParameterizedSingleton):
    """
    Caches the variant csv file in a dataframe. Implements the singleton
    pattern to ensure there is only one instance of the cached dataframe.
    We use an _init_once method rather than the normal __init__ method
    as required by the ParameterizedSingleton class.
    """

    def _init_once(self, data_folder_root: str):
        self._data_folder_root = data_folder_root
        self._lock = threading.Lock()

    @property
    def data_frame(self):
        if self._data_frame is None:
            with self._lock:
                if self._data_frame is None:
                    self._data_frame = pd.read_csv(
                        os.path.join(self._data_folder_root,
                                     VARIANT_TASK_TABLE_DEF.full_file_name))
        return self._data_frame


class VariantEffectSourceCache(DataCache):
    """
    Caches the variant csv file in a dataframe. Implements the singleton
    pattern to ensure there is only one instance of the cached dataframe.
    We use an _init_once method rather than the normal __init__ method
    as required by the ParameterizedSingleton class.
    """

    def _init_once(self, data_folder_root: str):
        super()._init_once(data_folder_root, VARIANT_EFFECT_SOURCE_TABLE_DEF)


class VariantFilterCache(ParameterizedSingleton):
    """
    Caches the variant csv file in a dataframe. Implements the singleton
    pattern to ensure there is only one instance of the cached dataframe.
    We use an _init_once method rather than the normal __init__ method
    as required by the ParameterizedSingleton class.
    """

    def _init_once(self, data_folder_root: str):
        self._data_folder_root = data_folder_root
        self._lock = threading.Lock()
        self._cache = dict()

    def get_data_frames(self, task_name: str) -> Dict:
        if task_name not in self._cache:
            with self._lock:
                if task_name not in self._cache:
                    folder = os.path.join(self._data_folder_root, DATA_FOLDER,
                                          task_name)
                    cache_dict = dict()
                    cache_dict["filter_df"] = pd.read_csv(
                        os.path.join(folder,
                                     VARIANT_FILTER_TABLE_DEF.file_name))
                    cache_dict["filter_gene_df"] = pd.read_csv(
                        os.path.join(folder,
                                     VARIANT_FILTER_GENE_TABLE_DEF.file_name))
                    cache_dict["filter_variant_df"] = pd.read_csv(
                        os.path.join(
                            folder,
                            VARIANT_FILTER_VARIANT_TABLE_DEF.file_name))
                    self._cache[task_name] = cache_dict
        return self._cache[task_name]


class VariantEffectSourceRepository:
    """
    Caches the variant csv file in a dataframe. Implements the singleton
    pattern to ensure there is only one instance of the cached dataframe.
    We use an _init_once method rather than the normal __init__ method
    as required by the ParameterizedSingleton class.
    """

    def __init__(self, session_context: RepoSessionContext,
                 variant_effect_score_repo):
        self._cache = VariantEffectSourceCache(
            session_context.data_folder_root)
        self._variant_effect_score_repo = variant_effect_score_repo

    def get_all(self) -> pd.DataFrame:
        return self._cache.data_frame.copy(deep=True)

    def get_by_task(self, task_name: str) -> pd.DataFrame:
        score_sources = self._variant_effect_score_repo.get_all_by_task(
            task_name)['SCORE_SOURCE'].unique()
        source_df = self._cache.data_frame.copy(deep=True)
        return source_df[source_df['CODE'] in score_sources]

    def get_by_code(self, codes: list[str]) -> pd.DataFrame:
        return filter_dataframe_by_list(
            self._cache.data_frame, codes,
            VARIANT_EFFECT_SOURCE_TABLE_DEF.pk_columns)


class VariantTaskRepository:
    """
    Caches the variant csv file in a dataframe. Implements the singleton
    pattern to ensure there is only one instance of the cached dataframe.
    We use an _init_once method rather than the normal __init__ method
    as required by the ParameterizedSingleton class.
    """

    def __init__(self, session_context: RepoSessionContext):
        self._cache = VariantTaskCache(session_context.data_folder_root)

    def get_all(self) -> pd.DataFrame:
        return self._cache.data_frame.copy(deep=True)


class VariantFilterRepository:
    """
    Caches the variant csv file in a dataframe. Implements the singleton
    pattern to ensure there is only one instance of the cached dataframe.
    We use an _init_once method rather than the normal __init__ method
    as required by the ParameterizedSingleton class.
    """

    def __init__(self, session_context: RepoSessionContext):
        self._cache = VariantFilterCache(session_context.data_folder_root)

    def get_by_task(self, task_name: str) -> Dict:
        return self._cache.get_data_frames(task_name)

    def get_by_task_filter_name(
            self, task_name: str, filter_name: str) -> VariantFilter:
        filter_dfs = self.get_by_task(task_name)
        filter = filter_dfs["filter_df"].query(f"NAME == '{filter_name}'")
        if len(filter) == 0:
            return None
        filter = filter.iloc[0]
        filter_code = filter["CODE"]
        filter_genes = filter_dfs["filter_gene_df"].query(
            f"FILTER_CODE == '{filter_code}'")
        filter_variants = filter_dfs["filter_variant_df"].query(
            f"FILTER_CODE == '{filter_code}'")
        if len(filter_genes) == 0 and len(filter_variants) == 0:
            return None
        return VariantFilter(filter, filter_genes, filter_variants)


def query_by_filter(query_df: pd.DataFrame,
                    filter: pd.Series,
                    filter_gene_df: pd.DataFrame,
                    filter_variant_df: pd.DataFrame) -> pd.DataFrame:

    if len(filter_gene_df) > 0:
        query_df = filter_dataframe_by_list(query_df, filter_gene_df,
                                            ['GENE_SYMBOL'], None,
                                            filter["INCLUDE_GENES"] ==
                                            'Y')
    if len(filter_variant_df) > 0:
        query_df = filter_dataframe_by_list(query_df, filter_variant_df,
                                            VARIANT_PK_COLUMNS, None,
                                            filter["INCLUDE_VARIANTS"] == 'Y')
    return query_df


class VariantRepository:
    """
    Caches the variant csv file in a dataframe. Implements the singleton
    pattern to ensure there is only one instance of the cached dataframe.
    We use an _init_once method rather than the normal __init__ method
    as required by the ParameterizedSingleton class.
    """

    def __init__(self, session_context: RepoSessionContext):
        self._cache = VariantCache(session_context.data_folder_root)

    def get_all(self) -> pd.DataFrame:
        return self._cache.data_frame.copy(deep=True)

    def get(self, qry: VEQueryCriteria) -> pd.DataFrame:
        """
        Fetches variants. The optional parameters are filter criteria used to
        limit the set of variants returned.

        Parameters
        ----------
        gene_symbols : list, optional
            List of gene symbols
        include_genes : bool, optional
            If gene_symbols is provided, indicates whether to limit variants
            to associated with those gene_symbols to exclude variants
            associated with the gene_symbols.
        variant_ids : DataFrame, optional
            List of variant ids
        include_variant_ids : bool, optional
            If variant_ids is provided, indicates whether to limit variants
            to the variant_ids provided or to fetch all variants but those
            in variant_ids
        allele_frequency_operator : str, optional
            If allele_frequency is provided, this is one of "=", ">",
            "<", ">=", "<=", "!=". i.e. limit variants to those whose
            allele_frequency is equal to, greater than, etc. the
            allele_frequency.
        allele_frequency : float, optional
            Used in conjunction to allele_frequency_operator to limit variants
            to those meeting a certain allele_frequency criteria.
        filter_name : str, optional
            The name of a system filter that can be used to limit the variants
            returned.
        """
        where_clause = build_dataframe_where_clause(
            {"ALLELE_FREQUENCY": [qry.allele_frequency_operator,
                                  qry.allele_frequency]})
        variant_df = self._cache.data_frame
        if where_clause != "":
            variant_df = variant_df.query(where_clause)
        if qry.gene_symbols is not None:
            variant_df = filter_dataframe_by_list(variant_df, qry.gene_symbols,
                                                  ['GENE_SYMBOL'],
                                                  in_list=qry.include_genes)
        if qry.variant_ids is not None:
            variant_df = filter_dataframe_by_list(variant_df, qry.variant_ids,
                                                  VARIANT_PK_COLUMNS,
                                                  qry.column_name_map,
                                                  qry.include_variant_ids)
        return variant_df


class VariantEffectLabelRepository:
    """
    Caches the variant csv file in a dataframe. Implements the singleton
    pattern to ensure there is only one instance of the cached dataframe.
    We use an _init_once method rather than the normal __init__ method
    as required by the ParameterizedSingleton class.
    """

    def __init__(self, session_context: RepoSessionContext,
                 variant_repo: VariantRepository,
                 filter_repo: VariantFilterRepository):
        self._cache = VariantEffectLabelCache(session_context.data_folder_root)
        self._filter_repo = filter_repo
        self._variant_repo = variant_repo

    def get_all_by_task(self, task_name: str) -> pd.DataFrame:
        label_df = self._cache.get_data_frame(task_name)
        variant_df = self._variant_repo.get(variant_ids=label_df[
            VARIANT_PK_COLUMNS], include_variant_ids=True)
        return label_df.merge(variant_df, on=VARIANT_PK_COLUMNS,
                              how="inner")

    def get(self, task_name: str,
            qry: VEQueryCriteria = None) -> pd.DataFrame:
        """
        Fetches variants. The optional parameters are filter criteria used to
        limit the set of variants returned.

        Parameters
        ----------
        task_name : str
            Retrieve variants associated with this task.
        gene_symbols : list, optional
            List of gene symbols
        include_genes : bool, optional
            If gene_symbols is provided, indicates whether to limit variants
            to associated with those gene_symbols to exclude variants
            associated with the gene_symbols.
        variant_ids : DataFrame, optional
            List of variant ids
        include_variant_ids : bool, optional
            If variant_ids is provided, indicates whether to limit variants
            to the variant_ids provided or to fetch all variants but those
            in variant_ids
        allele_frequency_operator : str, optional
            If allele_frequency is provided, this is one of "eq", "gt",
            "lt", "ge", "le". i.e. limit variants to those whose
            allele_frequency is equal to, greater than, etc. the
            allele_frequency.
        allele_frequency : float, optional
            Used in conjunction to allele_frequency_operator to limit variants
            to those meeting a certain allele_frequency criteria.
        filter_name : str, optional
            The name of a system filter that can be used to limit the variants
            returned.
        """
        label_df = self._cache.get_data_frame(task_name)
        if qry is not None:
            if qry.variant_ids is not None and len(qry.variant_ids) > 0:
                label_df = filter_dataframe_by_list(
                    label_df, qry.variant_ids,
                    VARIANT_PK_COLUMNS,
                    in_list=qry.include_variant_ids)
            variant_df = self._variant_repo.get(qry)
        else:
            variant_df = self._variant_repo.get_all()
        merge_df = label_df.merge(variant_df, how="inner",
                                  on=VARIANT_PK_COLUMNS)
        if qry is not None and qry.filter_name is not None:
            filter_dfs = \
                self._filter_repo.get_by_task_filter_name(task_name,
                                                          qry.filter_name)
            if filter_dfs is None:
                raise Exception("Invalid filter name: " + qry.filter_name)
            merge_df = query_by_filter(merge_df, filter_dfs.filter,
                                       filter_dfs.filter_genes,
                                       filter_dfs.filter_variants)
        return merge_df


class VariantEffectScoreRepository:
    """
    Caches the variant csv file in a dataframe. Implements the singleton
    pattern to ensure there is only one instance of the cached dataframe.
    We use an _init_once method rather than the normal __init__ method
    as required by the ParameterizedSingleton class.
    """

    def __init__(self, session_context: RepoSessionContext,
                 variant_repo: VariantRepository,
                 filter_repo: VariantFilterRepository):
        self._cache = VariantEffectScoreCache(session_context.data_folder_root)
        self._filter_repo = filter_repo
        self._variant_repo = variant_repo

    def get_all_by_task(self, task_name: str) -> pd.DataFrame:
        label_df = self._cache.get_data_frame(task_name)
        query_criteria = VEQueryCriteria(variant_ids=label_df[
            VARIANT_PK_COLUMNS], include_variant_ids=True)
        variant_df = self._variant_repo.get(query_criteria)
        return label_df.merge(variant_df, on=VARIANT_PK_COLUMNS,
                              how="inner")

    def get(self, task_name: str,
            variant_effect_sources: list[str] | str = None,
            include_variant_effect_sources: bool = True,
            qry: VEQueryCriteria = None) -> pd.DataFrame:
        """
        Fetches variants. The optional parameters are filter criteria used to
        limit the set of variants returned.

        Parameters
        ----------
        task_name : str
            Retrieve variants associated with this task.
        gene_symbols : list, optional
            List of gene symbols
        include_genes : bool, optional
            If gene_symbols is provided, indicates whether to limit variants
            to associated with those gene_symbols to exclude variants
            associated with the gene_symbols.
        variant_ids : DataFrame, optional
            List of variant ids
        include_variant_ids : bool, optional
            If variant_ids is provided, indicates whether to limit variants
            to the variant_ids provided or to fetch all variants but those
            in variant_ids
        allele_frequency_operator : str, optional
            If allele_frequency is provided, this is one of "eq", "gt",
            "lt", "ge", "le". i.e. limit variants to those whose
            allele_frequency is equal to, greater than, etc. the
            allele_frequency.
        allele_frequency : float, optional
            Used in conjunction to allele_frequency_operator to limit variants
            to those meeting a certain allele_frequency criteria.
        filter_name : str, optional
            The name of a system filter that can be used to limit the variants
            returned.
        """
        score_df = self._cache.get_data_frame(task_name)
        if (variant_effect_sources is not None and
                len(variant_effect_sources) > 0):
            score_df = filter_dataframe_by_list(score_df,
                                                variant_effect_sources,
                                                "SCORE_SOURCE",
                                                None,
                                                include_variant_effect_sources)
        if qry is None:
            return score_df
        if qry.variant_ids is not None:
            score_df = filter_dataframe_by_list(score_df, qry.variant_ids,
                                                VARIANT_PK_COLUMNS,
                                                qry.column_name_map,
                                                qry.include_variant_ids)
        variant_df = self._variant_repo.get(qry)
        merge_df = score_df.merge(variant_df, how="inner",
                                  on=VARIANT_PK_COLUMNS)
        if qry.filter_name is not None:
            filter_dfs = \
                self._filter_repo.get_by_task_filter_name(task_name,
                                                          qry.filter_name)
            if filter_dfs is None:
                raise Exception("Invalid filter name: " + qry.filter_name)
            merge_df = query_by_filter(merge_df, filter_dfs.filter,
                                       filter_dfs.filter_genes,
                                       filter_dfs.filter_variants)
        return merge_df[VARIANT_EFFECT_SCORE_TABLE_DEF.columns]


