"""CellPhoneDB analysis with caller-provided resources; no local path setup or database download.
The external library requires resource/output arguments and may perform I/O itself.
"""

def infer_interactions(database, metadata, counts, destination):
    from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
    return cpdb_statistical_analysis_method.call(
        cpdb_file_path=database,
        meta_file_path=metadata,
        counts_file_path=counts,
        counts_data="hgnc_symbol",
        score_interactions=True,
        iterations=1000,
        threshold=0.1,
        threads=25,
        debug_seed=42,
        result_precision=3,
        pvalue=0.05,
        subsampling=False,
        subsampling_log=False,
        subsampling_num_pc=100,
        subsampling_num_cells=1000,
        separator="|",
        debug=False,
        output_path=destination,
        output_suffix="prad",
    )
