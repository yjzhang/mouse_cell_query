# TODO: aggregating queries from cellmesh and cellmarker, and perhaps other databases as well...

def aggregate_pvals(*pvals):
    """
    Given a list of lists of pvals, this returns a list of tuples:
        (list_id, inner_id, adjusted_pval).
    """
    from statsmodels.stats import multitest
    adjusted_pvals = []
    output_list = []
    for list_id, pvs in enumerate(pvals):
        rejected, pvs_corrected = multitest.fdrcorrection(pvs)
        adjusted_pvals.append(pvs_corrected)
        output_list += [(list_id, i, pv) for i, pv in enumerate(pvs_corrected)]
    output_list.sort(key=lambda x: x[2])
    return output_list

def cellmarker_cellmesh_hypergeometric_test(genes):
    """
    Returns aggregated results from running both cellmarker and cellmesh on the input genes.
    Always returns a header as the first element of the list.

    Returns:
        list of tuples (id, cell_name, FDR-adjusted pval, overlapping gens, pmids)
    """
    import cellmarker
    import cellmesh
    cellmarker_results = cellmarker.hypergeometric_test(genes, return_cl=True)
    cellmesh_results = cellmesh.hypergeometric_test(genes, return_header=True)
    header = cellmesh_results[0]
    cellmesh_results = cellmesh_results[1:]
    pvals_cellmarker = [x[2] for x in cellmarker_results]
    pvals_cellmesh = [x[2] for x in cellmesh_results]
    combined_pvals = aggregate_pvals(pvals_cellmarker, pvals_cellmesh)
    results = [cellmarker_results, cellmesh_results]
    databases = ['cellmarker', 'cellmesh']
    # combine both results set
    header = header + ['Database']
    output = [header]
    # filter results that have already been seen?
    seen_cell_types = set()
    for list_id, inner_id, pv in combined_pvals:
        row = list(results[list_id][inner_id])
        if row[1].lower in seen_cell_types or row[1].lower().strip('s') in seen_cell_types:
            continue
        seen_cell_types.add(row[1].lower())
        seen_cell_types.add(row[1].strip('s').lower())
        row[2] = pv
        row.append(databases[list_id])
        output.append(row)
    return output

def cellmarker_cellmesh_and_query(genes):
    """
    Does a "and" type of aggregation
    """
    # TODO

def cellmarker_cellmesh_tfidf_hypergeometric_test(genes):
    """
    Returns aggregated results from running both cellmarker and cellmesh on the input genes.
    Always returns a header as the first element of the list.

    Returns:
        list of tuples (id, cell_name, FDR-adjusted pval, overlapping gens, pmids)
    """
    import cellmarker
    import cellmesh
    cellmarker_results = cellmarker.hypergeometric_test(genes, return_cl=True)
    cellmesh_results = cellmesh.hypergeometric_test(genes, return_header=True)
    cellmesh_tfidf_results = cellmesh.normed_hypergeometric_test(genes, return_header=False)
    header = cellmesh_results[0]
    cellmesh_results = cellmesh_results[1:]
    pvals_cellmarker = [x[2] for x in cellmarker_results]
    pvals_cellmesh = [x[2] for x in cellmesh_results]
    pvals_cellmesh_tfidf = [x[2] for x in cellmesh_tfidf_results]
    combined_pvals = aggregate_pvals(pvals_cellmarker, pvals_cellmesh, pvals_cellmesh_tfidf)
    results = [cellmarker_results, cellmesh_results, cellmesh_tfidf_results]
    databases = ['cellmarker', 'cellmesh', 'cellmesh_tfidf']
    # combine both results set
    header = header + ['Database']
    output = [header]
    # filter results that have already been seen?
    seen_cell_types = set()
    for list_id, inner_id, pv in combined_pvals:
        row = list(results[list_id][inner_id])
        if row[1].lower in seen_cell_types or row[1].lower().strip('s') in seen_cell_types:
            continue
        seen_cell_types.add(row[1].lower())
        seen_cell_types.add(row[1].strip('s').lower())
        row[2] = pv
        row.append(databases[list_id])
        output.append(row)
    return output

if __name__ == '__main__':
    genes = '''Cd79a
Ms4a1
Cd79b
Faim3
Tnfrsf13c
2010001M09Rik
H2-Ob
Bank1
Cr2
Cd19
Cd37
Fcer2a
H2-Oa
Cd22
Fcrl1
Fcrla
Cd74
Ltb
Blk
Pou2af1
Hvcn1'''.split()
    results = cellmarker_cellmesh_hypergeometric_test(genes)
    print(results)
    print('')
    print('cellmarker + cellmesh + cellmesh_tfidf:')
    results_tfidf = cellmarker_cellmesh_tfidf_hypergeometric_test(genes)
    print(results)
