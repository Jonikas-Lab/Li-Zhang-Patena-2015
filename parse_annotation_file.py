#!/usr/bin/env python2.7

# standard library
from __future__ import division
import sys
import unittest
import os
from collections import defaultdict


HEADER_FIELDS_v4 = ['Phytozome transcript name', 
                    'PFAM', 'Panther', 'KOG', 'KEGG ec', 'KEGG Orthology', 
                    'best arabidopsis TAIR10 hit name', 'best arabidopsis TAIR10 hit symbol', 'best arabidopsis TAIR10 hit defline', 
                    'best rice hit name', 'best rice hit symbol', 'best rice hit defline']

HEADER_FIELDS_v5 = ['Phytozome internal transcript ID', 
                    'Phytozome gene locus name', 'Phytozome transcript name', 'Phytozome protein name', 
                    'PFAM', 'Panther', 'KOG', 'KEGG ec', 'KEGG Orthology', 'Gene Ontology terms', 
                    'best arabidopsis TAIR10 hit name', 'best arabidopsis TAIR10 hit symbol', 'best arabidopsis TAIR10 hit defline', 
                    'best rice hit name', 'best rice hit symbol', 'best rice hit defline']


def parse_gene_annotation_file(gene_annotation_filename, standard_Phytozome_file=False, header_fields=None, gene_ID_column=0, 
                               genes_start_with=None, delete_columns=None, remove_empty_columns=True, strip_gene_fields_start=None,
                               pad_with_empty_fields=True, ignore_comments=False, verbosity_level=1):
    """ Parse tab-separated gene annotation file; return (gene:annotation_list, header_list) tuple.

    If standard_Phytozome_file is non-zero, use special case headers for standard Phytozome files:
        v4 annotation file (Creinhardtii_169_annotation_info.txt) if value is 4, v5 (Creinhardtii_236_readme.txt) if 5.
        (those files don't have included headers - they're hardcoded. And the formats sometimes change, so this isn't a guarantee.)
    Otherwise, header can be provided with header_fields - if None, the first line of the file is assumed to be the header.

    Use column gene_ID_column to determine gene IDs; optionally shorten the gene name by truncating it starting with the 
     strip_gene_fields_start value if found (if not None) - e.g. if value is '.t', Cre01.g123450.t2.1 would become Cre01.g123450.
    If genes_start_with is not None, make sure all gene IDs start with it.

    If pad_with_empty_fields is True, pad shorter lines to the max length. 
    If ignore_comments is True, skip lines starting with #.

    Print some info/warnings to stdout depending on verbosity_level (0 - nothing, 1 - some, 2 - max).
    """
    if not os.path.lexists(gene_annotation_filename):
        raise Exception("Couldn't find the %s gene annotation file!"%gene_annotation_filename)
    if verbosity_level>0:
        print " *** Parsing file %s for gene annotation info..."%gene_annotation_filename

    if delete_columns is None:
        delete_columns = []

    # special cases for the standard Phytozome annotation files
    if standard_Phytozome_file:
        strip_gene_fields_start = '.t'
        if standard_Phytozome_file == 4:
            gene_ID_column = 0
            header_fields = list(HEADER_FIELDS_v4)  # using list() to make a separate copy, since it'll be modified!)
            delete_columns = []
        elif standard_Phytozome_file == 5:
            gene_ID_column = 1
            header_fields = list(HEADER_FIELDS_v5)  # using list() to make a separate copy, since it'll be modified!)
            delete_columns = [0,2,3]
        else:
            raise Exception("Invalid value for standard_Phytozome_file type arg! %s given, 4/5 accepted."%standard_Phytozome_file)

    ### Parse the whole file into lists of tab-separated fields
    #    (could change to a generator, but most of the data has to stay in memory anyway in a different format, so probably no point)
    #    (and we need special treatment for the first line which may or may not be a header line...)
    data_by_row = []
    for line in open(gene_annotation_filename):
        if ignore_comments and line[0]=='#':    continue
        fields = line.strip().split('\t')
        data_by_row.append(fields)
    if verbosity_level>0:
        print "Parsed %s lines"%len(data_by_row)
        
    # if header is given, use list() to make a separate copy, since it will be modified
    header_fields = list(header_fields)
    # if header not given, assume the first line is a header
    if header_fields is None:
        header_fields = data_by_row[0]
        del data_by_row[0]
        if verbosity_level>0:
            print "Assuming the first line is a header: %s"%'\t'.join(header_fields)
    # if any of the other lines doesn't start with a Cre* gene ID, fail!
    if genes_start_with is not None:
        for row in data_by_row:
            if not row[0].startswith(genes_start_with):
                raise Exception("Can't parse file %s - found line that doesn't start "%gene_annotation_filename
                                +"with a %s gene ID!\n  \"%s\""%(genes_start_with, '\t'.join(row)))

    # check that all the data lengths line up (if they don't, don't throw an error
    mismatched_lengths = False
    data_lengths = set([len(row) for row in data_by_row])
    if not len(data_lengths)==1:
        if verbosity_level>1 or (not pad_with_empty_fields and verbosity_level>0):
            print "Warning: not all data rows have the same length! Lengths found: %s"%list(data_lengths)
        mismatched_lengths = True
    if len(header_fields) not in data_lengths:
        if verbosity_level>1 or (not pad_with_empty_fields and verbosity_level>0):
            print("Warning: header has a different number of fields than the data! "
                  +"Header length: %s. Data lengths: %s"%(len(header_fields),list(data_lengths)))
        mismatched_lengths = True
    if len(data_lengths)>1 and pad_with_empty_fields:
        max_length = max(max(data_lengths), len(header_fields))
        if verbosity_level>0:
            print "Data field numbers vary between rows - padding all lower-length data rows to length %s"%max_length
        for row in data_by_row:
            if len(row)<max_length:
                row += ['' for x in range(max_length-len(row))]
        if len(header_fields)<max_length:
            header_fields += ['?' for x in range(max_length-len(header_fields))]
        mismatched_lengths = False

    # remove empty columns (only if all the data lengths match!)
    if remove_empty_columns and mismatched_lengths:
        raise Exception("Cannot remove empty columns if different lines have different numbers of columns!")
    if remove_empty_columns and not mismatched_lengths:
        data_length = len(header_fields)
        columns_to_remove = []
        for pos in range(data_length):
            values = set([row[pos] for row in data_by_row])
            if len(values)==1:
                value = values.pop()
                if value.strip()=='':
                    if verbosity_level>0:
                        if header_fields:   print "Column %s (%s) is always empty - removing it."%(pos+1, header_fields[pos])
                        else:               print "Column %s is always empty - removing it."%(pos+1)
                    columns_to_remove.append(pos)
                else:
                    if verbosity_level>0:
                        print "Note: all the values in column %s are the same! (%s)"%(pos+1, value)
        for pos in sorted(columns_to_remove, reverse=True):
            if header_fields:  del header_fields[pos]
            for row in data_by_row: 
                del row[pos]

    data_length = len(header_fields)
    assert data_length >= max(len(row) for row in data_by_row)

    ### convert the list-format data into a by-gene dictionary
    # we're removing one column that was the gene ID, and some other columns - that's the final length
    final_data_length = data_length - 1 - len(delete_columns)
    data_by_gene = defaultdict(lambda: [set() for x in range(final_data_length)])
    for data in data_by_row:
        gene = data[gene_ID_column]
        for col in sorted(delete_columns + [gene_ID_column], reverse=True):
            del data[col]
        if strip_gene_fields_start is not None:
            gene = gene.split(strip_gene_fields_start)[0]
        # We frequently get multiple lines, for different transcripts!  Just concatenate all of them.
        for N,field in enumerate(data):
            data_by_gene[gene][N].add(field)

    # remove the first word from the header, since it should be "gene ID" or such; 
    #  change spaces to underscores in header fields for readability
    if header_fields:  
        for col in sorted(delete_columns + [gene_ID_column], reverse=True):
            del header_fields[col]
        header_fields = [s.replace(' ','_') for s in header_fields]

    # At the end, change the sets to comma-separated strings, and also remove empty strings
    for gene,data in data_by_gene.items():
        data_by_gene[gene] = [' | '.join([f for f in fields if f.strip()]) for fields in data]

    if verbosity_level>0:
        print " *** DONE Parsing gene annotation file"
    return dict(data_by_gene), header_fields


class Testing(unittest.TestCase):
    """ Runs unit-tests for this module. """

    def test__parse_gene_annotation_file_v5(self):
        # do it twice just to make sure the results are the same - earlier I had a bug where they weren't!
        for x in (1,2):
            data, header = parse_gene_annotation_file('test_data/INPUT_annotation_file_v5.txt', standard_Phytozome_file=5, 
                                                      verbosity_level=0)
            assert header == ['PFAM', 'Panther', 'KOG', 'KEGG_ec', 'KEGG_Orthology', 'Gene_Ontology_terms', 
                              'best_arabidopsis_TAIR10_hit_name', 'best_arabidopsis_TAIR10_hit_symbol', 
                              'best_arabidopsis_TAIR10_hit_defline']
            assert set(data.keys()) == set(['Cre01.g000150', 'Cre01.g000100', 'g8754'])
            assert data['Cre01.g000100'] == ['' for x in header]
            assert data['g8754'] == ['', 'PTHR13693:SF10,PTHR13693', '', '', '', '', 'AT5G04620.2', 'ATBIOF,BIOF', 'biotin F']
            assert data['Cre01.g000150'] == ['PF02535', 'PTHR11040,PTHR11040:SF30', 'KOG1558', 'KEGGstuff', 'KEGGstuff2', 
                                             'GO:0046873,GO:0030001,GO:0055085,GO:0016020', 'AT2G04032.1', 'ZIP7', 
                                             'zinc transporter 7 precursor']

    # MAYBE-TODO add a v4 unit-test too, and for other formats?


if __name__=='__main__':
    """ If module is run directly, run tests. """
    print "This is a module for import by other programs - it doesn't do anything on its own.  Running tests..."
    unittest.main()
