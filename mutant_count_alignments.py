#!/usr/bin/env python2.7
""" Take a deepseq alignment file; group the reads into insertional mutants.
Output a line-per-mutant file containing position info, gene annotation data (optional), total/perfect read count, number of distinct sequences, and some of the sequences/counts (optional), and a summary of read/mutant counts etc. 
Also output a line-per-gene file containing gene ID/name/position/annotation and the number and read-counts of mutants that were inserted into that gene (optional, NOT IMPLEMENTED YET).
Output files are in simple tab-separated plaintext format, with all header/summary lines starting with #.

Grouping reads into mutants is currently done by alignment position (other options such as sequence clustering or grouping 1bp-adjacent positions together may be implemented later). 

The input file should be a SAM-format deepseq alignment file created by bowtie, novoalign, or other deepseq aligner programs (tested mainly on bowtie), with optionally a metadata file created by my deepseq_alignment_wrapper.py or deepseq_preprocessing_wrapper.py scripts.
The program assumes the SAM file contains only unique matches (i.e. each read was reported as aligning to at most one genomic location).
It also only deals with single-end alignments at the moment.

 -- Weronika Patena, Jonikas Lab, Carnegie Institution, 2011

USAGE: mutant_count_alignments.py [options] infile1 [infile2 infile3 ...] outfile """

# basic library
import sys, os, time
import unittest
import pickle
# other packages
import HTSeq
# my modules
from general_utilities import write_header_data
from testing_utilities import run_functional_tests
import mutant_analysis_classes
import mutant_Carette   # TODO this should be optional...

def do_test_run():
    """ Test run: run script on test infile, compare output to reference file."""
    test_folder = "test_data_v0"
    aln_infile0 = "%s/INPUT_alignment0_old-format.sam"%test_folder
    aln_infile1 = "%s/INPUT_alignment1_genomic-unique.sam"%test_folder
    aln_infile2 = "%s/INPUT_alignment2_for-genes.sam"%test_folder
    aln_infile3 = "%s/INPUT_alignment3_for-merging.sam"%test_folder
    gff_genefile = "%s/INPUT_gene-data-1_all-cases.gff3"%test_folder
    dataset_to_remove = "%s/INPUT_mutants_to_remove.txt"%test_folder

    test_runs = [
                 ('cassette-end-5prime', "-e 5prime -r forward -n3 -L", [aln_infile1]),
                 ('cassette-end-3prime', "-e 3prime -r forward -n3 -L", [aln_infile1]),
                 ('read-direction-reverse', "-r reverse -e 5prime -n3 -L", [aln_infile1]),
                 ('unknown-as-match', "--treat_unknown_as_match -e 5prime -r forward -n3 -L", [aln_infile1]),
                 ('dont-count-cassette', "-l -e 5prime -r forward -n3 -L", [aln_infile1]),
                 ('ignore-cassette', "-c -e 5prime -r forward -n3 -L", [aln_infile1]),
                 ('separate-cassette', "-C -e 5prime -r forward -n3 -L", [aln_infile1]),
                 ('sorted-by-count', "-o read_count -e 5prime -r forward -n3 -L", [aln_infile1]),
                 ('with-gene-info_merged', "-e 5prime -r forward -g %s -n0"%gff_genefile, [aln_infile2]),
                 ('with-gene-info_unmerged', "-B -e 5prime -r forward -g %s -n0"%gff_genefile, [aln_infile2]),
                 ('multiple-infiles', "-e 5prime -r forward -n0 -L", [aln_infile1,aln_infile2]),
                 ('dont-merge-tandems', "-n0", [aln_infile3]),
                 ('merge-adjacent-none', "-n0 -Q -Y0", [aln_infile3]),
                 ('merge-adjacent1-r3', "-MQ -D1 -w3 -Y0 -n0", [aln_infile3]),
                 ('merge-adjacent1-r1', "-MQ -D1 -w1 -Y0 -n0", [aln_infile3]),
                 ('merge-adjacent2-r3', "-MQ -D2 -w3 -Y0 -n0", [aln_infile3]), 
                 ('remove-from-other-all', "-x %s -n0"%dataset_to_remove, [aln_infile2]), 
                 ('remove-from-other-min4', "-x %s -z4 -n0"%dataset_to_remove, [aln_infile2]), 
                 ('remove-from-other-perfect', "-x %s -p -z4 -n0"%dataset_to_remove, [aln_infile2]),
                 ('remove-not-other-all', "-X %s -n0"%dataset_to_remove, [aln_infile2]), 
                 ('remove-not-other-min4', "-X %s -Z4 -n0"%dataset_to_remove, [aln_infile2]), 
                 ('remove-not-other-perfect', "-X %s -P -Z4 -n0"%dataset_to_remove, [aln_infile2]),
                 ('old-infile-format', "-e 5prime -r forward -n3 -L", [aln_infile0]),
                ]
    # TODO add run-test for removing data from multiple files?
    # MAYBE-TODO add run-test for a metadata file with 5' and 3' read counts?
    # MAYBE-TODO add run-tests for other mutant-merging options?  But they all have pretty good unit-tests.
    # MAYBE-TODO add run-test for --gene_annotation_file?
    # MAYBE-TODO add run-test for --input_collapsed_to_unique?  Or is that unit-tested already?

    # convert tests into (testname, arg_and_infile_string) format, adding the options that are always used
    test_names_and_args = [('count-aln__'+testname, test_args+' -q '+' '.join(infiles)) 
                           for testname,test_args,infiles in test_runs]

    parser = define_option_parser()
    argument_converter = lambda parser,options,args: (args[:-1], args[-1], options)
    return run_functional_tests(test_names_and_args, parser, main, test_folder, 
                                argument_converter=argument_converter, append_to_outfilenames='.txt') 


class Testing(unittest.TestCase):
    """ Unit-tests this module. """

    def test__(self):
        print "NO UNIT-TESTS FOR THIS MODULE"


######### Main function code #########

def define_option_parser():
    """ Populates and returns an optparse option parser object, with __doc__ as the usage string."""
    from optparse import OptionParser
    parser = OptionParser(__doc__)

    # taken:     aAbBcC-De---g-h---jJ--lLmMn-o-pPqQr---tTuUvVwWxXyYzZ  
    # free:      ------d--EfF-G-HiI--kK-----N-O-----RsS--------------  

    ### test options
    parser.add_option('-t','--test_functionality', action='store_true', default=False, 
                      help="Run the built-in unit test suite (ignores all other options/arguments; default %default).")
    parser.add_option('-T','--test_run', action='store_true', default=False, 
                      help="Run on a test input file, check output against reference files. "
                          + "Ignores all other options/arguments. (default %default).")

    ### functionality options
    parser.add_option('-e', '--read_cassette_end', choices=mutant_analysis_classes.SEQ_ENDS, default='5prime', 
                      metavar='|'.join(mutant_analysis_classes.SEQ_ENDS), 
                      help="Which end of the cassette are the sequenced reads from? (default %default).")
    parser.add_option('-r','--read_direction', choices=mutant_analysis_classes.SEQ_DIRECTIONS, default='forward',
                      metavar='|'.join(mutant_analysis_classes.SEQ_DIRECTIONS), 
                      help="Is the read in the forward or reverse direction compared to the cassette? (default %default).")
    parser.add_option('--Carette', action="store_true", default=False, 
                      help="Is this data from the Carette protocol instead of the MmeI one? (default %default).")

    parser.add_option('-D', '--adjacent_max_distance', type='int', default=1, metavar='N',
                      help="Count/merge adjacent mutants only if they're at most N bases distant; (default %default)")
    parser.add_option('-M', '--merge_adjacent_mutants', action="store_true", default=False, 
                      help="Merge adjacent mutants if they satisfy the -W/-w constraints (default %default)")
    parser.add_option('-W', '--merge_adjacent_leave_mutants', default='auto', metavar='N',
                      help="For -M: how many mutants to leave unmerged: a number, or 'auto' (based on opposite-strand "
                          +"adjacent mutants); ignored if -w option is given. (default %default)")
    parser.add_option('-w', '--merge_adjacent_count_ratio', type='int', default=None, metavar='K',
                      help="For -M: merge mutants only if one has at least K times fewer reads than the other (default %default)")
    parser.add_option('-Q', '--merge_opposite_tandem_mutants', action="store_true", default=False,
                      help="Merge opposite-strand same-position mutants (tail-to-tail-tandems) "
                          +"if they satisfy the -Y/-y contraints. (default %default)")
    parser.add_option('-Y', '--merge_opposite_leave_mutants', default='auto', metavar='N',
                      help="For -Q: how many mutants to leave unmerged: a number, or 'auto' (based on opposite-strand "
                          +"adjacent mutants); ignored if -y option is given, (default %default)")
    parser.add_option('-y', '--merge_opposite_count_ratio', type='int', default=None, metavar='K',
                      help="For -Q: merge mutants only if the readcounts are within Kx of each other. (default %default)")
    parser.add_option('--merge_mutant_choice_method', choices=['by_ratio','random'], default='by_ratio', metavar='by_ratio|random',
                      help="For -M and -Q: if using non-zero -W/-Y (and NOT using -w/-y), should the mutants to merge be chosen "
                          +"randomly, or by highest/lowest readcount ratio? (default %default)")

    parser.add_option('-j', '--merge_in_cassette', action='store_true', default=False, 
                      help="For adjacent and tandem merging/counts: include mutants in cassette (default %default)")
    parser.add_option('-J', '--merge_in_other_chrom', action='store_true', default=False, 
                      help="For adjacent and tandem merging/counts: include mutants in non-cassette non-nuclear chromosomes "
                          +"(like chloroplast and mitochondrial) (default %default)")

    parser.add_option('-x', '--remove_mutants_from_file', metavar='FILE1,FILE2',
                      help='Remove all mutants present in FILE1,FILE2,etc from the datasets (see -z/-p for read count cutoff).'
                          +' Any number of files is allowed - must be comma-separated without spaces.')
    parser.add_option('-z', '--remove_from_file_readcount_min', type='int', default=1, metavar='M',
                      help='When applying -x, only remove mutants with at least N reads in FILE (default %default).')
    parser.add_option('-p', '--remove_from_file_min_is_perfect', action='store_true', default=False,
                      help='When applying -x with -z M, compare M to perfect readcount, not total. (default %default).')

    parser.add_option('-X', '--remove_mutants_not_from_file', metavar='FILE',
                      help='Remove all mutants NOT present in FILE from the datasets (see -Z/-P for read count cutoff).')
    parser.add_option('-Z', '--remove_not_from_file_readcount_min', type='int', default=1, metavar='M',
                      help='When applying -X, only remove mutants with at least N reads in FILE (default %default).')
    parser.add_option('-P', '--remove_not_from_file_min_is_perfect', action='store_true', default=False,
                      help='When applying -X with -z M, compare M to perfect readcount, not total. (default %default).')

    parser.add_option('-c', '--ignore_cassette', action='store_true', default=False,
                      help="Ignore reads aligning to cassette (just print total count in the header as removed) "
                          +"(default %default)")
    parser.add_option('-C', '--separate_cassette', action='store_true', default=False,
                      help="Like -c, but also add a *_cassette.txt file with ONLY reads aligning to cassette "
                          +"(and print total non-cassette count in the header of that file as removed) (default %default)")
    parser.add_option('-l', '--dont_count_cassette', action='store_true', default=False, 
                      help="Don't give separate cassette read/mutant totals in the header; (default %default)")
    parser.add_option('-L', '--dont_count_other', action='store_true', default=False, 
                      help="Don't give separate read/mutant totals in the header for 'strange' chromosomes "
                          +"(not cassette and not named chromosome* or scaffold*; (default %default)")

    # extremely minor functionality options, do we even care??
    parser.add_option('--treat_unknown_as_match', action="store_true", default=False, 
                      help="When counting perfect reads, treat undefined alignment regions as matches (default %default)")
    parser.add_option('--dont_treat_unknown_as_match', action="store_false", dest='treat_unknown_as_match',
                      help="Turn -u off.")
    # MAYBE-TODO add user-provided mutation-count cutoffs like in old deepseq_count_alignments.py, instead of just all reads and perfet reads?   Currently useless, since we're only allowing one mutation in bowtie.  parser.add_option('-m', '--mutation_cutoffs', default="1,3,10", metavar="<comma-separated-int-list>")

    ### input options
    parser.add_option('-u','--input_collapsed_to_unique', action='store_true', default=False, 
                      help="Use to get correct original total read counts if the data was collapsed to unique sequences "
                          +"using fastx_collapser before alignment (default %default).")
    parser.add_option('-U','--input_not_collapsed_to_unique', action='store_false', dest="input_collapsed_to_unique", 
                      help="Turn -c off.")
    parser.add_option('-m', '--input_metadata_file', default='AUTO', metavar='FILE', 
                      help="File containing preprocessing and alignment metadata (scripts/options used etc). "
                          +"Can be a filename, AUTO for <infile_basename>_info.txt (warning will be raised if not found), "
                          +"or NONE to not look for a metadata file at all. Default %default.")

    ### gene-finding and gene-annotation options 
    parser.add_option('-g', '--gene_position_reference_file', default=None, metavar='FILE', 
                      help="File to use to look up gene IDs based on chromosomal location (default %default)")
    parser.add_option('--detailed_gene_features', action="store_true", default=True,
                      help="Find out what part of the gene (UTR,intron,exon) a mutant hit, based on the -g file "
                          +"(default %default). May take a lot of memory - increase --N_detail_run_groups to fix that.")
    parser.add_option('--no_detailed_gene_features', action="store_false", dest='detailed_gene_features',
                      help="Turns --detailed_gene_features off.")
    parser.add_option('--N_detail_run_groups', type="int", default=5, metavar='N', 
                      help="How many passes to split reading the detailed_gene_features into (default %default) "
                          +"- may take a lot of memory (and CPU) if read in a single pass; too many passes waste CPU.")
    # MAYBE-TODO add a "flank" option (with variable size), to catch mutants that are in the flanks of genes? Do we care?
    # MAYBE-TODO add a "negative flank" option (with variable size), to ignore mutants that are in the start/end of genes?

    parser.add_option('-A', '--gene_annotation_file', default=None, metavar='FILE', 
                      help="Tab-separated file to use to look up gene names/descriptions from IDs (default %default)")
    parser.add_option('-a', '--annotation_file_standard_type', type='int', default=None, metavar='V', 
                      help="Use if file provided in -A is a standard Phytozome annotation file (missing a header): "
                          +"use value 4 for the one with chlamy v4.3 genome, or value 5 for v5 genome (default %default)")

    ### output format options
    parser.add_option('-n', '--N_sequences_per_group', type='int', default=2, metavar='N', 
                      help="How many most common sequences should be shown per group? (default %default)")
    parser.add_option('-o', '--sort_data_key', choices=['position','read_count','none'], default='position', 
                      metavar='position|read_count|none', help="Sort the output data: by alignment position, read count, "
                         +"or don't sort at all (default %default) - sorting may be slow for large datasets!")
    parser.add_option('-B', '--dont_merge_boundary_features', action='store_true', default=False,
                      help="In the summary, count all feature-boundary cases separately instead of together "
                          +"(default %default)")

    parser.add_option('-V', '--verbosity_level', action="store_true", default=1, 
                      help="How much information to print to STDOUT: 0 - nothing, 1 - summary only, "
                          +"2 - summary and progress reports. (Default %default).")
    parser.add_option('-q', '--quiet', action="store_const", const=0, dest='verbosity_level', help="Equivalent to -V 0.")
    parser.add_option('-v', '--verbose', action="store_const", const=2, dest='verbosity_level', help="Equivalent to -V 2.")

    return parser


def get_info_from_metadata_files(infiles, input_metadata_file, cassette_end, verbosity_level):
    """ Parse metadata files to get wrong-format and unaligned/multiple read counts (or 'unknown' if cannot be determined). 

    Returns a tuple of the following counts: discarded, wrong_start, no_cassette, non_aligned, unaligned, multiple_aligned.
    Any of them can be 'unknown', meaning it could not be determined. 
    It looks for the information in the metadata file for each infile; if it can't find or parse the metadata file for 
     any infile, the final count is unknown.
    Special treatment for the non_aligned value for old-format files, in which that count is not in metadata at all.
    """
    # if the option specified no metadata files, total discarded readcount cannot be determined
    if input_metadata_file == 'NONE':
        return 'unknown', 'unknown', 'unknown', 'unknown', 'unknown', 'unknown', 'unknown'
    # make sure the -m option has a value that will work with the number of infiles
    if len(infiles)>1 and not input_metadata_file=='AUTO':
        print "Warning: when multiple input files are given, the -m option must be NONE or AUTO - ignoring other value."
        return 'unknown', 'unknown', 'unknown', 'unknown', 'unknown', 'unknown', 'unknown'

    # get the read counts for each infile; only return the total at the end, if all values are found
    counts_per_file = {'discarded':[], 'wrong_start':[], 'no_cassette':[], 'other_end':[], 'unaligned':[], 'multiple_aligned':[]}
    file_formats = []

    for infile in infiles:
        file_format = '?'
        # take the actual input_metadata_file value as the file name, or infer it from the infile name if 'AUTO'
        if input_metadata_file == 'AUTO':
            infile_base = os.path.splitext(infile)[0]
            if infile_base.endswith('_genomic-unique'):     infile_base = infile_base[:-len('_genomic-unique')]
            curr_input_metadata_file = infile_base + '_info.txt'
            if verbosity_level>1:  
                print 'Automatically determining metadata input file name: %s'%curr_input_metadata_file
        else:
            curr_input_metadata_file = input_metadata_file
            if verbosity_level>1:  
                print 'Metadata input file name provided in options: %s'%curr_input_metadata_file
        # if file is missing, readcounts cannot be determined
        if not os.path.exists(curr_input_metadata_file):
            if verbosity_level>0:
                print 'Warning: metadata input file %s not found! Proceeding without it.'%curr_input_metadata_file
            counts_per_file['discarded'].append('unknown')
            counts_per_file['wrong_start'].append('unknown')
            counts_per_file['no_cassette'].append('unknown')
            counts_per_file['unaligned'].append('unknown')
            counts_per_file['multiple_aligned'].append('unknown')
            file_formats.append(file_format)
            continue

        # go through the metadata file to find the lines with various pieces of useful information
        #  (if the line isn't found, the related readcount cannot be determined)

        ### total discarded read count (there are two possible formats of that line, old and new)
        for line in open(curr_input_metadata_file):
            if line.startswith('## final "bad" reads'):         # new discarded-read line format
                counts_per_file['discarded'].append(int(line.split(':\t')[1].split(' (')[0]))
                file_format = 'new'
                break
            if line.startswith('## reads removed: '):           # old discarded-read line format
                counts_per_file['discarded'].append(int(line.split()[3]))
                file_format = 'old'
                break
        else:   # in a for-else loop the else is executed if the for wasn't ended with break, i.e. the line wasn't found
            if verbosity_level>0:
                print("Warning: metadata input file %s didn't contain discarded read count line! "%curr_input_metadata_file
                      +"Proceeding without it.")
            counts_per_file['discarded'].append('unknown')

        file_formats.append(file_format)
        ### all the other information only shows up in new-format files - in old-format just assume unknown
        #  (except the other_end value, which was always 0 in old-format because old-format lacked that dual-end processing)
        if file_format=='old':
                counts_per_file['wrong_start'].append('unknown')
                counts_per_file['no_cassette'].append('unknown')
                counts_per_file['other_end'].append(0)
                counts_per_file['unaligned'].append('unknown')
                counts_per_file['multiple_aligned'].append('unknown')
    
        else:
            ### wrong-start discarded read count (new-format only)
            for line in open(curr_input_metadata_file):
                if line.startswith('#  "bad" read count (wrong-start)'):
                    counts_per_file['wrong_start'].append(int(line.split(':\t')[1].split(' (')[0]))
                    break
            else:
                if verbosity_level>0:
                    print("Warning: metadata file %s didn't contain wrong-start readcount line! "%curr_input_metadata_file
                          +"Proceeding without it.")
                counts_per_file['wrong_start'].append('unknown')

            ### no-cassette discarded read count (new-format only)
            for line in open(curr_input_metadata_file):
                if line.startswith('#  "bad" read count (no-cassette)'):
                    counts_per_file['no_cassette'].append(int(line.split(':\t')[1].split(' (')[0]))
                    break
            else:
                if verbosity_level>0:
                    print("Warning: metadata file %s didn't contain no-cassette readcount line! "%curr_input_metadata_file
                          +"Proceeding without it.")
                counts_per_file['no_cassette'].append('unknown')

            ### other-end read count (new-format only, optional)
            # this line is not present if there weren't any other-end reads, so if it's not, we'll save 0 instead of unknown, 
            #  since it means a 0, not an unknown value (and old-format files didn't have this functionality so it's always 0 there)
            other_end = [x for x in mutant_analysis_classes.SEQ_ENDS if not x==cassette_end][0].replace("prime","'")
            for line in open(curr_input_metadata_file):
                if line.startswith('# "good" %s read count after cassette stripping'%other_end):
                    counts_per_file['other_end'].append(int(line.split(':\t')[1].split(' (')[0]))
                    break
            else:
                if verbosity_level>0:
                    print("Warning: metadata file %s didn't contain other-end readcount line! "%curr_input_metadata_file
                          +"Proceeding without it.")
                counts_per_file['other_end'].append(0)
            # for the purposes of this_end readcounts, other_end readcounts should be included in the discarded count 
            #  (UNLESS the discarded or other_end count is unknown, then just set discarded to that.)
            if 'unknown' in (counts_per_file['discarded'][-1], counts_per_file['other_end'][-1]):
                counts_per_file['discarded'][-1] = 'unknown'
            else:
                counts_per_file['discarded'][-1] += counts_per_file['other_end'][-1]

            ### unaligned read count (new-format only)
            for line in open(curr_input_metadata_file):
                if line.startswith('# unaligned:  '):
                    counts_per_file['unaligned'].append(int(line.split(':  ')[1].split(' (')[0]))
                    break
            else:
                if verbosity_level>0:
                    print("Warning: metadata file %s didn't contain unaligned readcount line! "%curr_input_metadata_file
                          +"Proceeding without it.")
                counts_per_file['unaligned'].append('unknown')

            ### multiple-aligned read count (new-format only)
            for line in open(curr_input_metadata_file):
                if line.startswith('# multiple-genomic:  '):
                    counts_per_file['multiple_aligned'].append(int(line.split(':  ')[1].split(' (')[0]))
                    break
            else:
                if verbosity_level>0:
                    print("Warning: metadata file %s "%curr_input_metadata_file 
                          +"didn't contain genomic-multiple readcount line! Proceeding without it.")
                counts_per_file['multiple_aligned'].append('unknown')

            # MAYBE-TODO also grab the "good" readcounts (from preprocessing and alignment) and make sure they match 
            #  the numbers we see in the dataset?

    # The calculation for total non-aligned is a bit complicated, because:
    #  for new-format files that information comes from the metadata (and split into unaligned/multiple subcategories), 
    #  but for the old-format files it comes from the data infile itself (without the category split. 
    # So for each old-format file, the general nonaligned count can be assumed to be 0 (and the real number will be 
    #  (added during actual infile processing), but the unaligned/multiple counts stay unknown
    counts_per_file['total_non_aligned'] = []
    for file_format, unaligned_count, multiple_count in zip(file_formats, counts_per_file['unaligned'], 
                                                            counts_per_file['multiple_aligned']):
        if file_format=='old':                                  counts_per_file['total_non_aligned'].append(0)
        elif 'unknown' in (unaligned_count, multiple_count):    counts_per_file['total_non_aligned'].append('unknown')
        else:                                                   counts_per_file['total_non_aligned'].append(\
                                                                                        unaligned_count+multiple_count)

    # do sums for each count category: of the length of the list is wrong or unknown is on the list, return unknown:
    #  if some metadata files were missing the discarded counts line, total discarded readcount cannot be determined
    #  if the discarded read counts for all infiles were found, return the sum as the total discarded read count
    total_counts = {}
    for key, list_of_counts in counts_per_file.iteritems():
        if len(list_of_counts)<len(infiles) or ('unknown' in list_of_counts):
            if verbosity_level>0:
                print "Warning: %s read count not found for some files! Ignoring all values, keeping 'unknown'."%key
            total_counts[key] = 'unknown'
        else:
            total_counts[key] = sum(list_of_counts)

    return (total_counts['discarded'], total_counts['wrong_start'], total_counts['no_cassette'], total_counts['other_end'], 
            total_counts['total_non_aligned'], total_counts['unaligned'], total_counts['multiple_aligned'])


# TODO merge this with a similar code bit in mutant_join_datasets.py to minimize code duplication
def save_dataset_files(dataset, outfile, verbosity_level=0, if_pickle=True, count_cassette=True, count_other=True, 
                       merge_boundary_features=True, sort_data_by='position', N_sequences_per_mutant=5, options="N/A"):
    """ Print summary and data to output file; optionally print summary to stdout; optionally pickle dataset to picklefile. 
    
    The options argument is only used to be printed in the header to make it clear how the file was generated - 
     it should be the applicable optparse options object if there is one, or a text message otherwise.
    """
    # print summary info to stdout if desired
    if verbosity_level>1: print "\nDATA SUMMARY:"
    if verbosity_level>0: dataset.print_summary(merge_boundary_features=merge_boundary_features, 
                                                count_cassette=count_cassette, count_other=count_other)
    # print full data to outfile
    if verbosity_level>1: print "printing output - time %s."%time.ctime()
    with open(outfile,'w') as OUTFILE:
        write_header_data(OUTFILE,options)
        OUTFILE.write("### SUMMARY:\n")
        dataset.print_summary(OUTFILE, line_prefix="#  ", header_prefix="## ", merge_boundary_features=merge_boundary_features,
                              count_cassette = count_cassette, count_other=count_other)
        OUTFILE.write("### HEADER AND DATA:\n")
        dataset.print_data(OUTPUT=OUTFILE, sort_data_by=sort_data_by, N_sequences=N_sequences_per_mutant, 
                           header_line=True, header_prefix='# ')
    # print pickled dataset to picklefile, if desired
    if if_pickle:
        outfile_basename = os.path.splitext(outfile)[0]
        pickled_outfile = outfile_basename + '.pickle'
        with open(pickled_outfile,'w') as PICKLEFILE:
            pickle.dump(dataset, PICKLEFILE, 0)


def main(infiles, outfile, options):
    """ Run the main functionality of the module (see module docstring for more information), excluding testing.
    The options argument should be generated by an optparse parser.
    """
    ### parse/process/reformat some options
    options.ignore_cassette |= options.separate_cassette
    options.count_cassette = not options.dont_count_cassette
    options.count_other = not options.dont_count_other
    options.merge_boundary_features = not options.dont_merge_boundary_features
    # MAYBE-TODO change outfile to a folder?  Since it'll have three things in it now...
    outfile_basename = os.path.splitext(outfile)[0]
    mutant_merging_outfile = outfile_basename + '_merging-info.txt'
    # MAYBE-TODO let -C take an optional argument to put the cassette files elsewhere?
    cassette_outfile = outfile_basename + '_cassette.txt'
    cassette_merging_outfile = outfile_basename + '_cassette_merging-info.txt'

    ### generate empty alignment set object with basic read position/orientation properties defined by options
    if options.Carette:     dataset_class = mutant_Carette.Insertional_mutant_pool_dataset_Carette
    else:                   dataset_class = mutant_analysis_classes.Insertional_mutant_pool_dataset
    all_alignment_data = dataset_class(options.read_cassette_end, options.read_direction=='reverse')
    if options.separate_cassette:
        cassette_alignment_data = dataset_class(options.read_cassette_end, options.read_direction=='reverse')
    # MAYBE-TODO refactor the whole bunch of "if options.separate_cassette:" clauses to avoid code duplication?

    ### parse preprocessing/alignment metadata file to get discarded/not-aligned/etc readcounts, pass to all_alignment_data
    #   (all_alignment_data initializes them to 'unkown', so if file is not given or can't be found/parsed, do nothing)
    N_discarded, N_wrong_start, N_no_cassette, N_other_end, N_non_aligned, N_unaligned, N_multiple = \
            get_info_from_metadata_files(infiles, options.input_metadata_file, options.read_cassette_end, options.verbosity_level)
    if 'unknown' not in (N_wrong_start, N_no_cassette):
        assert N_discarded == N_wrong_start+N_no_cassette+N_other_end,\
                "Discarded subtotals don't add up to discarded total! %s+%s+%s != %s"%(N_wrong_start,N_no_cassette,
                                                                                       N_other_end,N_discarded)
    all_alignment_data.summary.add_discarded_reads(N_discarded, N_wrong_start, N_no_cassette, N_other_end)
    if options.separate_cassette:
        cassette_alignment_data.summary.add_discarded_reads(N_discarded, N_wrong_start, N_no_cassette, N_other_end)
    all_alignment_data.summary.add_nonaligned_reads(N_non_aligned, N_unaligned, N_multiple)
    if options.separate_cassette:
        cassette_alignment_data.summary.add_nonaligned_reads(N_non_aligned, N_unaligned, N_multiple)
    # MAYBE-TODO also get the final total number of reads from the metadata infile and make sure it's the same 
    #   as the number of processed reads I get from all_alignment_data.print_summary()?

    ### parse input file and store data - the add_alignment_reader_to_data function here does pretty much all the work!
    for infile in infiles:
        # if this is a new-style *_genomic-unique.sam file and has a matching *_cassette.sam file, parse that file too
        part_infiles = [infile]
        if infile.endswith('_genomic-unique.sam'):
            cassette_file = infile[:-len('_genomic-unique.sam')] + '_cassette.sam'
            if os.path.exists(cassette_file):
                part_infiles.append(cassette_file)
        for part_infile in part_infiles:
            # initialize a parser for the SAM infile
            if options.verbosity_level>1: 
                print "parsing input file %s - time %s."%(part_infile, time.ctime())
            infile_reader = HTSeq.SAM_Reader(part_infile)
            # fill the new alignment set object with data from the infile parser
            all_alignment_data.add_alignment_reader_to_data(infile_reader, 
                                        uncollapse_read_counts = options.input_collapsed_to_unique, 
                                        ignore_cassette = options.ignore_cassette, cassette_only = False, 
                                        treat_unknown_as_match = options.treat_unknown_as_match)
            if options.separate_cassette:
                cassette_alignment_data.add_alignment_reader_to_data(infile_reader, 
                                        uncollapse_read_counts = options.input_collapsed_to_unique, 
                                        ignore_cassette = False, cassette_only = True, 
                                        treat_unknown_as_match = options.treat_unknown_as_match)

    ### optionally remove mutants based on another dataset - BEFORE adjacent mutant counting/merging
    # remove mutants that ARE present in another file (also do it for cassette mutants if those are separate)
    if options.remove_mutants_from_file:
        for other_file in options.remove_mutants_from_file.split(','):
            other_dataset = mutant_analysis_classes.read_mutant_file(other_file)
            all_alignment_data.remove_mutants_in_other_dataset(other_dataset, 
                     readcount_min=options.remove_from_file_readcount_min, perfect_reads=options.remove_from_file_min_is_perfect)
            if options.separate_cassette:
                cassette_alignment_data.remove_mutants_in_other_dataset(other_dataset, 
                     readcount_min=options.remove_from_file_readcount_min, perfect_reads=options.remove_from_file_min_is_perfect)
    # remove mutants that are NOT present in another file (also do it for cassette mutants if those are separate)
    # TODO should I implement using multiple files here too?
    if options.remove_mutants_not_from_file:
        other_dataset = mutant_analysis_classes.read_mutant_file(options.remove_mutants_not_from_file)
        all_alignment_data.remove_mutants_not_in_other_dataset(other_dataset, 
                 readcount_min=options.remove_not_from_file_readcount_min, perfect_reads=options.remove_not_from_file_min_is_perfect)
        if options.separate_cassette:
            cassette_alignment_data.remove_mutants_not_in_other_dataset(other_dataset, 
                 readcount_min=options.remove_not_from_file_readcount_min, perfect_reads=options.remove_not_from_file_min_is_perfect)

    ### optionally merge some mutant categories
    with open(mutant_merging_outfile, 'w') as MERGEFILE:
      with open(cassette_merging_outfile, 'w') as CASSETTE_MERGEFILE:
        # 1) adjacent same-strand mutants (since they're probably just artifacts of indels during deepseq/PCR)
        if options.merge_adjacent_mutants: 
            try:                                                leave_N_mutants = int(options.merge_adjacent_leave_mutants)
            except ValueError:                                  leave_N_mutants = options.merge_adjacent_leave_mutants
            if options.merge_adjacent_count_ratio is not None:  leave_N_mutants = 'use_ratio'
            all_alignment_data.merge_adjacent_mutants(merge_max_distance=options.adjacent_max_distance, 
                      leave_N_mutants=leave_N_mutants, min_count_ratio = options.merge_adjacent_count_ratio, 
                      leave_method=options.merge_mutant_choice_method, merge_cassette_chromosomes = options.merge_in_cassette, 
                      merge_other_chromosomes = options.merge_in_other_chrom, OUTPUT = MERGEFILE)
            if options.merge_in_cassette and options.separate_cassette:
                cassette_alignment_data.merge_adjacent_mutants(merge_max_distance = options.adjacent_max_distance, 
                          leave_N_mutants=leave_N_mutants, min_count_ratio = options.merge_adjacent_count_ratio, 
                          leave_method=options.merge_mutant_choice_method, merge_cassette_chromosomes = True, 
                          merge_other_chromosomes = False, OUTPUT = CASSETTE_MERGEFILE)
        # 2) opposite-strand same-position mutants (since they're probably just tail-to-tail cassette tandems)
        if options.merge_opposite_tandem_mutants: 
            try:                                                leave_N_mutants = int(options.merge_opposite_leave_mutants)
            except ValueError:                                  leave_N_mutants = options.merge_opposite_leave_mutants
            if options.merge_opposite_count_ratio is not None:  leave_N_mutants = 'use_ratio'
            all_alignment_data.merge_opposite_tandem_mutants(leave_N_mutants=leave_N_mutants, 
                          max_count_ratio = options.merge_opposite_count_ratio, leave_method=options.merge_mutant_choice_method, 
                          merge_cassette_chromosomes = options.merge_in_cassette, 
                          merge_other_chromosomes = options.merge_in_other_chrom, OUTPUT = MERGEFILE)
            if options.merge_in_cassette and options.separate_cassette:
                cassette_alignment_data.merge_opposite_tandem_mutants(leave_N_mutants=leave_N_mutants, 
                              max_count_ratio = options.merge_opposite_count_ratio, leave_method=options.merge_mutant_choice_method, 
                              merge_cassette_chromosomes = True, merge_other_chromosomes = False, OUTPUT = CASSETTE_MERGEFILE)
        # 3) count adjacent mutants, even if not doing any merging
        # Actually there's always a count done after each merge, but that count doesn't print the details to anything, 
        #  so just run it again here with detail-printing - MAYBE-TODO find a more efficient way instead of counting twice?
        all_alignment_data.count_adjacent_mutants(max_distance_to_print = options.adjacent_max_distance, 
                                  max_distance_to_count = 10000, count_cassette_chromosomes = options.merge_in_cassette, 
                                  count_other_chromosomes = options.merge_in_other_chrom, OUTPUT = MERGEFILE)
        if options.separate_cassette:
            cassette_alignment_data.count_adjacent_mutants(max_distance_to_print = options.adjacent_max_distance, 
                                  max_distance_to_count = 10, count_cassette_chromosomes = True, 
                                  count_other_chromosomes = False, OUTPUT = CASSETTE_MERGEFILE)
        # MAYBE-TODO make an option for max_distance_to_count?  I'm using a much lower one for cassette because it's so dense.
    # since there's no optional as/with statement, just remove the cassette_merging_outfile if unwanted
    if not options.separate_cassette:
        os.remove(cassette_merging_outfile)

    ### optionally parse gene position/info files and look up the genes for each mutant in the data
    if options.gene_position_reference_file is not None:
        genefile = options.gene_position_reference_file
        if options.verbosity_level>1: print "adding genes from file %s to mutant data - time %s."%(genefile, time.ctime())
        all_alignment_data.find_genes_for_mutants(genefile, detailed_features=options.detailed_gene_features, 
                                                  N_run_groups=options.N_detail_run_groups, 
                                                  verbosity_level=options.verbosity_level)

        # if we have gene info, optionally also add annotation
        if options.gene_annotation_file:
            if options.verbosity_level>1: 
                print "adding gene annotation from file %s - time %s."%(options.gene_annotation_file, time.ctime())
            all_alignment_data.add_gene_annotation(options.gene_annotation_file, 
                       if_standard_Phytozome_file=options.annotation_file_standard_type, print_info=(options.verbosity_level >= 2))

    ### output data to files
    save_dataset_files(all_alignment_data, outfile, options.verbosity_level, True, options.count_cassette, options.count_other, 
                         options.merge_boundary_features, options.sort_data_key, options.N_sequences_per_group, options)
    # TODO write some info about all the other files that go with this one (pickle, merging-info, *cassette*)

    if options.separate_cassette:
        save_dataset_files(cassette_alignment_data, cassette_outfile, 0, True, False, False, 
                             options.merge_boundary_features,  options.sort_data_key, options.N_sequences_per_group, options)
        # TODO probably write some extra bit of info about this being the cassette-file to MAINFILENAME


if __name__ == "__main__":
    """ Allows both running and importing of this file. """

    parser = define_option_parser()
    (options, args) = parser.parse_args()

    # if ran with -t option, do unit tests and quit
    if options.test_functionality:
        print("*** You used the -t option - ignoring all other options/arguments, running the built-in test suite. ***")
        print("\n * unit-tests for the mutant_analysis_classes.py module")
        # to run tests for another file, have to use TextTestRunner, not unittest.main -  make a test suite with 
        #   autodetection of all tests (see http://docs.python.org/library/unittest.html#unittest.TestLoader)
        test_suite_1 = unittest.defaultTestLoader.loadTestsFromModule(mutant_analysis_classes)
        unittest.TextTestRunner(verbosity=1).run(test_suite_1)
        # to run tests for current module, just run unittest.main, passing it only the filename 
        #   (by default it takes all of sys.argv and complains about options/arguments it can't recognize)
        print("\n * unit-tests for this module (%s)"%sys.argv[0])
        unittest.main(argv=[sys.argv[0]])   # unittest.main automatically runs sys.exit()

    if options.test_run:
        print("*** You used the -T option - ignoring all other options and running the built-in example test runs. ***")
        test_result = do_test_run()
        sys.exit(test_result)

    # otherwise parse the arguments and run main function
    if len(args)<2:
        parser.print_help()
        sys.exit("\nError: at least one infile and exactly one outfile are required!")
    outfile = args[-1]
    infiles = args[:-1]

    main(infiles, outfile, options)

