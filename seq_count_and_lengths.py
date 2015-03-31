#! /usr/bin/env python2.7
"""
Given any number of fasta or fastq files, print a list of seq lengths and counts of seqs with that length (or just the total seq count).  The file type is detected by filename extension.  Uses the Biopython package.
 --Weronika Patena, June 2011
"""

# standard library
import sys
from collections import defaultdict
import unittest
# other packages
from Bio import SeqIO
# my modules
from basic_seq_utilities import get_seq_count_from_collapsed_header, check_fasta_fastq_format
from general_utilities import add_dicts_of_ints


def seq_count_and_lengths(seq_iterator, count_only=False, input_collapsed_to_unique=False):
    """ Given an iterator over sequences, return N_seqs and a seq_len:seq_count dict (empty if count_only). 

    Sequence length is determined by len(seq) - will fail if len() doesn't work on the elements of seq_iterator. 
    If input_collapsed_to_unique, decode the read count from seq header instead of counting each seq as 1, 
     using basic_seq_utilities.get_seq_count_from_collapsed_header (see docstring for that).
    """
    total_count = 0
    seqlen_counter = defaultdict(lambda: 0)
    for seq in seq_iterator: 
        N_seqs = get_seq_count_from_collapsed_header(seq.name) if input_collapsed_to_unique else 1
        total_count += N_seqs
        if not count_only:
            seqlen_counter[len(seq)] += N_seqs
    return total_count, dict(seqlen_counter)


def _format_lengths(seqlen_dict, include_zeros=False, verbosity=1):
    """ Given a length:seqcount dictionary, format it for printing, return list of lines. """
    output_lines = []
    lengths = seqlen_dict.keys()
    if include_zeros:   lengths = range(min(lengths),max(lengths)+1)
    else:               lengths.sort()
    if verbosity>0:     
        output_lines.append("Total %s seqs\n"%sum(seqlen_dict.values()))
    if verbosity>0:     
        output_lines.append("length\tseq count\n")
    for l in lengths:
        output_lines.append("%s\t%s\n"%(l,seqlen_dict[l]))
    return output_lines


def main(infiles, total_seq_number_only=False, input_collapsed_to_unique=False, 
         include_zeros=False, verbosity=1, OUTPUT=sys.stdout):
    """ Given a list of fastq/fasta files, return total seq number, a length:N dict and formatted info (optionally print).
    
    If total_seq_number_only is True, only return/print total seq count.
    If input_collapsed_to_unique is True, program assumes infile was preprocessed with fastx_collapser, 
     and attempts to give original pre-collapsing seq_counts (based on headers).
    If include_zeros is False (default), only print non-zero seq counts; if True, print seq counts for 
     all lengths between min and max length, even if they're 0.
    Verbosity: if >1, print filetype and seqcount for each input file; if 0, don't print header or summary.
    Prints to stdout by default; to print to file, pass open file object as OUTPUT; to suppress printing, pass None."""

    # a counter with a default value of 0
    total_seqcount, total_seqlen_dict = 0, {}
    formatted_output = []
    # add the numbers from each file to total_seqlen_dict
    for infile in infiles:
        # detect filetype based on extension
        #  MAYBE-TODO add command-line options that force the format to fasta/fastq instead of checking by extension?
        seq_format = check_fasta_fastq_format(infile, verbosity>1)
        # note: just using plain "fastq" quality encoding, because we're not dealing with qualities so it doesn't matter
        with open(infile) as INFILE:
            file_seqcount, file_seqlen_dict = seq_count_and_lengths(SeqIO.parse(INFILE, seq_format), 
                                                                    total_seq_number_only, input_collapsed_to_unique)
        total_seqcount += file_seqcount
        total_seqlen_dict = add_dicts_of_ints(total_seqlen_dict, file_seqlen_dict)

    # format and print (optionally) and return the output
    if total_seq_number_only:
        formatted_output.append("Total %s seqs\n"%total_seqcount)
    else:
        formatted_output += _format_lengths(total_seqlen_dict, include_zeros, verbosity)
    if not OUTPUT is None:
        for line in formatted_output:   OUTPUT.write(line)
    return total_seqcount, total_seqlen_dict, formatted_output
            

class Testing(unittest.TestCase):
    """ Unit-tests this module. """

    def test__seq_count_and_lengths__empty(self):
        for empty_iterator in [[], (), set(), iter([])]:
            assert seq_count_and_lengths(empty_iterator, count_only=False) == (0, {})
            assert seq_count_and_lengths(empty_iterator, count_only=True) == (0, {})

    def test__seq_count_and_lengths__single_length(self):
        for LEN in [0, 1, 3, 10, 100, 1000, 10000]:
            for N in [1, 3, 10, 100, 1000, 10000]:
                assert seq_count_and_lengths(['A'*LEN]*N, count_only=False) == (N, {LEN:N})
                assert seq_count_and_lengths(['A'*LEN]*N, count_only=True)  == (N, {})

    def test__seq_count_and_lengths__multi_length(self):
        assert seq_count_and_lengths(['AAA', 'AAA', 'A'], count_only=False) == (3, {1:1, 3:2})
        assert seq_count_and_lengths(['AAA', 'AAA', 'A'], count_only=True)  == (3, {})

    def test__seq_count_and_lengths__input_collapsed_to_unique(self):
        class Fake_seq(): 
            def __init__(self, name, length):   self.name, self.length = name, length
            def __len__(self):                  return self.length
        # sequences without a .name attribute or with one that doesn't have an encoded readcount cause an exception
        self.assertRaises(AttributeError, seq_count_and_lengths, ['AAA'], input_collapsed_to_unique=True)
        self.assertRaises(ValueError, seq_count_and_lengths, [Fake_seq('AAA', 3)], input_collapsed_to_unique=True)
        # sequences with a .name and a length work
        seqs = [Fake_seq('1-10', 1), Fake_seq('2-20', 2), Fake_seq('3-30', 3)]
        assert seq_count_and_lengths(seqs, input_collapsed_to_unique=False)  == (3, {1:1, 2:1, 3:1})
        assert seq_count_and_lengths(seqs, input_collapsed_to_unique=True)  == (60, {1:10, 2:20, 3:30})

    # MAYBE-TODO add unit-test for _format_lengths?

    def test__main(self):
        # FUNCTION SIGNATURE: main(infiles, total_seq_number_only=False, input_collapsed_to_unique=False, 
        #                          include_zeros=False, verbosity=1, OUTPUT=sys.stdout)
        quiet = dict(verbosity=0, OUTPUT=None)
        test_fa, test_fq, with_counts = "_test_inputs/test.fa", "_test_inputs/test.fq", "_test_inputs/with-counts.fa"
        # wrong input file format (unrecognized extension)
        self.assertRaises(ValueError, main, ["_test_inputs/textcmp_file2.txt"], **quiet)
        # no input files = empty output
        assert main([], **quiet) == (0, {}, [])
        # Note: from here on I only check the first two output elements ([:2]) - the third is a list of 
        #   formatted output lines, MAYBE-TODO test that too?
        # single input file
        assert main([test_fa], **quiet)[:2] == (5, {0:2, 3:1, 5:1, 12:1})
        assert main([test_fa], total_seq_number_only=True, **quiet)[:2] == (5, {})
        assert main([test_fq], **quiet)[:2] == (7, {36:7})
        # two input files
        assert main([test_fq, test_fq], **quiet)[:2] == (14, {36:14})
        assert main([test_fa,test_fq], **quiet)[:2] == (12, {0:2, 3:1, 5:1, 12:1, 36:7})
        # input_collapsed_to_unique
        assert main([with_counts], input_collapsed_to_unique=False, **quiet)[:2] ==  (5, {0:2, 3:1, 5:1, 12:1})
        assert main([with_counts], input_collapsed_to_unique=True, **quiet)[:2] == (42, {0:11, 3:10, 5:20, 12:1})


if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(__doc__)
    # input format
    parser.add_option('-c','--input_collapsed_to_unique', action='store_true', default=False, 
                      help="Use to get correct total counts if the infile was collapsed to unique sequences using fastx_collapser, with original sequence counts encoded in the headers (a '>2-572' header means there were 572 identical sequences); default %default).")
    # TODO add options for explicit infile format specification (fasta or fastq) - some files have different extensions!
    # output format
    # MAYBE-TODO add -u option to keep track of only unique sequences?
    parser.add_option('-n','--total_seq_number_only', action='store_true', default=False, 
                      help="Output only the total seq number, without the seq length distribution; default %default).")
    parser.add_option('-z','--include_zeros', action='store_true', default=False, 
                      help="Include lengths with 0 counts in the output (output the whole range from lowest to highest length, instead of just the lengths present; default %default).")
    # -v and -q modify the same variable (verbosity) - default 1, -v makes it 2, -q makes it 0.
    parser.add_option('-v','--verbose', action='store_const', const=2, dest="verbosity", default=1, 
                      help="Print file type and seq count info about each processed file (default off).")
    parser.add_option('-q','--quiet', action='store_const', const=0, dest="verbosity", 
                      help="Don't print the header or the total number of counts (default off).")

    parser.add_option('-t','--test_functionality', action='store_true', default=False, 
                      help="Run the built-in unit test suite (ignores all other options/arguments; default %default).")

    (options, args) = parser.parse_args()

    # if run with -t option, do unit tests and quit
    if options.test_functionality:
        print("*** You used the -t option - ignoring all other options/arguments, running the built-in test suite. ***")
        # to run tests for current module, just run unittest.main, passing it only the filename 
        #   (by default it takes all of sys.argv and complains about options/arguments it can't recognize)
        unittest.main(argv=[sys.argv[0]])   # unittest.main automatically runs sys.exit()

    # otherwise parse arguments and run main function
    if not args:
        parser.print_help()
        sys.exit("\nError: At least one argument file required.")
    infiles = args

    main(infiles, options.total_seq_number_only, options.input_collapsed_to_unique, 
         options.include_zeros, options.verbosity)
