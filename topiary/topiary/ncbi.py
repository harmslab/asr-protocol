
__author__ = "Michael J. Harms"
__date__ = "2021-04-08"
__description__ = \
"""
functions interacting with ncbi databases and file types.
"""

from . import util

import pandas as pd
import numpy as np

# Modules for blasting, etc.
from Bio import SeqIO, Entrez
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
import Bio.Blast.Applications as apps
from Bio import pairwise2

Entrez.email = "DUMMY_EMAIL@DUMMY_URL.COM"

from tqdm.auto import tqdm

import re, sys, os, string, random, pickle, io, urllib, http
import warnings
import multiprocessing as mp

def read_blast_xml(filename):
    """
    Read BLAST XML format and return as a pandas data frame.

    filename: blast file name
    """

    # Read file.
    with open(filename, 'r') as f:
        blast_record = NCBIXML.read(f)

    # Prepare DataFrame fields.
    data = {'accession': [],
            'hit_def': [],
            'hit_id': [],
            'title': [],
            'length': [],
            'e_value': [],
            'sequence': [],
            'subject_start': [],
            'subject_end':[],
            'query_start':[],
            'query_end':[]}

    # Get alignments from blast result.
    for i, s in enumerate(blast_record.alignments):
        data['accession'].append(s.accession)
        data['hit_def'].append(s.hit_def)
        data['hit_id'].append(s.hit_id)
        data['title'].append(s.title)
        data['length'].append(s.length)
        data['e_value'].append(s.hsps[0].expect)
        data['sequence'].append(s.hsps[0].sbjct)
        data['subject_start'].append(s.hsps[0].sbjct_start)
        data['subject_end'].append(s.hsps[0].sbjct_end)
        data['query_start'].append(s.hsps[0].query_start)
        data['query_end'].append(s.hsps[0].query_end)

    # Port to DataFrame.
    return pd.DataFrame(data)

def local_blast(sequence,
                db,
                blast_program="blastp",
                keep_tmp=False,
                hitlist_size=100,
                e_value_cutoff=0.01,
                gapcosts=(11,1),
                **kwargs):
    """
    Perform a blast query using a sequence against a local blast data base.

    sequence: sequence as a string
    blast_program: blast program to use
    keep_tmp: whether or not to keep temporary blast output
    """

    recognized_functions = {"blastp":apps.NcbiblastpCommandline,
                            "blastn":apps.NcbiblastnCommandline,
                            "blastx":apps.NcbiblastxCommandline,
                            "tblastn":apps.NcbitblastnCommandline,
                            "tblastx":apps.NcbitblastxCommandline,
                            "psiblast":apps.NcbipsiblastCommandline,
                            "rpsblast":apps.NcbirpsblastCommandline,
                            "rpstblastn":apps.NcbirpstblastnCommandline,
                            "deltablast":apps.NcbideltablastCommandline}

    if not os.path.exists(f"{db}.psq"):
        err = f"{db} not found!\n"
        raise RuntimeError(err)

    # Make sure we can actually run the local blasting
    try:
        _local_blast = recognized_functions[blast_program]
    except KeyError:
        err = "\nblast_program '{}' not recognized.\n\nAllowed programs:\n".format(blast_program)
        for k in recognized_functions.keys():
            err += "    {}\n".format(k)
        raise ValueError(err)

    # Gap penalties
    gaps = '{} {}'.format(*gapcosts)

    # make a 10-character random string for temporary files
    tmp_file_root = "".join([random.choice(string.ascii_letters) for i in range(10)])
    input_file = "{}_blast_in.fasta".format(tmp_file_root)
    out_file = "{}_blast_out.xml".format(tmp_file_root)

    f = open(input_file,'w')
    f.write(f">{tmp_file_root}\n{sequence}\n")
    f.close()

    _local_blast(query=input_file,
                 cmd=blast_program,
                 db=db,
                 out=out_file,
                 outfmt=5,
                 max_target_seqs=hitlist_size,
                 threshold=e_value_cutoff,
                 gapopen=gapcosts[0],
                 gapextend=gapcosts[1],
                 **kwargs)()

    out_df = read_blast_xml(out_file)

    if not keep_tmp:
        os.remove(input_file)
        os.remove(out_file)

    return out_df

def parse_ncbi_line(line,accession=None,aliases=None):
    """
    Parse an ncbi line of the sort seen in the BLAST title field or on each
    line of a fasta file.

    accession: extract entry from line that matches acccession.  Ignores version
               (e.g. "1" in XXXXXXXXX.1).  If None, parse first entry on line.
    aliases: dictionary for standardizing protein names.  Key specifies what
             should be output, values degenerate names that map back to that
             key.  For example:
                 "S100A9":("S100-A9","S100 A9","S-100 A9")
             would replace "S100-A9", "S100 A9", and "S-100 A9" with "S100A9"

    Returns a dictionary with following keys:
        raw_line -> unprocessed line (input)
        line -> processed line (remove multiple titles)
        protein -> protein name (renamed using aliases)
        structure -> whether or not this is a structure (bool)
        low_quality -> whether or not this is low quality (bool)
        predicted -> whether or not this is predicted (bool)
        precursor -> whether or not this is a precursor (bool)
        isoform -> whether or not this is an isoform (bool)
        hypothetical -> whether or not this is hypothetical (bool)
        partial -> whether or not this is a partial sequence (bool)
    """

    out = {"raw_line":line}

    # Clean up leading and trailing spaces
    line = line.strip()

    # Clean up leading ">" if it's there (difference between fasta and blast
    # lines)
    if line.startswith(">"):
        line = line[1:]

    # Particularly in BLAST hits, there are often multiple entries separated
    # by ">" (basically, pasted together fasta entries). Only take one.  If
    # accession is specified, take that accesion.  If no accession is
    # specified, take the first one.

    # Split on ">"
    line_splitter = re.compile(">")
    entries_on_line = [s.strip() for s in line_splitter.split(line)]

    # Extract accession from entry (assumes xxx|acccession stuff).  Ignores
    # trailing accession version number XXXXXXXX.1 -> XXXXXXXX
    accession_dict = {}
    for e in entries_on_line:
        try:
            if e.startswith("pdb"):
                k = [v.split()[0] for v in e.split("|")[1:3]]
                k = "{}_{}".format(*k)
            else:
                k = e.split("|")[-2].split(".")[0].strip()
            accession_dict[k] = e
        except IndexError:
            warnings.warn(f"Could not parse line {e}. Skipping.\n")
            return None

    # If accession is specified, grab entry corresponding to that from the
    # line.  If not, grab first entry.
    if not accession:
        accession = list(accession_dict.keys())[0]

    accession = accession.split(".")[0].strip()
    try:
        line = accession_dict[accession]
    except KeyError:
        warn = f"accession '{accession}' not found in line '{line}'\n"
        warn += f"Using accession {accession} anyway."
        warnings.warn(warn)

    out["accession"] = accession

    # Record the processed line
    out["line"] = line

    # Grab meta stuff (predicted, isoform, etc.)
    meta = util.grab_line_meta_data(line)
    for m in meta:
        out[m] = meta[m]

    # Look for species name (thing within [ xxx ])
    species = None
    species_pattern = re.compile("\[.*?]")

    sm = species_pattern.search(line)
    if sm:
        out["species"] = sm.group(0)[1:-1]
    else:
        out["species"] = None

    # Use "aliases" to clean up the protein bit
    if aliases is None:
        aliases = {}
    try:
        aliases[""]
    except KeyError:
        aliases[""] = ("protein",
                       "product",
                       "PREDICTED:",
                       "Crystal structure of",
                       "hypothetical")
    for a in aliases:
        for p in aliases[a]:
            ap = re.compile(p,re.IGNORECASE)
            line = ap.sub(a,line)

    # Clean up any double spaces introduced into the line at this point
    line = re.sub("  "," ",line)

    # Protein name (takes between '| XXXXXXXXX [' ).
    out["protein"] = line.split("|")[-1].split("[")[0].strip()

    return out

def _entrez_download_thread(args):
    """
    Download ids from the NCBI. args is interpreted as:

    index: number of request, used for sorting results after multithreading
    ids: string with comma-separated list of ids to download
    num_tries_allowed: how many attempts should be made to actually do
                       a given query
    queue: multiprocessing queue in which to store results
    """

    index = args[0]
    ids = args[1]
    num_tries_allowed = args[2]
    queue = args[3]

    # While we haven't run out of tries...
    tries = 0
    while tries < num_tries_allowed:

        # Try to read output.
        try:
            handle = Entrez.efetch(db="protein",
                                   id=ids,
                                   rettype="fasta",
                                   retmode="text")
            out = handle.read()
            handle.close()

        # If some kind of http error or timeout, set out to None
        except (urllib.error.HTTPError,http.client.IncompleteRead):
            out = None

        # If out is None, try again. If not, break out of loop--success!
        if out is None:
            tries += 1
            continue
        else:
            break

    # We failed to download. Throw an error.
    if out is None:
        err = "could not download sequences\n"
        raise ValueError(err)

    # Record out and initial index in the output queue.
    queue.put((index,out))

def entrez_download(to_download,block_size=50,num_tries_allowed=5,num_threads=-1):
    """
    Download sequences off of entrez, catching errors. This is done in a
    multi-threaded fashion, allowing multiple requests to NCBI at the same
    time.

    to_download: list of ids to download
    block_size: download in chunks this size
    num_tries_allowed: number of times to try before giving up and throwing
                       an error.
    num_threads: number of threads to use. if -1, use all available.
    """

    # Figure out number of threads to use
    if num_threads < 0:
        try:
            num_threads = mp.cpu_count()
        except NotImplementedError:
            num_threads = os.cpu_count()
            if num_threads is None:
                warnings.warn("Could not determine number of cpus. Using single thread.\n")
                num_threads = 1


    # Figure out how many blocks we're going to download
    num_blocks = len(to_download) // block_size
    if len(to_download) % block_size > 0:
        num_blocks += 1
    print(f"Downloading {num_blocks} blocks of <={block_size} sequences... ")

    # queue will hold results from each download batch.
    queue = mp.Manager().Queue()
    with mp.Pool(num_threads) as pool:

        # This is a bit obscure. Build a list of args to pass to the pool. Each
        # tuple of args matches the args in _entrez_download_thread.
        # all_args has all len(df) reverse blast runs we want to do.
        all_args = []
        for i in range(0,len(to_download),block_size):
            ids = ",".join(to_download[i:(i+block_size)])
            all_args.append((i,ids,num_tries_allowed,queue))

        # Black magic. pool.imap() runs a function on elements in iterable,
        # filling threads as each job finishes. (Calls _entrez_download_thread
        # on every args tuple in all_args). tqdm gives us a status bar.
        # By wrapping pool.imap iterator in tqdm, we get a status bar that
        # updates as each thread finishes.
        list(tqdm(pool.imap(_entrez_download_thread,all_args),
                  total=len(all_args)))

    # Get results out of the queue. Sort by i to put in correct order
    results = []
    while not queue.empty():
        results.append(queue.get())
    results.sort()

    # Get final downloaded stuff
    total_download = []
    for r in results:
        total_download.append(r[1])

    print("Done.")

    return "".join(total_download)
