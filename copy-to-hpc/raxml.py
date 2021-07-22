#!/usr/bin/env python3
__description__ = \
"""
Wrap key features of raxml for asr work.
"""
import pastml.acr
from ete3 import Tree

import pandas as pd
import numpy as np
import subprocess, os, glob, re, sys, time, random, argparse, inspect, string
import shutil, warnings

RAXML_BINARY = "raxmlHPC-PTHREADS-SSE3"

# List of AAs as ordered in RAxML output
AA_LIST = ['ALA','ARG','ASN','ASP','CYS',
           'GLN','GLU','GLY','HIS','ILE',
           'LEU','LYS','MET','PHE','PRO',
           'SER','THR','TRP','TYR','VAL']

# Dictionary for converting three letter to single letter aa code
THREE_TO_ONE = dict([('ALA','A'),
                     ('CYS','C'),
                     ('ASP','D'),
                     ('GLU','E'),
                     ('PHE','F'),
                     ('GLY','G'),
                     ('HIS','H'),
                     ('HSE','H'),
                     ('HSD','H'),
                     ('ILE','I'),
                     ('LYS','K'),
                     ('LEU','L'),
                     ('MET','M'),
                     ('MSE','M'),
                     ('ASN','N'),
                     ('PRO','P'),
                     ('GLN','Q'),
                     ('ARG','R'),
                     ('SER','S'),
                     ('THR','T'),
                     ('VAL','V'),
                     ('TRP','W'),
                     ('TYR','Y')])

# lists of chemically similar amino acids. the same amino acid can occur
# in different lists.
CHEM_SIMILAR = [["A","C","N","Q","S","T"],
                ["A","V"],
                ["H","K","R"],
                ["D","E"],
                ["F","W","Y"],
                ["I","L","M","V"]]


def _gen_seed():
    """
    Generate a random string of 10 integers and return as a string for passing
    to raxml.
    """

    return "".join([f"{random.choice(range(10)):d}" for _ in range(10)])

def _create_new_dir(dir_name=None,dir_base=RAXML_BINARY):
    """
    Create a new directory.

    dir_name: if specified, name the directory this (ignoring dir_base)
    dir_base: if dir_name is not specified, name the directory dir_base_random_bit

    returns name of created directory
    """

    if dir_name is None:

        rand_name = "".join([random.choice(string.ascii_letters)
                              for _ in range(10)])
        dir_name = f"{dir_base}_{rand_name}"

    if os.path.exists(dir_name):
        err = f"{dir_name} already exists.\n"
        raise FileExistsError(err)

    os.mkdir(dir_name)

    return dir_name


def _parse_raxml_info_for_aic(info_file):
    """
    Open an info file from a tree generation run and get parameters important
    for calculating an AIC value in model testing.
    """

    with open(info_file,'r') as f:
        for line in f:
            if re.search("likelihood:",line):
                L = float(line.strip().split()[-1])

            if re.search("AIC-TEST\(BR-LEN\)",line):
                N = int(line.strip().split()[-1])

    return {"L":L,"N":N,"AIC":(2*N - 2*L)}


def _fix_raxml_tree(raxml_tree,out_file):
    """
    Clean up an raxml newick tree so it is readable by other software.

    raxml_tree: newick file dumped by raxml
    out_file: name of file to write out. (does not check for existance; will
              overwrite)
    """

    # Open raxml tree
    f = open(raxml_tree,"r")
    tree = f.read()
    f.close()

    # Deal with wacky support patterns in raxml output
    support_pattern = re.compile("\):.*?\[.*?\]")
    specific_matches = []
    for x in support_pattern.finditer(tree):
        m = x.group(0)
        support = m.split("[")[1][:-1]
        length = m.split(":")[1].split("[")[0]
        out = f"){support}:{length}"

        p = re.sub("\)","\\\)",m)
        p = re.sub("\[","\\\[",p)
        p = re.sub("\]","\\\]",p)

        specific_matches.append((re.compile(p),out))

    # Actually do substitutions
    for s in specific_matches:
        tree = s[0].sub(s[1],tree)

    # Write output file
    g = open(out_file,"w")
    g.write(tree)
    g.close()

def _copy_root(unrooted_newick,rooted_newick,output_newick,unrooted_tree_fmt=0):
    """
    Root the tree in an unrooted newick file using the root from a rooted
    newick file with the same taxa. Write to an output file.

    unrooted_newick: newick file containing an unrooted tree
    rooted_newick: newick file containing a rooted tree with the same taxa
                   as the unrooted tree
    output_newick: output file to write results
    unrooted_tree_fmt: what to preserve from unrooted tree. integer.
                       interpretation is done by ETE3 (table below
                       current as of v. 3.1.1).

     |        ======  ==============================================
     |        FORMAT  DESCRIPTION
     |        ======  ==============================================
     |        0        flexible with support values
     |        1        flexible with internal node names
     |        2        all branches + leaf names + internal supports
     |        3        all branches + all names
     |        4        leaf branches + leaf names
     |        5        internal and leaf branches + leaf names
     |        6        internal branches + leaf names
     |        7        leaf branches + all names
     |        8        all names
     |        9        leaf names
     |        100      topology only
     |        ======  ==============================================

    """

    # Load trees
    rooted_tree = Tree(rooted_newick)
    unrooted_tree = Tree(unrooted_newick,format=unrooted_tree_fmt)

    # Make sure they have the same taxa
    rooted_leaves = set(rooted_tree.get_leaf_names())
    unrooted_leaves = set(unrooted_tree.get_leaf_names())
    if rooted_leaves != unrooted_leaves:
        err = "both trees must have the exact same leaves\n"
        raise ValueError(err)

    left_leaves = []
    right_leaves = []

    root_node = None
    left_anc_node = None

    # Pre-order traverses root, left, right
    for node in rooted_tree.traverse("preorder"):

        # First iteration gets root node
        if root_node is None:
            root_node = node
            continue

        # Second iteration gets ancestor of all left. If this is a leaf,
        # there is a single left outgroup.
        if left_anc_node is None:
            left_anc_node = node
            if left_anc_node.is_leaf():
                left_leaves.append(left_anc_node.name)
            continue

        # Look for ancestor of all left
        if left_anc_node in node.get_ancestors():
            if node.is_leaf():
                left_leaves.append(node.name)

        # Or all right
        else:
            if node.is_leaf():
                right_leaves.append(node.name)

    # If we have single outgroups on the left or right, root on that
    if len(left_leaves) == 1:
        unrooted_tree.set_outgroup(left_leaves[0])
    elif len(right_leaves) == 1:
        unrooted_tree.set_outgroup(right_leaves[0])

    # Otherwise, try to root on last common ancestor of left or right. This
    # may throw error, depending on tree topology, so we try left first,
    # and then try right if that does not work.
    else:

        try:
            root_anc = unrooted_tree.get_common_ancestor(*left_leaves)
            unrooted_tree.set_outgroup(root_anc)
        except ete3.coretype.tree.TreeError:
            try:
                root_anc = unrooted_tree.get_common_ancestor(*right_leaves)
                unrooted_tree.set_outgroup(root_anc)
            except ete3.coretype.tree.TreeError:
                err = "could not root tree for asr!\n"
                raise ValueError(err)

    # Write out newly rooted tree
    unrooted_tree.write(outfile=output_newick)

def _get_ancestral_gaps(alignment_file,tree_file):
    """
    Get ancestral gaps from an alignment and raxml output tree file. Gaps are
    reconstructed by parsimony using the DOWNPASS algorithm as implemented in
    pastml.

    alignment_file: phy file used to generate ancestors in RAxML
    tree_file: output tree file with labeled internal nodes

    output: dictionary keying internal node names to lists of True (gap),
            False (no gap), and None (gapping unclear) for each site in that
            ancestor.
    """

    # Read the alignment file
    counter = 0
    leaf_names = []
    with open(alignment_file) as f:
        for line in f:

            # First line
            if counter == 0:
                col = line.split()
                num_taxa = int(col[0])
                num_sites = int(col[1])
                counter += 1

                char_matrix = np.zeros((num_taxa,num_sites),dtype=np.uint8)
                continue

            # Next, blank line
            if line.strip() == "":
                counter += 1
                continue

            # Alternating leaf id and sequence lines
            if counter % 2 == 0:
                leaf_names.append(line.strip())
                counter += 1
            else:
                index = (counter - 3)//2
                char_matrix[index,:] = np.array([c == "-" for c in line.strip()])
                counter += 1

    # Create a data frame where indexes are leaf names and columns are each gap
    out_dict = {}
    for i in range(char_matrix.shape[1]):
        out_dict[f"g{i}"] = char_matrix[:,i]
    gap_df = pd.DataFrame(out_dict)
    gap_df.index = leaf_names

    # Gaps, named as column names
    gap_names = list(gap_df.columns)

    # Load the tree, keeping the internal node names
    tree = Tree(tree_file,format=1)

    # Reconstruct gaps across tree by parsimony
    acr_result = pastml.acr.acr(tree,gap_df,prediction_method="DOWNPASS")

    # Create dictionary keying ancestor name to gap status across sequence
    gap_anc_dict = {}
    for node in tree.traverse('preorder'):
        if node.name in leaf_names:
            continue

        gap_anc_dict[node.name] = []

        # Go through gaps and grab from node feature
        for g in gap_names:
            state = node.__dict__[g]
            if len(state) == 1:
                if state == {0}:
                    state = False
                else:
                    state = True
            else:
                state = None

            gap_anc_dict[node.name].append(state)

    return gap_anc_dict

def _parse_raxml_anc_output(anc_file,phy_file,tree_file,alt_cutoff=0.25):
    """
    Parse raxml marginal ancestral state reconstruction output and put out in
    human-readable fashion.

    anc_file: ancestor file as written out by raxml
    phy_file: phylip alignment used to create ancestors
    tree_file: *output* tree file written by raxml
    alt_cutoff: cutoff (inclusive) for identifying plausible alternate states

    writes fasta file and csv file with ancestors

    returns pandas data frame with ancestral information
    """

    # Get gaps, reconstructed by parsimony
    gap_anc_dict = _get_ancestral_gaps(phy_file,tree_file)

    # Read pp file into a list of ancestors (anc_list) and posterior
    # probabilities for amino acids at each site
    anc_list = []
    anc_all_pp = []
    last_line = ""
    with open(anc_file) as f:
        for l in f:
            line = l.strip()

            # New ancestor if last line was blank
            if last_line == "":
                anc_list.append(line)
                anc_all_pp.append([])
            else:

                # If not a blank line, parse ancestor posterior probabilities
                if line != "":
                    pp = np.array([float(c) for c in line.split()])
                    anc_all_pp[-1].append(pp)

            last_line = line

    # Create data structures for output
    out = []
    for_df = {"anc":[],"site":[],"gap":[],
              "ml_state":[],"ml_pp":[],
              "alt_state":[],"alt_pp":[],
              "site_type":[],"entropy":[]}

    # Go through each ancestor
    for i in range(len(anc_list)):

        anc = anc_list[i]
        gap_list = gap_anc_dict[anc]

        # Get ml and next best seq and posterior probability. Entropy holds
        # the shannon entropy for all posterior probabilities across the
        # row.  If entropy is very high, we have very little confidence --
        # likely a gap.
        anc_ml_seq = []
        anc_ml_pp = []
        anc_alt_seq = []
        anc_alt_pp = []
        entropy = []

        # Go through sites
        for j in range(len(anc_all_pp[i])):

            # Get posterior probability vector for site
            pp = anc_all_pp[i][j].copy()

            # Get entropy, ignoring positions with pp = 0
            tmp_pp = pp[pp != 0.0]
            entropy.append(np.sum(-tmp_pp*np.log(tmp_pp)))

            # Get ml posterior probability and sequence
            anc_ml_pp.append(np.max(pp))
            anc_ml_seq.append(THREE_TO_ONE[AA_LIST[np.argmax(pp)]])

            # Set max to zero to get next best
            pp[np.argmax(pp)] = 0.0

            # Get second best posterior probability and sequence
            anc_alt_pp.append(np.max(pp))
            anc_alt_seq.append(THREE_TO_ONE[AA_LIST[np.argmax(pp)]])

            # Add gap!
            if gap_anc_dict[anc][j] == True:
                anc_ml_pp[-1] = 1.0
                anc_ml_seq[-1] = "-"
                anc_alt_pp[-1] = 1.0
                anc_alt_seq[-1] = "-"

        # Convert lists read into numpy arrays
        ml_pp = np.array(anc_ml_pp)
        ml_seq = np.array(anc_ml_seq)
        alt_pp = np.array(anc_alt_pp)
        alt_seq = np.array(anc_alt_seq)

        # Create alt all sequence. This sequence is ML at all states except
        # those with pp >= some alternate cutoff (usually 0.25)
        alt_mask = alt_pp >= alt_cutoff
        alt_all_seq = ml_seq.copy()
        alt_all_seq[alt_mask] = alt_seq[alt_mask]
        alt_all_pp = ml_pp.copy()
        alt_all_pp[alt_mask] = alt_pp[alt_mask]

        # Write out information for data frame
        for_df["anc"].extend([i for _ in range(len(anc_all_pp[i]))])
        for_df["site"].extend([_ for _ in range(len(anc_all_pp[i]))])
        for_df["gap"].extend(gap_anc_dict[anc])
        for_df["ml_state"].extend(list(ml_seq))
        for_df["ml_pp"].extend(list(ml_pp))
        for_df["alt_state"].extend(list(alt_seq))
        for_df["alt_pp"].extend(list(alt_pp))
        for_df["entropy"].extend(list(entropy))

        # Classify sites according to semi-intelligent criteria. If gapping
        # unclear by parsimony, call as "possible gap."  If it is a gap, call
        # as "gap." If second best pp is above cutoff, call
        # ambiguous.  Look at ml and alt aa to see if this is a chemically
        # similar or dissimilar ambiguity.  Otherwise, call 'good'
        site_type = []
        for j in range(len(anc_all_pp[i])):

            if gap_anc_dict[anc][j] is None:
                site_type.append("possible gap")
                warnings.warn(f"ancestor {anc} has unresolved gap(s)\n")

            elif gap_anc_dict[anc][j]:
                site_type.append("gap")

            elif alt_pp[j] >= alt_cutoff:

                found_match = False
                for ml_set in CHEM_SIMILAR:
                    if ml_seq[j] in ml_set:
                        if alt_seq[j] in ml_set:
                            found_match = True
                            continue
                if found_match:
                    site_type.append("ambig_similar")
                else:
                    site_type.append("ambig_dissimilar")
            else:
                site_type.append("good")

        for_df["site_type"].extend(site_type)

        # Calculate average pp, log sum pp, and number of ambiguous sites for ML sequence
        avgPP = np.mean(ml_pp)
        lnPP = np.sum(np.log(ml_pp))
        num_ambig = np.sum(alt_mask)

        # Write ML sequence to fasta
        out.append(f">{anc}|avgPP:{avgPP:.3f}|lnPP:{lnPP:.3f}|num_ambig:{num_ambig}\n")
        out.append("".join(ml_seq))
        out.append("\n")

        # Calculate average pp, log sum pp, and number of ambiguous sites for alt sequence
        avgPP = np.mean(alt_all_pp)
        lnPP = np.sum(np.log(alt_all_pp))
        num_ambig = np.sum(alt_mask)

        # Write alt sequence to fasta
        out.append(f">{anc}_altAll|avgPP:{avgPP:.3f}|lnPP:{lnPP:.3f}|num_ambig:{num_ambig}\n")
        out.append("".join(alt_all_seq))
        out.append("\n")

    # Write final fasta file
    f = open("ancestors.fasta","w")
    f.write("".join(out))
    f.close()

    # Write ancestor df to csv
    df = pd.DataFrame(for_df)
    df.to_csv("ancestors.csv")

    # Return df
    return df

def _copy_input_file(input_file,dir_name,make_input_dir=True):
    """
    copy an input file into a directory in a stereotyped way.

    If make_input_dir is specified, copy input_file into dir_name/00_input,
    creating 00_input if necessary.  If make_input_dir is not specified,
    copy in the file as dir_name/input_{input_file}.

    input_file: file to copy in
    dir_name: copy into dir_name
    make_input_dir: (bool) make input directory 00_input or not.
    """

    file_alone = os.path.split(input_file)[-1]

    # If we are putting this into an input subdirectory
    if make_input_dir:
        input_dir = os.path.join(dir_name,"00_input")
        if not os.path.exists(input_dir):
            os.mkdir(input_dir)
        file_alone = os.path.join("00_input",file_alone)

    # If we are not making an input directory, append input_ to front
    else:
        if not file_alone.startswith("input"):
            file_alone = f"input_{file_alone}"

    shutil.copy(input_file,os.path.join(dir_name,file_alone))

    return file_alone

def _run_raxml(algorithm=None,
               alignment_file=None,
               tree_file=None,
               model=None,
               name=None,
               seed=True,
               threads=1,
               raxml_binary=RAXML_BINARY,
               other_args=[]):
    """
    Run raxml. Creates a working directory, copies in the relevant files, runs
    there, and then returns to the previous directory.

    algorithm: letter passed to the -f flag
    alignment_file: alignment file in .phy format
    tree_file: tree file in .newick format
    model: model in format recognized by -m
    name: name to pass via -n. If specified, this will also be the name of the
          working directory.
    seed: true/false. if true, pass a randomly generated seed to raxml.
    threads: number of threads to use
    raxml_binary: raxml binary to use
    other_args: list of arguments to pass to raxml
    """

    dir_name = _create_new_dir(dir_name=name,dir_base=raxml_binary)

    if alignment_file is not None:
        alignment_file = _copy_input_file(alignment_file,
                                          dir_name,
                                          make_input_dir=True)
    if tree_file is not None:
        tree_file = _copy_input_file(tree_file,
                                     dir_name,
                                     make_input_dir=True)

    os.chdir(dir_name)

    cmd = [raxml_binary]

    if algorithm is not None:
        cmd.extend(["-f",algorithm])

    if alignment_file is not None:
        cmd.extend(["-s",alignment_file])

    if tree_file is not None:
        cmd.extend(["-t",tree_file])

    if model is not None:
        cmd.extend(["-m",model])

    if name is not None:
        cmd.extend(["-n",name])

    if seed:
        cmd.extend(["-p",_gen_seed()])

    cmd.extend(["-T",f"{threads:d}"])

    for a in other_args:
        cmd.append(a)

    full_cmd = " ".join(cmd)
    print(f"Running '{full_cmd}'")
    sys.stdout.flush()

    ret = subprocess.run(cmd,stdout=subprocess.PIPE)
    if ret.returncode != 0:
        err = f"ERROR: raxml returned {ret.returncode}\n\n"
        err += "------------------------------------------------------------\n"
        err += " raxml output \n"
        err += "------------------------------------------------------------\n"
        err += "\n\n"

        err += "".join([line for line in ret.stdout.decode()])

        raise RuntimeError(err)

    os.chdir("../")

def _generate_parsimony_tree(alignment_file,
                             name="parsimony-tree",
                             threads=1,
                             raxml_binary=RAXML_BINARY):
    """
    Generate a parsimony tree from an alignment.

    alignment_file: alignment file in .phy format
    name: name to give directory and pass via -n to raxml.
    threads: number of threads to use
    raxml_binary: raxml binary to use
    """

    _run_raxml(alignment_file=alignment_file,
               name=name,
               seed=True,
               model="PROTGAMMAJTT",
               threads=threads,
               raxml_binary=raxml_binary,
               other_args=["-y"])

def find_best_model(alignment_file,
                    tree_file=None,
                    model_matrices=["DAYHOFF","DCMUT","JTT","MTREV","WAG",
                                    "RTREV","CPREV","VT","BLOSUM62","MTMAM",
                                    "LG","MTART","MTZOA","PMB","HIVB","HIVW",
                                    "JTTDCMUT","FLU","STMTREV","LG4M","LG4X",
                                    "GTR"],
                    model_rates=["CAT","GAMMA"],
                    model_freqs=["","F","X"],
                    output=None,
                    threads=1,
                    raxml_binary=RAXML_BINARY):
    """
    Find the best phylogentic model to use for tree and ancestor reconstruction
    given an alignment and (possibly) a tree.

    alignment_file: alignment file in .phy format
    tree_file: tree file in newick format. If not specified, parsimony tree
               is generated and used
    model_matrices: list of model matrices to check
    model_rates: ways to treat model rates
    model_freqs: ways to treat model freqs.
    output: directory for output. it none, will generate random name
    threads: number of threads to use
    raxml_binary: raxml binary to use
    """

    # Models that might concievably, but do not actually, exist
    NOT_ALLOWED_MODELS = ["PROTCATGTRF","PROTGAMMAGTRF"]

    # Make sure alignment file exists
    if not os.path.exists(alignment_file):
        err = f"alignment file {alignment_file} does not exist\n"
        raise ValueError(err)

    # Create output directory
    output = _create_new_dir(dir_name=output,dir_base="find_best_model")

    # Copy files into input directory
    alignment_file = _copy_input_file(alignment_file,
                                      output,
                                      make_input_dir=True)
    if tree_file is not None:
        tree_file = _copy_input_file(tree_file,
                                     output,
                                     make_input_dir=True)

    # Move into the output directory
    os.chdir(output)

    # Generate a parsimony tree if not was specified
    if tree_file is None:
        _generate_parsimony_tree(alignment_file,name="01_make-parsimony-tree",
                                 threads=threads,raxml_binary=raxml_binary)
        tree_file = "02_parsimony-tree.newick"
        shutil.copy(os.path.join("01_make-parsimony-tree",
                                 "RAxML_parsimonyTree.01_make-parsimony-tree"),
                    tree_file)

    # Dictionary to hold stats for each model
    out = {"model":[],
           "L":[],
           "N":[],
           "AIC":[]}

    # Go over all combos of the requested matrices, rates, and freqs.
    for matrix in model_matrices:
        for rate in model_rates:
            for freq in model_freqs:

                # Check for incompatible matrix/freq/rate combos
                if matrix in ["LG4M","LG4X"]:
                    if freq != "":
                        continue
                    if rate == "CAT":
                        continue

                # Construct model string
                model = f"PROT{rate}{matrix}{freq}"

                # Skip models that are not actually implemented in raxml
                if model in NOT_ALLOWED_MODELS:
                    continue

                # Optimize branch lengths etc. on the existing tree
                _run_raxml(algorithm="e",
                           alignment_file=alignment_file,
                           tree_file=tree_file,
                           model=model,
                           name="tmp",
                           threads=threads,
                           raxml_binary=raxml_binary)

                # Grab the info file from this run
                os.chdir("tmp")
                info_out = None
                info_files = glob.glob("*_info*")
                if len(info_files) != 1:
                    err = "could not find info file in output"
                    raise RuntimeError(err)
                info_out = info_files[0]

                # Get results from the info file
                result = _parse_raxml_info_for_aic(info_out)
                out["model"].append(model)
                out["L"].append(result["L"])
                out["N"].append(result["N"])
                out["AIC"].append(result["AIC"])

                # Get out of temporary directory and nuke
                os.chdir("..")
                shutil.rmtree("tmp")

    # Create a csv file sorted best to worst aic
    df = pd.DataFrame(out)
    min_aic = np.min(df.AIC)
    df["p"] = np.exp((min_aic - df.AIC)/2)
    indexer = np.argsort(df.p)[::-1]
    df = df.iloc[indexer,:]
    df.to_csv("model-comparison.csv")

    # Get best model
    best_model = df.model.iloc[0]

    # Write model to a file
    f = open("best-model.txt","w")
    f.write(best_model)
    f.close()

    # Print best model to stdout
    print(f"\n\nBest model: {best_model}\nAIC Prob:{df.p.iloc[0]}\n\n")

    os.chdir("../")

def generate_ml_tree(alignment_file,
                     model,
                     tree_file=None,
                     output=None,
                     threads=1,
                     raxml_binary=RAXML_BINARY):
    """
    Generate maximum likelihood tree with SH supports from an alignment given
    a substitution model.

    alignment_file: alignment in .phy format
    model: model (e.g. PROTGAMMAJTTF).  spit out in best-model.txt from
           "find_best_model" function/mode
    tree_file: tree_file in newick format. If not specified, a parsimony tree
               will be generated. used as starting point.
    output: name out output directory.
    threads: number of threads to use
    raxml_binary: what raxml binary to use
    """

    # Make sure alignment file exists
    if not os.path.exists(alignment_file):
        err = f"alignment file {alignment_file} does not exist\n"
        raise ValueError(err)

    # Create output directory
    output = _create_new_dir(dir_name=output,dir_base="ml_tree")

    # Copy files into input directory
    alignment_file = _copy_input_file(alignment_file,
                                      output,
                                      make_input_dir=True)
    if tree_file is not None:
        tree_file = _copy_input_file(tree_file,
                                     output,
                                     make_input_dir=True)
    # Move into directory
    os.chdir(output)

    # Generate a parsimony tree if not was specified
    if tree_file is None:
        _generate_parsimony_tree(alignment_file,name="01_make-parsimony-tree",
                                 threads=threads,raxml_binary=raxml_binary)
        tree_file = "02_parsimony-tree.newick"
        shutil.copy("01_make-parsimony-tree/RAxML_parsimonyTree.01_make-parsimony-tree",
                    tree_file)

    # Run raxml to create tree
    _run_raxml(alignment_file=alignment_file,
               tree_file=tree_file,
               model=model,
               name="03_make-ml-tree",
               threads=threads,
               raxml_binary=raxml_binary)
    tree_file = "04_ml-tree.newick"
    shutil.copy("03_make-ml-tree/RAxML_result.03_make-ml-tree",tree_file)

    # Get SH supports for this tree
    _run_raxml(algorithm="J",
               alignment_file=alignment_file,
               tree_file=tree_file,
               model=model,
               name="05_add-supports",
               seed=True,
               threads=threads,
               raxml_binary=raxml_binary)

    tree_file = "06_tree-with-supports.newick"
    shutil.copy("05_add-supports/RAxML_fastTreeSH_Support.05_add-supports",
                tree_file)

    # Clean up supports and make them figtree readable.
    _fix_raxml_tree(tree_file,"07_final-tree.newick")

    os.chdir("../")


def generate_ancestors(alignment_file,
                       tree_file,
                       model,
                       output=None,
                       threads=1,
                       raxml_binary=RAXML_BINARY,
                       root_tree=False,
                       alt_cutoff=0.25):

    # Make sure alignment file exists
    if not os.path.exists(alignment_file):
        err = f"alignment file {alignment_file} does not exist\n"
        raise ValueError(err)

    # Make sure tree file exists
    if not os.path.exists(tree_file):
        err = f"tree file {tree_file} does not exist\n"
        raise ValueError(err)

    # Create output directory
    output = _create_new_dir(dir_name=output,dir_base="marginal_anc")

    # Copy files into input directory
    alignment_file = _copy_input_file(alignment_file,
                                      output,
                                      make_input_dir=True)
    tree_file = _copy_input_file(tree_file,
                                 output,
                                 make_input_dir=True)

    os.chdir(output)

    if root_tree:
        _run_raxml(algorithm="I",
                   alignment_file=alignment_file,
                   tree_file=tree_file,
                   model=model,
                   seed=True,
                   name="01_root-tree",
                   threads=threads,
                   raxml_binary=raxml_binary)
        tree_file = "01_root-tree/RAxML_rootedTree.01_root-tree"

    # Record this file name as one with roots.
    rooted_tree_file = tree_file

    # Optimize branch lengths for this tree
    _run_raxml(algorithm="e",
               alignment_file=alignment_file,
               tree_file=tree_file,
               model=model,
               name="02_optimize-branch-lengths",
               threads=threads,
               raxml_binary=raxml_binary)

    tree_file_with_bl = "02_optimize-branch-lengths/RAxML_result.02_optimize-branch-lengths"

    # Get SH supports for this tree
    _run_raxml(algorithm="J",
               alignment_file=alignment_file,
               tree_file=tree_file,
               model=model,
               name="03_add-supports",
               seed=True,
               threads=threads,
               raxml_binary=raxml_binary)

    tree_file_with_supports = "03_add-supports/RAxML_fastTreeSH_Support.03_add-supports"

    tree_file_for_asr = "04_tree-for-asr.newick"
    _copy_root(tree_file_with_bl,
               rooted_tree_file,
               tree_file_for_asr,
               unrooted_tree_fmt=0)

    _run_raxml(algorithm="A",
               alignment_file=alignment_file,
               tree_file=tree_file_for_asr,
               model=model,
               seed=True,
               name="05_marginal-anc",
               threads=threads,
               raxml_binary=raxml_binary)

    _parse_raxml_anc_output("05_marginal-anc/RAxML_marginalAncestralProbabilities.05_marginal-anc",
                            alignment_file,
                            "05_marginal-anc/RAxML_nodeLabelledRootedTree.05_marginal-anc",
                            alt_cutoff=alt_cutoff)

    os.chdir("../")


def main(argv):

    parser = argparse.ArgumentParser(prog="raxml.py",description='Run RAxML.')

    parser.add_argument("mode",nargs=1,help="mode. must be ml, model, or anc")

    # Arguments shared by all modes
    parser.add_argument("-a","--alignment",nargs=1,help="alignment file (.phy)")
    parser.add_argument("-t","--tree",nargs=1,help="tree file (newick)")
    parser.add_argument("-m","--model",nargs=1,help="phylogenetic model")
    parser.add_argument("-o","--output",nargs=1,help="name for outputs")
    parser.add_argument("-T","--threads",nargs=1,help="number of threads",type=int,
                        default=1)
    parser.add_argument("-b","--binary",nargs=1,help="raxml binary",
                        default=RAXML_BINARY)


    parser.add_argument("--anc-alt-pp",nargs=1,
                        help="anc mode: alternate ancestor posterior probability cutoff",
                        type=float,
                        default=0.25)

    parser.add_argument("--anc-root-tree",
                        help="anc mode: should raxml root tree before ancestor reconstruction?",
                        action='store_true')

    args = parser.parse_args(argv)

    modes = {"model":find_best_model,
             "ml":generate_ml_tree,
             "anc":generate_ancestors}

    # Figure out what function to run
    mode = args.mode[0]

    try:
        func = modes[mode]
    except KeyError:
        err = f"\n\nmode {mode} not recognized. should be one of:\n"
        for m in modes.keys():
            err += f"    {m}\n"
        err += "\n\n"
        raise ValueError(err)

    # Get kwarg arguments specified
    kwargs = {}
    tree_file = args.tree
    if tree_file is not None:
        kwargs["tree_file"] = tree_file[0]

    alignment_file = args.alignment
    if alignment_file is not None:
        kwargs["alignment_file"] = alignment_file[0]

    model = args.model
    if model is not None:
        kwargs["model"] = model[0]

    output = args.output
    if output is not None:
        kwargs["output"] = output[0]

    threads = args.threads
    if threads is not None:
        if type(threads) is int:
            kwargs["threads"] = threads
        else:
            kwargs["threads"] = threads[0]

    binary = args.binary
    if binary is None:
        binary = RAXML_BINARY
    kwargs["raxml_binary"] = binary

    # Ancestor-specific arguments
    if mode == "anc":
        kwargs["alt_cutoff"] = args.anc_alt_pp
        kwargs["root_tree"] = args.anc_root_tree

    # Figure out what kwargs are required and allowed for the function being
    # called
    required = []
    param = inspect.signature(func).parameters
    for p in param:
        if param[p].default is param[p].empty:
            required.append(p)
    allowed = list(param)

    extra = []
    for k in kwargs:
        try:
            required.remove(k)
        except ValueError:
            pass

        if k not in allowed:
            extra.append(k)

    # Check for sanity of specified args
    err = ""
    if len(extra):
        err += f"\n\nThe following arguments are not compatible with mode {mode}:\n"
        for e in extra:
            err += f"    {e}\n"

    if len(required):
        err += f"\n\nThe following arguments are required for mode {mode}:\n"
        for r in required:
            err += f"    {r}\n"

    if err != "":
        err += "\n\n"
        raise ValueError(err)

    # Call actual function
    func(**kwargs)



if __name__ == "__main__":
    main(sys.argv[1:])
