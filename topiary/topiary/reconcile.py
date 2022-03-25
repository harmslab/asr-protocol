
import topiary
import numpy as np
import os, shutil, re

def _annotate_tree_with_calls(df,tree,work_on_copy=True):
    """
    Annotate the leaves of an ete3 tree with information extracted from a
    topiary dataframe.

    df: topiary data frame
    tree: tree (newick, dendropy, or ete3)
    work_on_copy: whether or not to work on copy of existing tree.

    returns: ete3 with annotations.
    """
    # Copy tree -- do not operate on input tree directly
    tree = topiary.util.load_tree(tree)
    if work_on_copy:
        tree = tree.copy(method="deepcopy")

    # Create dictionaries mapping uid to species, paralog, ott, and call.
    out_dict = {}
    for i in range(len(df)):

        uid = df.uid.iloc[i]
        species = df.species.iloc[i]
        ott = f”{int(np.round(df.ott.iloc[i],0)):d}”
        try:
            paralog = df.paralog.iloc[i]
        except AttributeError:
            paralog = df.protein.iloc[i]
        call = f"{ott}|{paralog}"

        out_dict[uid] = {"species":species,
                         "paralog":paralog,
                         "ott":ott,
                         "call":call}

    # Go through tree and assign leaves their calls etc.
    for node in tree.get_leaves():
        try:
            for k in out_dict[node.name]:
                node.add_feature(k,out_dict[node.name][k])
        except KeyError:
            err = f"uid {node.name} in tree but not data frame!\n"
            raise ValueError(err)

    return tree

def _create_generax_input(df,gene_tree):
    """
    Take a dataframe and generate the data structures necessary for a GeneRax
    run.

    df: topiary data frame.
    gene_tree: gene tree with uid taxon names.

    returns consistently named gene_tree, species_tree, and link_dict
    """

    # Only look at sequences flagged to keep
    df = df.loc[df.keep].copy() # copy is to avoid assign-to-copy warning

    # Annotate gene tree with uid, ott, etc. and arbitrarily resolve any
    # polytomies.
    gene_tree = _annotate_tree_with_calls(df,gene_tree)
    gene_tree.resolve_polytomy()

    uid_in_gene_tree = []
    link_dict = {}
    for l in gene_tree.get_leaves():

        # For generating mapping file, record which uid are associated with what ott
        try:
            link_dict[l.ott].append(l.name)
        except KeyError:
            link_dict[l.ott] = [l.name]

        # Make sure this uid is actually in the df once and only once
        if np.sum(df.uid == l.name) != 1:
            err = f"tree taxon {l.name} either missing or duplicated in data frame\n"
            raise ValueError(err)

        # Record that we saw this uid
        uid_in_gene_tree.append(l.name)

    # Make df only have uid seen (will automatically trim down to only ott
    # of interest)
    mask = np.array([u in uid_in_gene_tree for u in df.uid],dtype=np.bool)
    df = df.loc[mask]

    # Get species tree.
    species_tree = topiary.get_species_tree(df)

    # Resolve polytomies and make sure all branch lenghts/supports have values
    species_tree.resolve_polytomy()
    for n in species_tree.traverse():
        if n.dist != 1:
            n.dist = 1
        if n.support != 1:
            n.support = 1

    return gene_tree, species_tree, link_dict

def _write_generax_input(df,gene_tree,species_tree,link_dict,model,out_dir):
    """
    Write out files for running generax. The contents of this directory can be
    run on the command line by:

    df: topiary data frame
    gene_tree: gene tree returned from _create_generax_input
    species_tree: species tree returned from _create_generax_input
    link_dict: link_dict returned from _create_generax_input
    model: phylogenetic model to use (should match ml tree)
    out_dir: output directory to write files
    """

    # Construct the control file for generax
    control_out = []
    control_out.append("[FAMILIES]")
    control_out.append("- reconcile")
    control_out.append("starting_gene_tree = gene_tree.newick")
    control_out.append("alignment = alignment.phy")
    control_out.append("mapping = mapping.link")
    control_out.append(f"subst_model = {model}")

    # Write out control file
    f = open(os.path.join(out_dir,"control.txt"),"w")
    f.write("\n".join(control_out))
    f.close()

    # Write out .phy file
    topiary.write_phy(df,os.path.join(out_dir,"alignment.phy"),
                      seq_column="alignment")

    # Write out gene tree
    gene_tree.write(outfile=os.path.join(out_dir,"gene_tree.newick"),
                    format=5)

    # Write out species tree
    species_tree.write(outfile=os.path.join(out_dir,"species_tree.newick"))

    # Write out link file
    f = open(os.path.join(out_dir,"mapping.link"),"w")
    for k in link_dict:
        f.write(f"{k}:")
        f.write(";".join(link_dict[k]))
        f.write("\n")
    f.close()

    # Command that would run generax
    f = open(os.path.join(out_dir,"run.sh"),"w")
    f.write(re.sub("    ","",
    """#!/bin/bash
    num_nodes="${1}"
    if [[ ! "${num_nodes}" ]]; then
        start_cmd=""
    else
        start_cmd="mpirun -np ${num_nodes} "
    fi
    ${start_cmd}generax -f control.txt -s species_tree.newick -p result -r UndatedDL

    """))

    f.close()

def _write_script(script_string,out_file,strip_indent=1):
    """
    Write a script given a string, stripping leading indents due to them being
    in python code. (Yes, this is janky).

    script_string: string with script
    out_file: where to write out script
    strip_indent: strip indents to this level (assumes 4 space indent)
    """

    f = open(out_file,"w")
    for l in script_string.split("\n"):
        if l.strip() == "":
            continue

        f.write(re.sub("    ","",l,count=strip_indent))
        f.write("\n")
    f.close()


def setup_generax(df,gene_tree,model,out_dir,dir_with_bootstraps=None):

    if os.path.isdir(out_dir):
        err = f"out_dir '{out_dir}' already exists."
        raise FileExistsError(err)

    os.mkdir(out_dir)
    template_dir = os.path.join(out_dir,"ml")
    os.mkdir(template_dir)

    # Create generax data structures
    gene_tree, species_tree, link_dict = _create_generax_input(df,gene_tree)

    # Actually write out generax output
    _write_generax_input(df,gene_tree,species_tree,link_dict,model,template_dir)

    if dir_with_bootstraps is not None:

        bs_tree_file = os.path.join(dir_with_bootstraps,"alignment.raxml.bootstraps")
        bs_trees = []
        with open(bs_tree_file) as f:
            for line in f:
                bs_trees.append(line.strip())


        for i in range(len(bs_trees)):

            # Create the bootstrap directory from the template directory
            bs_dir = os.path.join(out_dir,f"bs_{i+1:04}")
            shutil.copytree(template_dir,bs_dir)

            # Copy in bootstrap msa
            bs_msa = os.path.join(dir_with_bootstraps,
                                  f"alignment.raxml.bootstrapMSA.{i+1}.phy")
            shutil.copy(bs_msa,os.path.join(bs_dir,"alignment.phy"))

            # Copy in the bootstrap tree
            f = open(os.path.join(bs_dir,"gene_tree.newick"),"w")
            f.write(bs_trees[i])
            f.close()

def write_hpc_templates(out_dir):
    """
    Janky function that dumps files useful for hpc runs (slurm) into input
    directory.

    out_dir: output directory
    """

    this_str = \
    """#!/bin/bash -l
    #SBATCH --account=harmslab      ### change this to your actual account for charging
    #SBATCH --job-name=generax      ### job name
    #SBATCH --output=hostname.out   ### file in which to store job stdout
    #SBATCH --error=hostname.err    ### file in which to store job stderr
    #SBATCH --partition=short       ### can be short long fat longfat
    #SBATCH --time=01-00:00:00      ### Run for 1 day
    #SBATCH --ntasks=3              ### 3 mpi tasks

    module load openmpi

    # Loop over all directories specified on the command line
    for dir in "$@"; do

        # Go into the directory and launch with 3 threads. (This 3 should match
        # --ntasks above)
        cd ${dir}
        bash run.sh 3
        cd ../

    done
    """
    _write_script(this_str,os.path.join(out_dir,"_generax_bootstrap.srun"))

    this_str = \
    """
    #!/bin/bash

    # How many replicates to run per job. Lower number runs faster. If you have
    # 1000 bootstrap replicates and select 5, this will run 1000/5 = 200 jobs,
    # each with 5 replicates.
    num_per_node="${1}"
    if [[ ! "${num_per_node}" ]]; then
        echo "launch_generax-bs.sh num_to_run_per_node"
        exit
    fi

    # Start with command string planning to run on ml
    cmd_string="ml "
    counter=0

    # Go over bootstrap directories
    for bs_rep in bs_*; do

        # Grab the name of the bootstrap directory
        cmd_string="${cmd_string} ${bs_rep}"
        counter=$(($counter + 1))

        # If our counter is big enough, launch the job with these directories and
        # reset the counter.
        if [[ ${counter} -ge num_per_node ]]; then
            echo ${cmd_string}
            qsub _generax_bootstrap.srun ${cmd_string}

            cmd_string=""
            counter=0
        fi
    done
    """
    _write_script(this_str,os.path.join(out_dir,"00_launch_generax_bootstrap.sh"))


    this_str = \
    """
    #!/bin/bash

    rm -f bs-trees.newick
    for dir in bs_*; do
        cat ${dir}/result/results/reconcile/geneTree.newick >> bs-trees.newick
        echo "" >> bs-trees.newick
    done

    raxml-ng.dev --support --tree ml/result/results/reconcile/geneTree.newick --bs-trees bs-trees.newick --redo

    cp ml/result/results/reconcile/geneTree.newick reconciled-tree.newick
    cp ml/result/results/reconcile/geneTree.newick.raxml.support reconciled-tree-with-supports.newick
    """
    _write_script(this_str,os.path.join(out_dir,"01_assemble_generax_bootstrap.sh"))
