"""Translation and ORF prediction functions."""

import multiprocessing.pool
from needletail import parse_fastx_file

from rolypoly.utils.various import run_command_comp
from typing import Union
from pathlib import Path

def translate_6frx_seqkit(
    input_file: str,
    output_file: str,
    threads: int,
    min_orf_length: int = 0,
) -> None:
    """Translate nucleotide sequences in all 6 reading frames using seqkit.

    Args:
        input_file (str): Path to input nucleotide FASTA file
        output_file (str): Path to output amino acid FASTA file
        threads (int): Number of CPU threads to use

    Note:
        Requires seqkit to be installed and available in PATH.
        The output sequences are formatted with 20000bp line width.
    """
    # import subprocess as sp

    # command = f"seqkit translate -x -F --clean --min-len {min_orf_length} -w 0 -f 6 {input_file} --id-regexp '(\\*)' --clean  --threads {threads} > {output_file}"
    run_command_comp(base_cmd="seqkit translate",
                     assign_operator="=",
                     prefix_style="double",
                     params={
                         "allow-unknown-codon": True,
                         "clean": True, 
                         "min-len": min_orf_length,
                         "line-width": 0,
                         "frame": 6,
                         "id-regexp": "(\\*)",
                         "clean": True,
                         "threads": threads},
                     positional_args=[f"{input_file} --out-file {output_file}"],
                     positional_args_location="end")
    # sp.run(command, shell=True, check=True)


def translate_with_bbmap(input_file: str, output_file: str, threads: int) -> None:
    """Translate nucleotide sequences using BBMap's callgenes.sh

    Args:
        input_file (str): Path to input nucleotide FASTA file
        output_file (str): Path to output amino acid FASTA file
        threads (int): Number of CPU threads to use

    Note:
        - Requires BBMap to be installed and available in PATH (should be done via bbmapy)
        - Generates both protein sequences (.faa) and gene annotations (.gff)
        - The GFF output file is named by replacing .faa with .gff
    """
    import subprocess as sp

    gff_o = output_file.replace(".faa", ".gff")
    command = (
        f"callgenes.sh threads={threads} in={input_file} outa={output_file} out={gff_o}"
    )
    sp.run(command, shell=True, check=True)


def pyro_predict_orfs(
    input_file: str,
    output_file: str,
    threads: int,
    min_gene_length: int = 30,
    genetic_code: int = 11,  # NOT USED YET # TODO: add SUPPORT for this.
    model: str = "1",  # NOT USED YET/at all.
) -> None:
    """Predict and translate Open Reading Frames using Pyrodigal.

    Uses either Pyrodigal-GV (optimized for viruses) or standard Pyrodigal
    to predict and translate ORFs from nucleotide sequences.

    Args:
        input_file (str): Path to input nucleotide FASTA file
        output_file (str): Path to output amino acid FASTA file
        threads (int): Number of CPU threads to use
        genetic_code (int, optional): Genetic code table to use (Standard/Bacterial) (NOT USED YET).

    Note:
        - Creates both protein sequences (.faa) and gene annotations (.gff)
        - genetic_code is 11 for standard/bacterial
    """
    import pyrodigal_gv as pyro_gv

    sequences = []
    ids = []
    for record in parse_fastx_file(input_file):
        sequences.append((record.seq)) # type: ignore
        ids.append((record.id)) # type: ignore

    gene_finder = pyro_gv.ViralGeneFinder(
        meta=True,
        min_gene=min_gene_length,
        
    )  # a single gv gene finder object

    with multiprocessing.pool.Pool(processes=threads) as pool:
        orfs = pool.map(gene_finder.find_genes, sequences)

    with open(output_file, "w") as dst:
        for i, orf in enumerate(orfs):
            orf.write_translations(dst, sequence_id=ids[i], width=111110)

    with open(output_file.replace(".faa", ".gff"), "w") as dst:
        for i, orf in enumerate(orfs):
            orf.write_gff(dst, sequence_id=ids[i], full_id=True) 
            
            
def predict_orfs_orffinder(
    input_fasta: Union[str, Path],
    output_file: Union[str, Path],
    min_orf_length: int,
    genetic_code: int,
    start_codon: int = 1,
    strand: str = "both",
    outfmt: int = 1,
    ignore_nested: bool = False,
) -> None:
    run_command_comp(
        "ORFfinder",
        params={
            "in": str(input_fasta),
            "out": str(output_file),
            "ml": min_orf_length, # orfinder automatically replaces values below 30 to 30.
            "s": start_codon, # ORF start codon to use, 0 is atg only, 1 atg + alt start codons
            "g": genetic_code,
            "n": "false" if ignore_nested else "true", # do not ignore nested ORFs
            "strand": strand, # both is plus and minus.
            "outfmt": outfmt, # 1 is fasta, 3 is feature table
        },prefix_style="single"
    )
