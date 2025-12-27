"""
Genome downloading using pre-computed rRNA to genome mappings.

Replaces ncbi-datasets + taxonkit dependencies with direct FTP downloads
using pre-computed taxonomy-based mappings. Prefers transcript files over
genome files and constrains matches to the same genus.

Usage:
    # Fetch genomes for a list of taxids
    fetch_genomes_from_mapping(
        taxids=[9606, 10090],
        mapping_path="data/contam/rrna/rrna_to_genome_mapping.parquet",
        output_file="genomes.fasta",
        prefer_transcript=True
    )
    
    # Drop-in replacement for old fetch_genomes() # internally it calls the one above...
    fetch_genomes_from_stats_file(
        stats_file="bbmap_stats.txt",
        mapping_path="mapping.parquet",
        output_file="genomes.fasta",
        max_genomes=5
    )

Test:
    python src/tests/test_genome_fetch_mapping.py
"""

import gzip
import os
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
from urllib.parse import urlparse
from urllib.request import urlretrieve

import polars as pl
from rich.progress import Progress, SpinnerColumn, BarColumn, TextColumn, TimeElapsedColumn, TaskID

from rolypoly.utils.bio.sequences import remove_duplicates
from rolypoly.utils.logging.loggit import get_logger


def load_rrna_genome_mapping(mapping_path: str) -> pl.DataFrame:
    """Load the pre-computed rRNA to genome mapping table.
    
    Args:
        mapping_path: Path to the parquet file containing mappings
        
    Returns:
        Polars DataFrame with mapping information
    """
    if not Path(mapping_path).exists():
        raise FileNotFoundError(
            f"Mapping file not found: {mapping_path}\n"
            "Please run the rrna_genome_mapping_taxonomy notebook to generate it."
        )
    
    return pl.read_parquet(mapping_path)


def get_ftp_path_for_taxid(
    taxid: int, 
    mapping_df: pl.DataFrame
) -> Optional[str]:
    """Get the best FTP path for a given taxonomy ID.
    
    Args:
        taxid: NCBI taxonomy ID
        mapping_df: The mapping DataFrame
        
    Returns:
        FTP path string, or None if no mapping found
    """
    # Query the mapping for this taxid
    result = mapping_df.filter(pl.col("query_tax_id") == taxid)
    
    if result.height == 0:
        return None
    
    row = result.row(0, named=True)
    
    # Priority: self > ancestor > leaf
    if row.get("self_has_data"):
        # For self, we need to look up in the original assembly data
        # This would require joining with the assembly summary
        pass
    
    # Try ancestor first if it has data
    if row.get("ancestor_ftp_path"):
        return row["ancestor_ftp_path"]
    
    # Fall back to leaf relative
    if row.get("leaf_ftp_path"):
        return row["leaf_ftp_path"]
    
    return None


def download_from_ftp_path(
    ftp_path: str,
    output_dir: Path,
    prefer_transcript: bool = True,
    overwrite: bool = False,
    logger=None
) -> Optional[Path]:
    """Download sequence files from an NCBI FTP path.
    
    Args:
        ftp_path: Base FTP path (e.g., https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/...)
        output_dir: Directory to save downloaded files
        _cds_from_genomic.fna.gz
        prefer_transcript: If True, try to download *_cds_from_genomic.fna.gz first
        overwrite: If True, re-download even if file exists
        
    Returns:
        Path to the downloaded file, or None if download failed
    """
    logger = get_logger(logger)
    
    if not ftp_path:
        return None
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Extract the assembly name from the FTP path
    # e.g., GCF_000001405.40_GRCh38.p14 from the URL
    assembly_name = Path(urlparse(ftp_path).path).name
    
    # File types to try, in order of preference
    if prefer_transcript:
        file_suffixes = [
            "_cds_from_genomic.fna.gz",  # Transcript sequences
            "_genomic.fna.gz",           # Genome sequences
        ]
    else:
        file_suffixes = [
            "_genomic.fna.gz",           # Genome sequences
            "_cds_from_genomic.fna.gz",  # Transcript sequences
        ]
    
    for suffix in file_suffixes:
        filename = f"{assembly_name}{suffix}"
        url = f"{ftp_path}/{filename}"
        output_file = output_dir / filename
        
        # Check if file already exists
        if output_file.exists() and not overwrite:
            logger.info(f"File already exists: {output_file}")
            return output_file
        
        try:
            logger.info(f"Downloading: {url}")
            urlretrieve(url, output_file)
            
            # Verify the download (basic check)
            if output_file.exists() and output_file.stat().st_size > 0:
                logger.info(f"Successfully downloaded: {output_file}")
                return output_file
            else:
                logger.warning(f"Downloaded file is empty, trying next option")
                output_file.unlink(missing_ok=True)
                
        except Exception as e:
            logger.warning(f"Could not download {filename}: {e}")
            output_file.unlink(missing_ok=True)
            continue
    
    logger.error(f"No files could be downloaded from: {ftp_path}")
    return None


def fetch_genomes_from_mapping(
    taxids: List[int],
    mapping_path: str,
    output_file: str,
    temp_dir: Optional[str] = None,
    prefer_transcript: bool = True,
    threads: int = 1,
    overwrite: bool = False,
    exclude_viral: bool = True,
    clean_headers = True, # currently part of the exclude virus if else, for now no need to take it out.
    logger=None
) -> None:
    """Fetch genome/transcript sequences for a list of taxids using pre-computed mappings.
    
    This replaces the ncbi-datasets + taxonkit workflow with a direct FTP download
    approach using pre-computed rRNA to genome mappings.
    
    Args:
        taxids: List of NCBI taxonomy IDs
        mapping_path: Path to the rRNA genome mapping parquet file
        output_file: Output fasta file path
        temp_dir: Temporary directory for downloads (default: creates one in current dir)
        prefer_transcript: If True, prefer transcript files over genome files
        threads: Number of parallel downloads
        overwrite: If True, re-download existing files
        exclude_viral: If True, filter out viral sequences from final output
        remove_lcl: clean the "lcl|" prefix from fasta headers (can cause issues if kept)

    """
    logger = get_logger(logger)
    
    # Setup directories
    if temp_dir is None:
        temp_dir = "tmp_genome_downloads"
    temp_path = Path(temp_dir)
    temp_path.mkdir(parents=True, exist_ok=True)
    
    # Load mapping
    logger.info(f"Loading mapping from: {mapping_path}")
    mapping_df = load_rrna_genome_mapping(mapping_path)
    
    # Track downloaded files
    downloaded_files: List[Path] = []
    unmapped_taxids: List[int] = []
    
    # Worker function for parallel processing
    def process_taxid(taxid: int, progress: Optional[Progress] = None, task: Optional[TaskID] = None) -> Tuple[int, Optional[Path]]:
        """Process a single taxid and return (taxid, downloaded_path)."""
        ftp_path = get_ftp_path_for_taxid(taxid, mapping_df)
        
        if ftp_path is None:
            if progress and task is not None:
                progress.update(task, advance=1)
            return (taxid, None)
        
        downloaded = download_from_ftp_path(
            ftp_path=ftp_path,
            output_dir=temp_path,
            prefer_transcript=prefer_transcript,
            overwrite=overwrite,
            logger=logger
        )
        
        if progress and task is not None:
            progress.update(task, advance=1)
        
        return (taxid, downloaded)
    
    # Process each taxid
    logger.info(f"Fetching sequences for {len(taxids)} taxids using {threads} thread(s)...")
    
    # Create progress bar
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        TextColumn("({task.completed}/{task.total})"),
        TimeElapsedColumn(),
    ) as progress:
        task = progress.add_task("[cyan]Downloading genomes...", total=len(taxids))
        
        if threads > 1:
            # Parallel processing
            with ThreadPoolExecutor(max_workers=threads) as executor:
                futures = {executor.submit(process_taxid, taxid, progress, task): taxid for taxid in taxids}
                
                for future in as_completed(futures):
                    taxid = futures[future]
                    try:
                        result_taxid, downloaded = future.result()
                        
                        if downloaded is None:
                            unmapped_taxids.append(result_taxid)
                            logger.warning(f"No mapping found for taxid: {result_taxid}")
                        else:
                            downloaded_files.append(downloaded)
                            
                    except Exception as e:
                        logger.error(f"Error processing taxid {taxid}: {e}")
                        unmapped_taxids.append(taxid)
                        progress.update(task, advance=1)
        else:
            # Sequential processing
            for taxid in taxids:
                result_taxid, downloaded = process_taxid(taxid, progress, task)
                
                if downloaded is None:
                    unmapped_taxids.append(result_taxid)
                    logger.warning(f"No mapping found for taxid: {result_taxid}")
                else:
                    downloaded_files.append(downloaded)
    # Report statistics
    logger.info(f"Downloaded {len(downloaded_files)} genome files")
    if unmapped_taxids:
        logger.warning(f"Could not map {len(unmapped_taxids)} taxids")
    
    if not downloaded_files:
        logger.error("No genome files were downloaded. Exiting.")
        return
    
    # Decompress and concatenate files
    logger.info("Processing downloaded files...")
    temp_concat = temp_path / "concat_genomes.fasta"
    
    with open(temp_concat, "w") as out_f:
        for i, gz_file in enumerate(downloaded_files, 1):
            logger.debug(f"Decompressing {gz_file.name} ({i}/{len(downloaded_files)})")
            try:
                with gzip.open(gz_file, "rt") as in_f:
                    shutil.copyfileobj(in_f, out_f)
            except Exception as e:
                logger.error(f"Error processing {gz_file}: {e}")
    
    # Deduplicate sequences
    logger.info("Deduplicating sequences...")
    temp_dedup = temp_path / "dedup_genomes.fasta"
    
    try:
        remove_duplicates(
            input_file=str(temp_concat),
            output_file=str(temp_dedup),
            by="seq",
            logger=logger
        )
    except Exception as e:
        logger.error(f"Error during deduplication: {e}")
        shutil.copy(temp_concat, temp_dedup)
    
    # Filter viral sequences if requested
    if exclude_viral:
        logger.info("Filtering out headers with 'phage', 'virus', and 'viral' in headers...")
        from rolypoly.utils.bio.sequences import filter_fasta_by_headers, clean_fasta_headers
        
        temp_dedup1 = str(temp_dedup) + "1"
        filter_fasta_by_headers(
            str(temp_dedup),
            ["virus", "viral", "phage"],
            temp_dedup1,
            wrap=True,
            invert=True,
        )

        clean_fasta_headers(
            fasta_file=temp_dedup1,
            drop_from_space=True,
            output_file=output_file,
            strip_prefix="lcl|",
            strip_suffix="bla bla bla"
        )

    else:
        shutil.copy(temp_dedup, output_file)
    
    logger.info(f"Final output written to: {output_file}")
    
    # Cleanup
    if not overwrite:  # Keep temp files if we might need them
        shutil.rmtree(temp_path, ignore_errors=True)


def fetch_genomes_from_stats_file(
    stats_file: str,
    mapping_path: str,
    output_file: str,
    max_genomes: int = 5,
    logger=None,
    **kwargs
) -> None:
    """Fetch genomes based on a BBMap stats file, using the mapping table.
    
    This is a drop-in replacement for the original fetch_genomes function.
    
    Args:
        stats_file: BBMap stats file with taxonomic information
        mapping_path: Path to the rRNA genome mapping parquet file
        output_file: Output fasta file path
        max_genomes: Maximum number of genomes to fetch
        **kwargs: Additional arguments passed to fetch_genomes_from_mapping
    """
    logger = get_logger(logger)
    
    # Parse the stats file to extract taxon names
    logger.info(f"Parsing stats file: {stats_file}")
    
    # Read the BBMap stats file (skip first 4 lines, then parse)
    with open(stats_file, "r") as f:
        lines = f.readlines()[4:]
    
    taxons: Set[str] = set()
    for line in lines:
        taxon = line.split(sep=";")[-1]
        taxon = taxon.split(sep="\t")[0]
        
        # Skip unwanted taxons - in the new approach this shouldn't be necessary but keeping for safety
        if any(
            word in taxon.lower()
            for word in [
                "meta",
                "uncultured",
                "unidentified",
                "synthetic",
                "construct",
                "coli",
            ]
        ):
            continue
        
        # Clean up taxon name
        if "PREDICTED: " in taxon:
            taxon = taxon.split(sep=": ")[1]
        
        # Remove rRNA-specific suffixes
        for suffix in [" 16S", " 18S", " 28S", " small subunit ", " large subunit "]:
            taxon = taxon.split(sep=suffix)[0]
        
        taxons.add(taxon)
        
        if len(taxons) >= max_genomes:
            break
    
    logger.info(f"Found {len(taxons)} taxons to fetch")
    
    # Now we need to map taxon names to taxids
    # We can do this by loading the NCBI names.dmp file
    # OR we can use the mapping table which already has names
    
    # Load the mapping to get name -> taxid mapping
    mapping_df = load_rrna_genome_mapping(mapping_path)
    
    # Create a reverse lookup: name -> taxid
    name_to_taxid = {}
    for row in mapping_df.select(["query_tax_id", "query_name"]).unique().iter_rows(named=True):
        name_to_taxid[row["query_name"]] = row["query_tax_id"]
    
    # Also check ancestor and leaf names
    for row in mapping_df.select(["ancestor_tax_id", "ancestor_name"]).filter(
        pl.col("ancestor_name").is_not_null()
    ).unique().iter_rows(named=True):
        name_to_taxid[row["ancestor_name"]] = row["ancestor_tax_id"]
    
    for row in mapping_df.select(["leaf_tax_id", "leaf_name"]).filter(
        pl.col("leaf_name").is_not_null()
    ).unique().iter_rows(named=True):
        name_to_taxid[row["leaf_name"]] = row["leaf_tax_id"]
    
    # Map taxon names to taxids
    taxids: List[int] = []
    unmapped_names: List[str] = []
    
    for taxon in taxons:
        if taxon in name_to_taxid:
            taxids.append(name_to_taxid[taxon])
        else:
            # Try partial matching (genus name)
            genus = taxon.split()[0] if " " in taxon else taxon
            if genus in name_to_taxid:
                taxids.append(name_to_taxid[genus])
            else:
                unmapped_names.append(taxon)
    
    if unmapped_names:
        logger.warning(f"Could not map {len(unmapped_names)} taxon names to taxids")
        logger.warning(f"Examples: {unmapped_names[:5]}")
    
    # Fetch genomes for the mapped taxids
    fetch_genomes_from_mapping(
        taxids=taxids,
        mapping_path=mapping_path,
        output_file=output_file,
        **kwargs
    )
