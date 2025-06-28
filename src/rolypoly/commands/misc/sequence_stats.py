from pathlib import Path
import hashlib
import json
from collections import Counter

import polars as pl
import rich_click as click
from rich.console import Console

from rolypoly.utils.citation_reminder import remind_citations
from rolypoly.utils.fax import (
    read_fasta_df,
)

# show all columns
pl.Config().set_tbl_cols(-1)

# global console
console = Console()


@click.command()
@click.option(
    "-i",
    "--input",
    required=True,
    type=click.Path(exists=True),
    help="Input file (fasta, fa, fna, faa)",
)
@click.option(
    "-agg",
    "--aggregate",
    default=False,
    type=bool,
    help="aggregate statistics across all sequences",
)
@click.option(
    "-o",
    "--output",
    default="rp_sequence_stats.txt",
    type=click.Path(exists=False),
    help="Output path",
)
@click.option(
    "--log-file",
    default="command.log",
    type=click.Path(exists=False),
    help="Path to log file",
    hidden=True,
)
@click.option(
    "--log-level",
    hidden=True,
    default="INFO",
    help="Log level",
)
@click.option(
    "--min_length",
    hidden=True,
    default=None,
    type=int,
    help="minimum sequence length to consider",
)
@click.option(
    "--max_length",
    hidden=True,
    default=None,
    type=int,
    help="maximum sequence length to consider",
)
@click.option(
    "--format",
    default="tsv",
    type=click.Choice(case_sensitive=False, choices=["csv", "tsv", "md", "parquet"]),
    help="output format, either a parquet/csv/tsv file with the data or a markdown file with summary statistics",
)
@click.option(
    "-f",
    "--fields",
    type=click.Choice(case_sensitive=False, choices=["length", "gc_content", "n_count", "hash", "kmer_freq"]),
    multiple=True,
    default=["length", "gc_content", "n_count", "hash"],
    help="""
              comma-separated list of fields to include.  
              Available:
              length - mandatory
              gc_content - percentage of GC nucleotides
              n_count - total number of Ns 
              hash - md5 hash of the sequence
              kmer_freq - k-mer frequencies (k=3 by default)
              """,
)
def sequence_stats(
    input,
    aggregate,
    output,
    log_file,
    log_level,
    min_length,
    max_length,
    format,
    fields,
):
    """Calculate sequence statistics using Polars expressions"""
    from rolypoly.utils.loggit import log_start_info, setup_logging

    logger = setup_logging(log_file, log_level)
    log_start_info(logger, locals())

    output_path = Path(output)

    # Read sequences into DataFrame
    df = read_fasta_df(input)
    total_seqs = len(df)
    df = df.with_columns(pl.col("sequence").str.len_chars().alias("length"))
    
    logger.info(f"Read {total_seqs} sequences from {input}")

    # Apply length filters
    if min_length:
        df = df.filter(pl.col("length") >= min_length)
        logger.info(f"Applied minimum length filter: {min_length}")
    if max_length:
        df = df.filter(pl.col("length") <= max_length)
        logger.info(f"Applied maximum length filter: {max_length}")

    filtered_seqs = len(df)
    logger.info(f"After filtering: {filtered_seqs} sequences")

    # Define available fields and their dependencies
    field_options = {
        "length": {"desc": "Sequence length"},
        "gc_content": {"desc": "GC content percentage"},
        "n_count": {"desc": "Count of Ns in sequence"},
        "hash": {"desc": "Sequence hash (MD5)"},
        "kmer_freq": {"desc": "K-mer frequencies", "complex": True},
    }

    # Parse fields - fields is already a tuple from multiple=True
    selected_fields = list(fields) if fields else ["length", "gc_content", "n_count", "hash"]
    
    # Always include length for summaries
    if "length" not in selected_fields:
        selected_fields.append("length")

    logger.info(f"Selected fields: {', '.join(selected_fields)}")

    # Helper function for kmer calculation
    def calculate_kmers(sequence, k=3):
        """Calculate k-mer frequencies for a sequence"""
        if len(sequence) < k:
            return {}
        kmers = [sequence[i:i+k] for i in range(len(sequence) - k + 1)]
        return dict(Counter(kmers))

    # Helper function for MD5 hash
    def calculate_hash(sequence):
        """Calculate MD5 hash of sequence"""
        return hashlib.md5(sequence.encode()).hexdigest()

    # Calculate selected statistics
    if "gc_content" in selected_fields:
        df = df.with_columns(
            (pl.col("sequence").str.count_matches(r"[GCgc]") / pl.col("length") * 100.0)
            .alias("gc_content")
        )
    
    if "n_count" in selected_fields:
        df = df.with_columns(
            pl.col("sequence").str.count_matches(r"[Nn]").alias("n_count")
        )
    
    if "hash" in selected_fields:
        df = df.with_columns(
            pl.col("sequence").map_elements(calculate_hash, return_dtype=pl.String).alias("hash")
        )

    if "kmer_freq" in selected_fields:
        df = df.with_columns(
            pl.col("sequence").map_elements(
                lambda x: json.dumps(calculate_kmers(x, k=3)), 
                return_dtype=pl.String
            ).alias("kmer_frequencies")
        )

    # Select only the fields we want to output
    output_columns = ["header"] + selected_fields
    if "kmer_freq" in selected_fields:
        # Replace kmer_freq with kmer_frequencies in output columns
        output_columns = [col if col != "kmer_freq" else "kmer_frequencies" for col in output_columns]
    
    df_output = df.select([col for col in output_columns if col in df.columns])

    if aggregate:
        # Create aggregated statistics
        agg_stats = {}
        
        if "length" in selected_fields:
            length_stats = df.select([
                pl.col("length").min().alias("min_length"),
                pl.col("length").max().alias("max_length"),
                pl.col("length").mean().alias("mean_length"),
                pl.col("length").median().alias("median_length"),
                pl.col("length").std().alias("std_length"),
                pl.col("length").sum().alias("total_length")
            ]).to_dicts()[0]
            agg_stats.update(length_stats)
        
        if "gc_content" in selected_fields:
            gc_stats = df.select([
                pl.col("gc_content").min().alias("min_gc"),
                pl.col("gc_content").max().alias("max_gc"),
                pl.col("gc_content").mean().alias("mean_gc"),
                pl.col("gc_content").median().alias("median_gc"),
                pl.col("gc_content").std().alias("std_gc")
            ]).to_dicts()[0]
            agg_stats.update(gc_stats)
        
        if "n_count" in selected_fields:
            n_stats = df.select([
                pl.col("n_count").min().alias("min_n_count"),
                pl.col("n_count").max().alias("max_n_count"),
                pl.col("n_count").mean().alias("mean_n_count"),
                pl.col("n_count").sum().alias("total_n_count")
            ]).to_dicts()[0]
            agg_stats.update(n_stats)
        
        agg_stats["total_sequences"] = filtered_seqs
        agg_stats["sequences_before_filter"] = total_seqs
        
        # Convert to DataFrame for consistent output
        df_output = pl.DataFrame([agg_stats])

    # Output results
    if format.lower() == "parquet":
        df_output.write_parquet(output_path)
        logger.info(f"Results written to {output_path} (parquet format)")
    
    elif format.lower() == "csv":
        df_output.write_csv(output_path)
        logger.info(f"Results written to {output_path} (CSV format)")
    
    elif format.lower() == "tsv":
        df_output.write_csv(output_path, separator="\t")
        logger.info(f"Results written to {output_path} (TSV format)")
    
    elif format.lower() == "md":
        # Create markdown summary
        md_content = f"# Sequence Statistics Report\n\n"
        md_content += f"**Input file:** {input}\n"
        md_content += f"**Total sequences:** {total_seqs}\n"
        md_content += f"**Sequences after filtering:** {filtered_seqs}\n\n"
        
        if aggregate:
            md_content += "## Aggregate Statistics\n\n"
            for key, value in agg_stats.items():
                if isinstance(value, float):
                    md_content += f"- **{key}:** {value:.2f}\n"
                else:
                    md_content += f"- **{key}:** {value}\n"
        else:
            md_content += "## Summary Statistics\n\n"
            # Add basic summary stats even when not aggregating
            if "length" in selected_fields:
                length_summary = df.select([
                    pl.col("length").min().alias("min"),
                    pl.col("length").max().alias("max"),
                    pl.col("length").mean().alias("mean"),
                    pl.col("length").median().alias("median")
                ]).to_dicts()[0]
                
                md_content += "### Length Statistics\n"
                for stat, value in length_summary.items():
                    md_content += f"- **{stat}:** {value:.2f}\n"
                md_content += "\n"
            
            md_content += f"\n## First 10 sequences\n\n"
            md_content += df_output.head(10).to_pandas().to_markdown(index=False)
        
        with open(output_path, 'w') as f:
            f.write(md_content)
        logger.info(f"Markdown report written to {output_path}")

    # Display summary to console
    console.print(f"\n[bold green]✓[/bold green] Processed {filtered_seqs} sequences")
    console.print(f"[bold blue]Output:[/bold blue] {output_path}")
    
    if not aggregate and format.lower() != "md":
        console.print("\n[bold]First 5 rows:[/bold]")
        console.print(df_output.head(5))

    # Remind about citations
    tools = ["polars"]  # Tools used in this analysis
    remind_citations(tools)
    
    logger.info("Sequence statistics calculation completed successfully")
