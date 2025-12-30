import shutil
from pathlib import Path

from rolypoly.utils.bio.interval_ops import mask_nuc_range_from_sam


def test_mask_nuc_range_from_sam(tmp_path: Path):
    """Copy provided SAM/FASTA to tmp, run masking, and assert Ns present."""
    sam_src = Path(
        "tests/test_sam_for_masking.sam"
    )
    fasta_src = Path(
        "tests/test_fasta_for_masking.fasta"
    )

    assert sam_src.exists(), f"SAM input not found: {sam_src}"
    assert fasta_src.exists(), f"FASTA input not found: {fasta_src}"

    sam = tmp_path / "in.sam"
    fasta = tmp_path / "in.fasta"
    out = tmp_path / "out.fasta"

    shutil.copy(sam_src, sam)
    shutil.copy(fasta_src, fasta)

    mask_nuc_range_from_sam(str(fasta), str(sam), str(out))

    assert out.exists(), "Output FASTA was not created"

    # Check at least one sequence contains masked 'N' bases
    with open(out, "r") as fh:
        seqs = []
        cur = []
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                    cur = []
            else:
                cur.append(line.strip())
        if cur:
            seqs.append("".join(cur))

    assert any("N" in s for s in seqs), "No masked Ns found in output sequences" #NOT GOOD THIS ASSUMES THE INPUT DIDN'T HAVE Ns TO BEGIN WITH.
if __name__ == "__main__":
    from pathlib import Path
    import tempfile

    test_mask_nuc_range_from_sam(Path("./tests/"))
    print("Test completed successfully.")