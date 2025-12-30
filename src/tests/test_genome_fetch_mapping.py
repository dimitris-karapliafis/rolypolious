"""
Simple tests for genome_fetch_from_mapping module.
Run with: python src/tests/test_genome_fetch_mapping.py
"""

import os
from pathlib import Path

import polars as pl

from rolypoly.utils.bio.genome_fetch import (
    download_from_ftp_path,
    get_ftp_path_for_taxid,
)

# Get data directory
DATA_DIR = os.environ.get(
    "ROLYPOLY_DATA",
    "<REPO_PATH>/data",
)
MAPPING_PATH = f"{DATA_DIR}/contam/rrna/rrna_to_genome_mapping.parquet"


def run_tests():
    """Run all tests."""
    passed = 0
    failed = 0
    skipped = 0
    
    print("=" * 70)
    print("Testing genome_fetch_from_mapping module")
    print("=" * 70)
    
    # Test 1: Load mapping
    print("\n[1/5] Test loading mapping file...")
    try:
        if not Path(MAPPING_PATH).exists():
            print(f"⚠️  SKIP: Mapping file not found: {MAPPING_PATH}")
            skipped += 1
        else:
            mapping = pl.read_parquet(MAPPING_PATH)
            assert isinstance(mapping, pl.DataFrame)
            assert mapping.height > 0
            assert "query_tax_id" in mapping.columns
            assert "relationship" in mapping.columns
            assert "ftp_path" in mapping.columns
            print(f"✅ PASS: Loaded {mapping.height:,} mappings")
            passed += 1
            
            # Test 2: Get FTP path for known taxid (with relationship info)
            print("\n[2/5] Test getting FTP path for known taxid...")
            # Find a taxid that has self-mapping
            self_mapped = mapping.filter(pl.col("relationship") == "self")
            if self_mapped.height > 0:
                taxid = self_mapped.row(0, named=True)["query_tax_id"]
                ftp_path, relationship, ref_name = get_ftp_path_for_taxid(taxid, mapping)
                if ftp_path:
                    assert "ftp.ncbi.nlm.nih.gov" in ftp_path or relationship is not None
                    assert relationship == "self"
                    print(f"✅ PASS: Found FTP path for taxid {taxid}")
                    print(f"   Path: {ftp_path}")
                    print(f"   Relationship: {relationship}")
                    passed += 1
                else:
                    print(f"⚠️  SKIP: No FTP path found for taxid {taxid}")
                    skipped += 1
            else:
                print(f"⚠️  SKIP: No self-mapped taxids found in mapping")
                skipped += 1
            
            # Test 3: Missing taxid
            print("\n[3/5] Test handling of missing taxid...")
            taxid = 99999999
            ftp_path, relationship, ref_name = get_ftp_path_for_taxid(taxid, mapping)
            assert ftp_path is None and relationship is None
            print(f"✅ PASS: Correctly returned None for missing taxid")
            passed += 1
            
            # Test 4: Mapping coverage
            print("\n[4/5] Test mapping coverage...")
            total = mapping.height
            with_ftp = mapping.filter(
                pl.col("ftp_path").is_not_null()
            ).height
            coverage = with_ftp / total * 100
            assert coverage > 20, f"Coverage too low: {coverage:.1f}%"
            print(f"✅ PASS: Coverage is {with_ftp:,}/{total:,} ({coverage:.1f}%)")
            passed += 1
            
            # Test 5: Relationship distribution
            print("\n[5/5] Test relationship distribution...")
            by_relationship = mapping.group_by("relationship").agg(pl.len().alias("count"))
            print("Relationship distribution:")
            for row in by_relationship.iter_rows(named=True):
                rel = row['relationship'] or 'None'
                count = row['count']
                pct = count / total * 100
                print(f"  {rel:15s}: {count:8,} ({pct:5.1f}%)")
            print(f"✅ PASS: Relationship data is present")
            passed += 1
            
    except Exception as e:
        print(f"❌ FAIL: {e}")
        import traceback
        traceback.print_exc()
        failed += 1
    
    # Summary
    print("\n" + "=" * 70)
    print(f"Results: {passed} passed, {failed} failed, {skipped} skipped")
    print("=" * 70)
    
    return failed == 0


if __name__ == "__main__":
    success = run_tests()
    exit(0 if success else 1)
