"""RNA structure analysis functionality."""

import os
import tempfile

import polars as pl


@pl.api.register_expr_namespace("rna")
class RNAStructureExpr:
    """Polars expression namespace for RNA structure prediction."""
    
    def __init__(self, expr: pl.Expr):
        self._expr = expr

    def predict_structure(self) -> pl.Expr:
        """Predict RNA structure and minimum free energy"""
        import RNA

        def _predict(seq: str) -> dict:
            if len(seq) > 10000:
                return {"structure": None, "mfe": None}
            try:
                ss_string, mfe = RNA.fold_compound(seq).mfe()
                return {"structure": ss_string, "mfe": mfe}
            except Exception:
                return {"structure": None, "mfe": None}

        return self._expr.map_elements(
            _predict,
            return_dtype=pl.Struct({"structure": pl.String, "mfe": pl.Float64}),
        )

    def predict_structure_with_tool(self, tool: str = "ViennaRNA") -> pl.Expr:
        """Predict RNA structure using specified tool"""

        def _predict_vienna(seq: str) -> dict:
            if len(seq) > 10000:
                return {"structure": None, "mfe": None}
            try:
                import RNA

                ss_string, mfe = RNA.fold_compound(seq).mfe()
                return {"structure": ss_string, "mfe": mfe}
            except Exception:
                return {"structure": None, "mfe": None}

        def _predict_linearfold(seq: str) -> dict:
            import subprocess

            if len(seq) > 10000:
                return {"structure": None, "mfe": None}
            try:
                # Convert T to U for RNA folding
                seq = seq.replace("T", "U")
                with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_in:
                    # Write sequence directly without FASTA header
                    temp_in.write(seq)
                    temp_in.flush()
                    temp_in_name = temp_in.name

                # Create temp file for output
                with tempfile.NamedTemporaryFile(mode="r+", delete=False) as temp_out:
                    temp_out_name = temp_out.name

                # Run LinearFold with default beam size (100)
                subprocess.run(
                    ["linearfold", temp_in_name],
                    stdout=open(temp_out_name, "w"),
                    stderr=subprocess.PIPE,
                )

                # Parse the results
                with open(temp_out_name, "r") as f:
                    result = f.read().strip().split("\n")
                    if len(result) >= 2:
                        structure, mfe_str = result[1].split()
                        mfe = float(mfe_str.replace("(", "").replace(")", ""))

                        # Clean up temp files
                        os.unlink(temp_in_name)
                        os.unlink(temp_out_name)

                        return {"structure": structure, "mfe": mfe}

                # Clean up temp files
                os.unlink(temp_in_name)
                if os.path.exists(temp_out_name):
                    os.unlink(temp_out_name)

                return {"structure": None, "mfe": None}

            except Exception:
                return {"structure": None, "mfe": None}

        if tool.lower() == "linearfold":
            return self._expr.map_elements(
                _predict_linearfold,
                return_dtype=pl.Struct({"structure": pl.String, "mfe": pl.Float64}),
            )
        else:  # Default to ViennaRNA
            return self._expr.map_elements(
                _predict_vienna,
                return_dtype=pl.Struct({"structure": pl.String, "mfe": pl.Float64}),
            ) 