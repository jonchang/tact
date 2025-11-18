#!/usr/bin/env python3
"""
Benchmarking script for TACT performance across different Python implementations.

This script measures the execution time of TACT commands using uv to manage Python versions:
- uv (pypy) - PyPy 3.11
- uv (python3.14) - Python 3.14 without JIT
- uv (python3.14-jit) - Python 3.14 with experimental JIT (PYTHON_JIT=1)

uv automatically downloads and manages Python versions, ensuring Python 3.14 has JIT enabled. The JIT version sets PYTHON_JIT=1
environment variable and verifies JIT is enabled using sys._jit.is_enabled().

Usage:
    python benchmark.py [--dataset carangaria|percomorphaceae] [--runs N] [--warmup N]
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from statistics import mean, stdev
from typing import Any, Dict, List, Optional, Tuple


class BenchmarkRunner:
    """Runs benchmarks across different Python implementations."""

    def __init__(
        self,
        dataset: str = "carangaria",
        runs: int = 5,
        warmup: int = 1,
        examples_dir: Optional[Path] = None,
    ):
        """Initialize the benchmark runner.

        Args:
            dataset: Dataset to use ('carangaria' or 'percomorphaceae')
            runs: Number of benchmark runs to perform
            warmup: Number of warmup runs before benchmarking
            examples_dir: Path to examples directory (defaults to ./examples)
        """
        self.dataset = dataset.lower()
        self.runs = runs
        self.warmup = warmup
        self.examples_dir = examples_dir or Path(__file__).parent / "examples"
        self.temp_dir = None

        # Python implementations to test (using uv to manage Python versions)
        # uv will download and manage these Python versions automatically
        self.implementations = {
            "uv (pypy)": {
                "python_version": "pypy3.11",
                "install_cmd": ["uv", "python", "install", "pypy3.11"],
                "run_cmd": ["uv", "run", "--python", "pypy3.11", "--with", "tact"],
                "env": None,  # No special environment variables
            },
            "uv (python3.14)": {
                "python_version": "3.14",
                "install_cmd": ["uv", "python", "install", "3.14"],
                "run_cmd": ["uv", "run", "--python", "3.14", "--with", "tact"],
                "env": None,  # No JIT enabled
            },
            "uv (python3.14-jit)": {
                "python_version": "3.14",
                "install_cmd": ["uv", "python", "install", "3.14"],
                "run_cmd": ["uv", "run", "--python", "3.14", "--with", "tact"],
                "env": {"PYTHON_JIT": "1"},  # Enable JIT
            },
        }

    def setup_dataset(self) -> Tuple[Path, Path, Path]:
        """Set up the dataset files and return paths to backbone, taxonomy, and output.

        Returns:
            Tuple of (backbone_path, taxonomy_path, output_base)
        """
        # Create temporary directory for outputs
        self.temp_dir = Path(tempfile.mkdtemp(prefix="tact_benchmark_"))
        output_base = self.temp_dir / "result"

        if self.dataset == "carangaria":
            csv_file = self.examples_dir / "Carangaria.csv"
            tre_file = self.examples_dir / "Carangaria.tre"
            taxonomy_file = self.temp_dir / "Carangaria.taxonomy.tre"
        elif self.dataset == "percomorphaceae":
            # For Percomorphaceae, look for files in examples directory
            # or allow user to provide via --examples-dir
            csv_file = self.examples_dir / "Percomorphaceae.csv"
            tre_file = self.examples_dir / "Percomorphaceae.tre"
            taxonomy_file = self.temp_dir / "Percomorphaceae.taxonomy.tre"
            if not csv_file.exists() or not tre_file.exists():
                raise FileNotFoundError(
                    f"Percomorphaceae dataset files not found in {self.examples_dir}. "
                    "Please provide Percomorphaceae.csv and Percomorphaceae.tre files."
                )
        else:
            raise ValueError(f"Unknown dataset: {self.dataset}")

        if not csv_file.exists():
            raise FileNotFoundError(f"CSV file not found: {csv_file}")
        if not tre_file.exists():
            raise FileNotFoundError(f"Tree file not found: {tre_file}")

        return tre_file, csv_file, taxonomy_file, output_base

    def build_taxonomy(
        self, csv_file: Path, taxonomy_file: Path, run_cmd: List[str], env: Optional[Dict[str, str]] = None
    ) -> float:
        """Build the taxonomy tree from CSV.

        Args:
            csv_file: Path to CSV file
            taxonomy_file: Path to output taxonomy file
            run_cmd: Command to run (e.g., ['uv', 'run', '--python', '3.11', '--with', 'tact'])
            env: Optional environment variables to set

        Returns:
            Execution time in seconds
        """
        cmd = run_cmd + [
            "tact_build_taxonomic_tree",
            str(csv_file),
            "--output",
            str(taxonomy_file),
        ]

        # Merge with current environment
        process_env = os.environ.copy()
        if env:
            process_env.update(env)

        start = time.perf_counter()
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
            env=process_env,
        )
        elapsed = time.perf_counter() - start

        if not taxonomy_file.exists():
            raise RuntimeError(f"Taxonomy file not created: {taxonomy_file}\n{result.stderr}")

        return elapsed

    def run_tact_add_taxa(
        self,
        backbone: Path,
        taxonomy: Path,
        output_base: Path,
        run_cmd: List[str],
        env: Optional[Dict[str, str]] = None,
    ) -> float:
        """Run tact_add_taxa command.

        Args:
            backbone: Path to backbone tree
            taxonomy: Path to taxonomy tree
            output_base: Base name for output files
            run_cmd: Command to run (e.g., ['uv', 'run', '--python', '3.11', '--with', 'tact'])
            env: Optional environment variables to set

        Returns:
            Execution time in seconds
        """
        cmd = run_cmd + [
            "tact_add_taxa",
            "--backbone",
            str(backbone),
            "--taxonomy",
            str(taxonomy),
            "--output",
            str(output_base),
        ]

        # Merge with current environment
        process_env = os.environ.copy()
        if env:
            process_env.update(env)

        start = time.perf_counter()
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
            env=process_env,
        )
        elapsed = time.perf_counter() - start

        return elapsed

    def ensure_python_version(self, name: str, impl_config: Dict[str, Any]) -> bool:
        """Ensure a Python version is installed via uv and verify it works.

        Args:
            name: Name of the implementation
            impl_config: Implementation configuration dict with install_cmd and run_cmd

        Returns:
            True if available or successfully installed and verified, False otherwise
        """
        # Check if uv is available
        try:
            subprocess.run(
                ["uv", "--version"],
                capture_output=True,
                text=True,
                check=True,
            )
        except (subprocess.CalledProcessError, FileNotFoundError):
            print(f"Error: uv is not installed. Please install uv first:")
            print("  curl -LsSf https://astral.sh/uv/install.sh | sh")
            return False

        # Install Python version if needed
        print(f"Ensuring Python version is available: {impl_config['python_version']}")
        try:
            result = subprocess.run(
                impl_config["install_cmd"],
                capture_output=True,
                text=True,
                check=True,
            )
        except subprocess.CalledProcessError as e:
            print(f"Failed to install Python version {impl_config['python_version']}: {e.stderr}")
            return False

        # Verify the interpreter can run by checking if tact is available
        print(f"Verifying {name} can run tact commands...", end=" ", flush=True)
        try:
            # Test with --help to verify the command works
            verify_cmd = impl_config["run_cmd"] + ["tact_add_taxa", "--help"]
            process_env = os.environ.copy()
            if impl_config.get("env"):
                process_env.update(impl_config["env"])
            result = subprocess.run(
                verify_cmd,
                capture_output=True,
                text=True,
                check=True,
                env=process_env,
            )
            print("✓")
        except subprocess.CalledProcessError as e:
            print(f"✗ Failed: {e.stderr}")
            return False

        # For Python 3.14 with JIT, verify JIT is enabled
        if name == "uv (python3.14-jit)":
            print(f"Verifying Python 3.14 JIT is enabled...", end=" ", flush=True)
            try:
                code = "import sys; jit = getattr(sys, '_jit', None); print('True' if jit and jit.is_enabled() else 'False')"
                process_env = os.environ.copy()
                process_env["PYTHON_JIT"] = "1"
                result = subprocess.run(
                    ["uv", "run", "--python", "3.14", "python", "-c", code],
                    capture_output=True,
                    text=True,
                    check=True,
                    env=process_env,
                )
                if result.stdout.strip() != "True":
                    print(f"✗ JIT not enabled. Output: {result.stdout.strip()}")
                    return False
                print("✓")
            except subprocess.CalledProcessError as e:
                print(f"✗ Failed to verify JIT: {e.stderr}")
                return False

        return True

    def benchmark_implementation(
        self,
        name: str,
        run_cmd: List[str],
        backbone: Path,
        taxonomy: Path,
        output_base: Path,
        env: Optional[Dict[str, str]] = None,
    ) -> Dict[str, float]:
        """Benchmark a single Python implementation.

        Args:
            name: Name of the implementation
            run_cmd: Command to run (e.g., ['uv', 'run', '--python', '3.11', '--with', 'tact'])
            backbone: Path to backbone tree
            taxonomy: Path to taxonomy tree
            output_base: Base name for output files
            env: Optional environment variables to set

        Returns:
            Dictionary with timing statistics
        """
        print(f"\n{'='*60}")
        print(f"Benchmarking: {name}")
        print(f"{'='*60}")

        times = []

        # Warmup runs
        if self.warmup > 0:
            print(f"Running {self.warmup} warmup run(s)...")
            for i in range(self.warmup):
                print(f"  Warmup {i+1}/{self.warmup}...", end=" ", flush=True)
                try:
                    self.run_tact_add_taxa(backbone, taxonomy, output_base, run_cmd, env)
                    print("✓")
                except subprocess.CalledProcessError as e:
                    print(f"✗ Failed: {e}")
                    return {"error": str(e)}

        # Actual benchmark runs
        print(f"\nRunning {self.runs} benchmark run(s)...")
        for i in range(self.runs):
            print(f"  Run {i+1}/{self.runs}...", end=" ", flush=True)
            try:
                elapsed = self.run_tact_add_taxa(backbone, taxonomy, output_base, run_cmd, env)
                times.append(elapsed)
                print(f"✓ ({elapsed:.2f}s)")
            except subprocess.CalledProcessError as e:
                print(f"✗ Failed: {e.stderr}")
                return {"error": str(e)}

        if not times:
            return {"error": "No successful runs"}

        # Calculate statistics
        mean_time = mean(times)
        min_time = min(times)
        max_time = max(times)
        std_time = stdev(times) if len(times) > 1 else 0.0

        return {
            "mean": mean_time,
            "min": min_time,
            "max": max_time,
            "std": std_time,
            "runs": times,
        }

    def run(self) -> Dict[str, Dict[str, float]]:
        """Run all benchmarks.

        Returns:
            Dictionary mapping implementation names to their timing statistics
        """
        print(f"Benchmarking TACT with dataset: {self.dataset}")
        print(f"Runs: {self.runs}, Warmup: {self.warmup}")

        # Set up dataset
        backbone, csv_file, taxonomy_file, output_base = self.setup_dataset()

        results = {}

        # Ensure all Python versions are installed
        print("\nEnsuring Python versions are available...")
        for name, impl_config in self.implementations.items():
            if not self.ensure_python_version(name, impl_config):
                print(f"\nSkipping {name} (failed to install/verify)")
                results[name] = {"error": "Implementation not available"}
                continue

        # Build taxonomy once (we'll reuse it for all implementations)
        print("\nBuilding taxonomy tree...")
        # Use python3.14 (without JIT) for building taxonomy (should be fast and available)
        try:
            build_cmd = self.implementations["uv (python3.14)"]["run_cmd"]
            build_env = self.implementations["uv (python3.14)"]["env"]
            build_time = self.build_taxonomy(csv_file, taxonomy_file, build_cmd, build_env)
            print(f"Taxonomy built in {build_time:.2f}s")
        except subprocess.CalledProcessError as e:
            print(f"Failed to build taxonomy: {e.stderr}")
            return {"error": "Failed to build taxonomy"}

        # Benchmark each implementation
        for name, impl_config in self.implementations.items():
            if name in results and "error" in results[name]:
                # Already marked as unavailable
                continue

            # Create a unique output base for each implementation
            safe_name = name.replace(" ", "_").replace("(", "").replace(")", "").replace(".", "_")
            impl_output_base = output_base.parent / f"{output_base.name}_{safe_name}"

            result = self.benchmark_implementation(
                name,
                impl_config["run_cmd"],
                backbone,
                taxonomy_file,
                impl_output_base,
                impl_config["env"],
            )
            results[name] = result

        return results

    def cleanup(self):
        """Clean up temporary files."""
        if self.temp_dir and self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)
            print(f"\nCleaned up temporary directory: {self.temp_dir}")

    def format_results(self, results: Dict[str, Dict[str, float]]) -> str:
        """Format results as a markdown table.

        Args:
            results: Dictionary of benchmark results

        Returns:
            Formatted markdown table string
        """
        # Find baseline (pypy) for relative calculations
        baseline = None
        if "uv (pypy)" in results and "mean" in results["uv (pypy)"]:
            baseline = results["uv (pypy)"]["mean"]
        elif "uv (python3.14)" in results and "mean" in results["uv (python3.14)"]:
            baseline = results["uv (python3.14)"]["mean"]

        lines = []
        lines.append("| python | Mean [s] | Min [s] | Max [s] | Relative |")
        lines.append("|:---|---:|---:|---:|---:|")

        # Order implementations for consistent output
        order = ["uv (pypy)", "uv (python3.14)", "uv (python3.14-jit)"]
        for name in order:
            if name not in results:
                continue

            result = results[name]
            if "error" in result:
                error_msg = result["error"][:50]  # Truncate long error messages
                lines.append(f"| `{name}` | Error: {error_msg} | | | |")
                continue

            mean_time = result["mean"]
            min_time = result["min"]
            max_time = result["max"]
            std_time = result["std"]

            # Calculate relative performance
            if baseline:
                relative = mean_time / baseline
                relative_std = std_time / baseline if std_time > 0 else 0
                relative_str = f"{relative:.2f} ± {relative_std:.2f}" if relative_std > 0 else f"{relative:.2f}"
            else:
                relative_str = "—"

            # Format mean with std dev (matching hyperfine format)
            mean_str = f"{mean_time:.3f} ± {std_time:.3f}" if std_time > 0 else f"{mean_time:.3f}"

            lines.append(
                f"| `{name}` | {mean_str} | {min_time:.3f} | {max_time:.3f} | {relative_str} |"
            )

        return "\n".join(lines)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.cleanup()


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Benchmark TACT performance across different Python implementations"
    )
    parser.add_argument(
        "--dataset",
        choices=["carangaria", "percomorphaceae"],
        default="carangaria",
        help="Dataset to use for benchmarking (default: carangaria)",
    )
    parser.add_argument(
        "--runs",
        type=int,
        default=5,
        help="Number of benchmark runs (default: 5)",
    )
    parser.add_argument(
        "--warmup",
        type=int,
        default=1,
        help="Number of warmup runs (default: 1)",
    )
    parser.add_argument(
        "--examples-dir",
        type=Path,
        help="Path to examples directory (default: ./examples)",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Output results as JSON instead of markdown table",
    )

    args = parser.parse_args()

    with BenchmarkRunner(
        dataset=args.dataset,
        runs=args.runs,
        warmup=args.warmup,
        examples_dir=args.examples_dir,
    ) as runner:
        results = runner.run()

        if args.json:
            print("\n" + json.dumps(results, indent=2))
        else:
            print("\n" + "="*60)
            print("RESULTS")
            print("="*60)
            print(runner.format_results(results))


if __name__ == "__main__":
    main()
