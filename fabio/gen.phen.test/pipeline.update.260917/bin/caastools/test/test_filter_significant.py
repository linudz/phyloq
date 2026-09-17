#!/usr/bin/env python3

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


CAASTOOLS_DIR = Path(__file__).resolve().parents[1]
CT = CAASTOOLS_DIR / "ct"

sys.path.insert(0, str(CAASTOOLS_DIR))

from modules.boot import infer_resampled_group_sizes


class FilterSignificantTests(unittest.TestCase):

    def write_bootstrap_inputs(self, directory):
        directory = Path(directory)
        alignment = directory / "toy.phy"
        trait_config = directory / "trait.cfg"
        resampled_traits = directory / "resampled.tsv"

        alignment.write_text(
            "8 2\n"
            "sp1  AA\n"
            "sp2  AA\n"
            "sp3  AA\n"
            "sp4  AA\n"
            "sp5  BB\n"
            "sp6  BB\n"
            "sp7  BC\n"
            "sp8  BD\n"
        )

        trait_config.write_text(
            "sp1\t1\n"
            "sp2\t1\n"
            "sp3\t1\n"
            "sp4\t1\n"
            "sp5\t0\n"
            "sp6\t0\n"
            "sp7\t0\n"
            "sp8\t0\n"
        )

        resampled_traits.write_text(
            "b_1\tsp1,sp2,sp3,sp4\tsp5,sp6,sp7,sp8\n"
            "b_2\tsp1,sp2,sp5,sp6\tsp3,sp4,sp7,sp8\n"
        )

        return alignment, trait_config, resampled_traits

    def run_bootstrap(self, alignment, trait_config, resampled_traits, output, *arguments):
        command = [
            sys.executable,
            str(CT),
            "bootstrap",
            "-a", str(alignment),
            "-t", str(trait_config),
            "-s", str(resampled_traits),
            "-o", str(output),
            "--fmt", "phylip-relaxed",
        ]
        command.extend(arguments)

        return subprocess.run(
            command,
            cwd=CAASTOOLS_DIR,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )

    def test_default_no_preserves_unfiltered_output(self):
        with tempfile.TemporaryDirectory() as directory:
            alignment, trait_config, resampled_traits = self.write_bootstrap_inputs(directory)
            default_output = Path(directory) / "default.bootstrap.tsv"
            explicit_no_output = Path(directory) / "explicit-no.bootstrap.tsv"

            default_run = self.run_bootstrap(
                alignment, trait_config, resampled_traits, default_output
            )
            explicit_no_run = self.run_bootstrap(
                alignment, trait_config, resampled_traits, explicit_no_output,
                "--filter_significant", "no"
            )

            self.assertEqual(default_run.returncode, 0, default_run.stdout)
            self.assertEqual(explicit_no_run.returncode, 0, explicit_no_run.stdout)
            self.assertEqual(default_output.read_text(), explicit_no_output.read_text())
            self.assertEqual(len(default_output.read_text().splitlines()), 2)

    def test_fs_alias_filters_positions_above_threshold(self):
        with tempfile.TemporaryDirectory() as directory:
            alignment, trait_config, resampled_traits = self.write_bootstrap_inputs(directory)
            filtered_output = Path(directory) / "filtered.bootstrap.tsv"

            filtered_run = self.run_bootstrap(
                alignment, trait_config, resampled_traits, filtered_output,
                "-fs", "0.02"
            )

            self.assertEqual(filtered_run.returncode, 0, filtered_run.stdout)
            output_lines = filtered_output.read_text().splitlines()
            self.assertEqual(len(output_lines), 1)
            self.assertTrue(output_lines[0].startswith("toy@1\t"))
            output_fields = output_lines[0].split("\t")
            self.assertEqual(len(output_fields), 7)
            self.assertLessEqual(float(output_fields[5]), 0.02)
            self.assertEqual(str(trait_config), output_fields[6])
            self.assertIn("Resampled group sizes: FG=4, BG=4", filtered_run.stdout)
            self.assertIn("Positions before filtering: 2", filtered_run.stdout)
            self.assertIn("Positions retained: 1", filtered_run.stdout)

    def test_invalid_threshold_returns_an_error(self):
        with tempfile.TemporaryDirectory() as directory:
            alignment, trait_config, resampled_traits = self.write_bootstrap_inputs(directory)

            for threshold in ("1.5", "nan", "not-a-number"):
                with self.subTest(threshold=threshold):
                    output = Path(directory) / (threshold + ".bootstrap.tsv")
                    invalid_run = self.run_bootstrap(
                        alignment, trait_config, resampled_traits, output,
                        "--filter_significant", threshold
                    )

                    self.assertNotEqual(invalid_run.returncode, 0)
                    self.assertIn("--filter_significant must", invalid_run.stdout)

    def test_mixed_group_sizes_are_rejected(self):
        resampled_traits = SimpleNamespace(
            cycles=2,
            trait2fg={"b_1": ["a", "b"], "b_2": ["a", "b", "c"]},
            trait2bg={"b_1": ["c", "d"], "b_2": ["d", "e"]},
        )

        with self.assertRaisesRegex(ValueError, "same FG/BG sizes"):
            infer_resampled_group_sizes(resampled_traits)


if __name__ == "__main__":
    unittest.main()
