#!/usr/bin/env python3

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


CAASTOOLS_DIR = Path(__file__).resolve().parents[1]
CT = CAASTOOLS_DIR / "ct"


class MinimumObservedSpeciesTests(unittest.TestCase):

    def write_bootstrap_inputs(self, directory):
        directory = Path(directory)
        alignment = directory / "coverage.phy"
        trait_config = directory / "trait.cfg"
        resampled_traits = directory / "resampled.tsv"

        # Position zero has only one observed FG and one observed BG species.
        # Position one has three observed species in each group. Additional
        # alignment species keep the overall gap fraction below the standard
        # 0.5 filter without contributing to the resampled groups.
        alignment.write_text(
            "16 2\n"
            "sp1   AA\n"
            "sp2   -A\n"
            "sp3   -A\n"
            "sp4   --\n"
            "sp5   BB\n"
            "sp6   -B\n"
            "sp7   -B\n"
            "sp8   --\n"
            "sp9   CC\n"
            "sp10  CC\n"
            "sp11  CC\n"
            "sp12  CC\n"
            "sp13  DD\n"
            "sp14  DD\n"
            "sp15  DD\n"
            "sp16  DD\n"
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

    @staticmethod
    def counts_by_position(output):
        return {
            line.split("\t")[0]: int(line.split("\t")[1])
            for line in output.read_text().splitlines()
        }

    def test_default_one_preserves_low_coverage_cycle(self):
        with tempfile.TemporaryDirectory() as directory:
            alignment, trait_config, resampled_traits = self.write_bootstrap_inputs(directory)
            output = Path(directory) / "default.bootstrap.tsv"

            run = self.run_bootstrap(alignment, trait_config, resampled_traits, output)

            self.assertEqual(run.returncode, 0, run.stdout)
            self.assertEqual(
                self.counts_by_position(output),
                {"coverage@0": 1, "coverage@1": 1},
            )

    def test_three_per_group_excludes_only_low_coverage_cycle(self):
        with tempfile.TemporaryDirectory() as directory:
            alignment, trait_config, resampled_traits = self.write_bootstrap_inputs(directory)
            output = Path(directory) / "minimum-three.bootstrap.tsv"

            run = self.run_bootstrap(
                alignment,
                trait_config,
                resampled_traits,
                output,
                "--min-fg-observed", "3",
                "--min-bg-observed", "3",
            )

            self.assertEqual(run.returncode, 0, run.stdout)
            self.assertEqual(
                self.counts_by_position(output),
                {"coverage@0": 0, "coverage@1": 1},
            )
            self.assertIn("Foreground: 3 of 4", run.stdout)
            self.assertIn("Background: 3 of 4", run.stdout)

    def test_invalid_minimum_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            alignment, trait_config, resampled_traits = self.write_bootstrap_inputs(directory)

            for argument, value in (("--min_fg_observed", "0"), ("--min_bg_observed", "5")):
                with self.subTest(argument=argument, value=value):
                    output = Path(directory) / (argument.lstrip("-") + ".tsv")
                    run = self.run_bootstrap(
                        alignment,
                        trait_config,
                        resampled_traits,
                        output,
                        argument, value,
                    )

                    self.assertNotEqual(run.returncode, 0)
                    self.assertIn(argument, run.stdout)


if __name__ == "__main__":
    unittest.main()
