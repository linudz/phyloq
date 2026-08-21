#!/usr/bin/env python3

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


CAASTOOLS_DIR = Path(__file__).resolve().parents[1]
CT = CAASTOOLS_DIR / "ct"

sys.path.insert(0, str(CAASTOOLS_DIR))

from modules.init_bootstrap import pool_discovery_hypotheses


FG_SPECIES = [
    "Aotus_trivirgatus",
    "Ateles_geoffroyi",
    "Cercocebus_torquatus",
    "Eulemur_fulvus",
    "Eulemur_macaco",
    "Lemur_catta",
    "Nomascus_concolor",
]
BG_SPECIES = [
    "Alouatta_seniculus",
    "Eulemur_rubriventer",
    "Macaca_thibetana",
    "Macaca_tonkeana",
    "Mico_humeralifer",
    "Trachypithecus_pileatus",
]


class PooledDiscoveryTests(unittest.TestCase):

    def write_inputs(self, directory):
        directory = Path(directory)
        pool_file = directory / "full-pools.cfg"
        alignment = directory / "pooled.phy"

        pool_file.write_text(
            "".join(species + "\t1\n" for species in FG_SPECIES)
            + "".join(species + "\t0\n" for species in BG_SPECIES)
        )

        amino_acids = {
            **{species: ("A" if index < 4 else "G") for index, species in enumerate(FG_SPECIES)},
            **{species: ("V" if index < 4 else "I") for index, species in enumerate(BG_SPECIES)},
        }
        alignment.write_text(
            "13 1\n"
            + "".join(species + "  " + amino_acids[species] + "\n" for species in FG_SPECIES + BG_SPECIES)
        )
        return pool_file, alignment

    def test_max_generates_all_525_unique_four_by_four_comparisons(self):
        with tempfile.TemporaryDirectory() as directory:
            pool_file, _ = self.write_inputs(directory)
            pooled = pool_discovery_hypotheses(pool_file, comparisons="max", seed=260811)

            self.assertEqual(pooled.cycles, 525)
            self.assertEqual(pooled.pool_maximum_comparisons, 525)
            self.assertEqual(pooled.discovery_fg, sorted(FG_SPECIES))
            self.assertEqual(pooled.discovery_bg, sorted(BG_SPECIES))
            comparisons = {
                (tuple(pooled.trait2fg[hypothesis]), tuple(pooled.trait2bg[hypothesis]))
                for hypothesis in pooled.trait2fg
            }
            self.assertEqual(len(comparisons), 525)

    def test_fixed_count_is_reproducible_and_seed_dependent(self):
        with tempfile.TemporaryDirectory() as directory:
            pool_file, _ = self.write_inputs(directory)
            first = pool_discovery_hypotheses(pool_file, comparisons=100, seed=260811)
            repeated = pool_discovery_hypotheses(pool_file, comparisons="100", seed=260811)
            different = pool_discovery_hypotheses(pool_file, comparisons=100, seed=260812)

            first_pairs = list(zip(first.trait2fg.values(), first.trait2bg.values()))
            repeated_pairs = list(zip(repeated.trait2fg.values(), repeated.trait2bg.values()))
            different_pairs = list(zip(different.trait2fg.values(), different.trait2bg.values()))
            self.assertEqual(first_pairs, repeated_pairs)
            self.assertNotEqual(first_pairs, different_pairs)
            self.assertEqual(len(set((tuple(fg), tuple(bg)) for fg, bg in first_pairs)), 100)

    def test_invalid_pooling_requests_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            pool_file, _ = self.write_inputs(directory)
            cases = [
                ({"comparisons": 526}, "only 525 unique"),
                ({"comparisons": 0}, "at least 1"),
                ({"comparisons": "invalid"}, "must be 'max'"),
                ({"fg_size": 8}, "exceeds the complete FG pool"),
                ({"bg_size": 7}, "exceeds the complete BG pool"),
            ]
            for arguments, message in cases:
                with self.subTest(arguments=arguments):
                    with self.assertRaisesRegex(ValueError, message):
                        pool_discovery_hypotheses(pool_file, **arguments)

    def test_cli_builds_and_saves_100_hypotheses_then_emits_events(self):
        with tempfile.TemporaryDirectory() as directory:
            pool_file, alignment = self.write_inputs(directory)
            output = Path(directory) / "pooled.caas.tsv"
            event_output = Path(directory) / "pooled.caas.events.tsv"
            hypotheses_output = Path(directory) / "pooled.caas.hypotheses.tsv"

            command = [
                sys.executable,
                str(CT),
                "pooled-discovery",
                "-a", str(alignment),
                "-t", str(pool_file),
                "-o", str(output),
                "--event-output", str(event_output),
                "--fmt", "phylip-relaxed",
                "--fg-size", "4",
                "--bg-size", "4",
                "--comparisons", "100",
                "--seed", "260811",
                "--min-fg-observed", "3",
                "--min-bg-observed", "3",
            ]
            run = subprocess.run(
                command,
                cwd=CAASTOOLS_DIR,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )

            self.assertEqual(run.returncode, 0, run.stdout)
            hypothesis_lines = hypotheses_output.read_text().splitlines()
            self.assertEqual(len(hypothesis_lines), 100)
            self.assertEqual(hypothesis_lines[0].split("\t")[0], "b_1")
            self.assertEqual(hypothesis_lines[-1].split("\t")[0], "b_100")
            self.assertEqual(len({tuple(line.split("\t")[1:]) for line in hypothesis_lines}), 100)

            legacy_fields = output.read_text().strip().split("\t")
            self.assertEqual(legacy_fields[2], "100")
            with event_output.open() as handle:
                event = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(event["fg_discovery_count"], "7")
            self.assertEqual(event["bg_discovery_count"], "6")
            self.assertEqual(event["event_pattern"], "4")
            self.assertIn("Possible unique comparisons: 525", run.stdout)
            self.assertIn("Selected comparisons: 100", run.stdout)
            self.assertIn("Hypotheses output: " + str(hypotheses_output), run.stdout)

            reused_output = Path(directory) / "pooled.reused.caas.tsv"
            reused_event_output = Path(directory) / "pooled.reused.caas.events.tsv"
            reused_command = list(command)
            reused_command[reused_command.index(str(output))] = str(reused_output)
            reused_command[reused_command.index(str(event_output))] = str(reused_event_output)
            reused_command.extend([
                "-s", str(hypotheses_output),
                "--hypotheses-output", "none",
            ])
            reused_run = subprocess.run(
                reused_command,
                cwd=CAASTOOLS_DIR,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )

            self.assertEqual(reused_run.returncode, 0, reused_run.stdout)
            self.assertEqual(output.read_text(), reused_output.read_text())
            self.assertEqual(event_output.read_text(), reused_event_output.read_text())
            self.assertIn("Saved hypothesis reuse", reused_run.stdout)


if __name__ == "__main__":
    unittest.main()
