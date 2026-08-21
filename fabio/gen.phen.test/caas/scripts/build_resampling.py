#!/usr/bin/env python3

"""Build one CAAStools resampling table from binary comparison configs."""

import argparse
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
DEFAULT_CONFIG_DIR = PROJECT_DIR.parent / "lq.table2.4species.comparisons" / "configurations" / "caas"
DEFAULT_OUTPUT = PROJECT_DIR / "inputs" / "longevity.100-comparisons.resampling.tsv"
DEFAULT_TEMPLATE = PROJECT_DIR / "inputs" / "longevity.template.001.caas.cfg"


def parse_config(config_path):
    foreground = []
    background = []
    observed_species = set()

    for line_number, line in enumerate(config_path.read_text().splitlines(), start=1):
        if not line.strip():
            continue

        fields = line.split("\t")
        if len(fields) != 2 or fields[1] not in {"0", "1"}:
            raise ValueError(
                str(config_path) + ":" + str(line_number) +
                " must contain a species and a binary state separated by one tab"
            )

        species, state = fields
        if species in observed_species:
            raise ValueError(str(config_path) + " contains duplicate species " + species)
        observed_species.add(species)

        if state == "1":
            foreground.append(species)
        else:
            background.append(species)

    if len(foreground) != 4 or len(background) != 4:
        raise ValueError(
            str(config_path) + " must contain exactly 4 FG and 4 BG species; found " +
            str(len(foreground)) + "/" + str(len(background))
        )

    return foreground, background


def build_resampling(config_dir, output_path, template_path, expected_cycles):
    config_paths = sorted(config_dir.glob("*.caas.cfg"))

    if len(config_paths) != expected_cycles:
        raise ValueError(
            "expected " + str(expected_cycles) + " comparison configs in " +
            str(config_dir) + ", found " + str(len(config_paths))
        )

    output_lines = []
    observed_comparisons = set()

    for config_path in config_paths:
        comparison_id = config_path.name.replace(".caas.cfg", "")
        if not comparison_id.isdigit():
            raise ValueError("comparison config names must begin with a numeric ID: " + str(config_path))

        foreground, background = parse_config(config_path)
        comparison = (tuple(sorted(foreground)), tuple(sorted(background)))

        if comparison in observed_comparisons:
            raise ValueError("duplicate FG/BG comparison found in " + str(config_path))
        observed_comparisons.add(comparison)

        cycle_id = "b_" + str(int(comparison_id))
        output_lines.append(
            "\t".join([cycle_id, ",".join(foreground), ",".join(background)])
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    template_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(output_lines) + "\n")
    template_path.write_text(config_paths[0].read_text())

    return len(output_lines)


def main():
    parser = argparse.ArgumentParser(
        description="Convert CAAStools binary comparison configs into one resampling table."
    )
    parser.add_argument("--config-dir", type=Path, default=DEFAULT_CONFIG_DIR)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--template-output", type=Path, default=DEFAULT_TEMPLATE)
    parser.add_argument("--expected-cycles", type=int, default=100)
    arguments = parser.parse_args()

    cycles = build_resampling(
        config_dir=arguments.config_dir,
        output_path=arguments.output,
        template_path=arguments.template_output,
        expected_cycles=arguments.expected_cycles,
    )

    print("CAAStools resampling table:", arguments.output)
    print("Bootstrap template config:", arguments.template_output)
    print("Cycles:", cycles)


if __name__ == "__main__":
    main()
