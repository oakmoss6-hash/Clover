"""Clover command-line and default configuration."""

from __future__ import annotations

import getopt
import json
import sys


SHORT_OPTIONS = "-I:-L:-D:-V:-H:-T:-P:-O:-h"
LONG_OPTIONS = [
    "help",
    "low",
    "align",
    "no-fast",
    "no-tag",
    "stat",
    "export-cluster-cores=",
    "input=",
    "reconstruct",
    "reconstruct-backend=",
    "reconstruct-backbone=",
    "consensus-output=",
]

RECONSTRUCTION_BACKENDS = ("wfa", "edlib", "nw")
RECONSTRUCTION_BACKBONES = ("core", "support_length", "max_span")


config_dict = {
    "read_len": 152,
    "end_tree_len": 15,
    "other_tree_len": 15,
    "other_tree_nums": 2,
    "thd_tree_loc": 40,
    "four_tree_loc": 40,
    "Vertical_drift": 2,
    "Horizontal_drift": 3,
    "tree_threshold": 10,
    "now_clust_threshold": 8,
    "tag_nums": 1,
    "processes_nums": 0,
    "Cluster_size_threshold": 1,
    "h_index_nums": 0,
    "e_index_nums": 0,
    "read_len_min": 0,
    "align_fuc": False,
    "mmr_mode": False,
    "Virtual_mode": True,
    "fast_mode": True,
    "tag_mode": False,
    "Statistical_model": False,
    "same_tree_len": True,
    "now_align_alg": False,
    "core_export_path": None,
    # Optional multi-read reconstruction.
    "reconstruct": False,
    "reconstruct_backend": "wfa",
    "reconstruct_backbone": "core",
    "consensus_output_path": None,
}


def load_json(path):
    """Load a JSON config while ignoring full-line // comments."""
    lines = []
    with open(path) as source:
        for row in source:
            if not row.strip().startswith("//"):
                lines.append(row)
    return json.loads("\n".join(lines))


def generate_vertical_drifts_list(x):
    return list(range(-x, x + 1))


def _parse_options(argv):
    return getopt.getopt(argv, SHORT_OPTIONS, LONG_OPTIONS)


def out_put_config():
    """Apply command-line overrides to Clover's historical config dict."""
    options, _ = _parse_options(sys.argv[1:])

    for opt_name, opt_value in options:
        if opt_name in {"-h", "--help"}:
            print("Please see README.md")
        elif opt_name in {"-I", "--input"}:
            config_dict["input_path"] = opt_value
        elif opt_name == "-L":
            config_dict["read_len"] = int(opt_value)
        elif opt_name == "-D":
            config_dict["end_tree_len"] = int(opt_value)
            if config_dict["same_tree_len"]:
                config_dict["other_tree_len"] = int(opt_value)
        elif opt_name == "-V":
            config_dict["Vertical_drift"] = generate_vertical_drifts_list(
                int(opt_value)
            )
        elif opt_name == "-H":
            config_dict["Horizontal_drift"] = int(opt_value)
        elif opt_name == "-T":
            config_dict["tag_nums"] = int(opt_value)
            config_dict["tag_mode"] = True
        elif opt_name == "-P":
            config_dict["processes_nums"] = int(opt_value)
        elif opt_name == "-O":
            config_dict["output_file"] = opt_value + ".txt"
        elif opt_name == "--align":
            config_dict["align_fuc"] = True
        elif opt_name == "--no-fast":
            config_dict["fast_mode"] = False
        elif opt_name == "--no-tag":
            config_dict["Virtual_mode"] = False
        elif opt_name == "--stat":
            config_dict["Statistical_model"] = True
        elif opt_name == "--export-cluster-cores":
            if not opt_value:
                raise ValueError(
                    "--export-cluster-cores requires a non-empty path"
                )
            config_dict["core_export_path"] = opt_value
        elif opt_name == "--reconstruct":
            config_dict["reconstruct"] = True
        elif opt_name == "--reconstruct-backend":
            backend = opt_value.lower()
            if backend not in RECONSTRUCTION_BACKENDS:
                raise ValueError(
                    "--reconstruct-backend must be one of: "
                    + ", ".join(RECONSTRUCTION_BACKENDS)
                )
            config_dict["reconstruct_backend"] = backend
        elif opt_name == "--reconstruct-backbone":
            backbone = opt_value.lower()
            if backbone not in RECONSTRUCTION_BACKBONES:
                raise ValueError(
                    "--reconstruct-backbone must be one of: "
                    + ", ".join(RECONSTRUCTION_BACKBONES)
                )
            config_dict["reconstruct_backbone"] = backbone
        elif opt_name == "--consensus-output":
            if not opt_value:
                raise ValueError("--consensus-output requires a path")
            config_dict["consensus_output_path"] = opt_value
        elif opt_name == "--low":
            config_dict["mmr_mode"] = True
            config_dict["fast_mode"] = False
            config_dict["align_fuc"] = False
            config_dict["Statistical_model"] = False
            config_dict["Virtual_mode"] = False

    if config_dict.get("reconstruct"):
        if config_dict.get("align_fuc"):
            raise ValueError(
                "old --align and new --reconstruct are mutually exclusive"
            )
        # Reconstruction is intentionally streaming: raw reads are not
        # preloaded into the parent's data_dict.
        config_dict["fast_mode"] = False
        if not config_dict.get("consensus_output_path"):
            config_dict["consensus_output_path"] = "clover_consensus.tsv"

    if config_dict["read_len_min"] == 0:
        config_dict["read_len_min"] = config_dict["read_len"] - 5

    # Preserve Clover's historical effective drift behavior.
    if isinstance(config_dict["Vertical_drift"], int):
        config_dict["Vertical_drift"] = generate_vertical_drifts_list(
            config_dict["Horizontal_drift"]
        )

    config_dict["tag"] = r"""
      ___           ___       ___           ___           ___           ___
     /\  \         /\__\     /\  \         /\__\         /\  \         /\  \
    /::\  \       /:/  /    /::\  \       /:/  /        /::\  \       /::\  \
   /:/\:\  \     /:/  /    /:/\:\  \     /:/  /        /:/\:\  \     /:/\:\  \
  /:/  \:\  \   /:/  /    /:/  \:\  \   /:/__/  ___   /::\~\:\  \   /::\~\:\  \
 /:/__/ \:\__\ /:/__/    /:/__/ \:\__\  |:|  | /\__\ /:/\:\ \:\__\ /:/\:\ \:\__\
 \:\  \  \/__/ \:\  \    \:\  \ /:/  /  |:|  |/:/  / \:\~\:\ \/__/ \/_|::\/:/  /
  \:\  \        \:\  \    \:\  /:/  /   |:|__/:/  /   \:\ \:\__\      |:|::/  /
   \:\  \        \:\  \    \:\/:/  /     \::::/__/     \:\ \/__/      |:|\/__/
    \:\__\        \:\__\    \::/  /       ~~~~          \:\__\        |:|  |
     \/__/         \/__/     \/__/                       \/__/         \|__|
    """
    return config_dict
